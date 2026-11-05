// SwiftDistSeqs: byte-identical, fast re-implementation of mothur 1.48.6 dist.seqs
// for aligned DNA sequences (calc = onegap | eachgap | nogaps, countends, cutoff,
// output = column | lt | square).
//
// Usage (gate style):
//   distseqs --fasta IN --out OUT [--calc C] [--countends T|F] [--cutoff X]
//            [--output column|lt|square] [--processors P]
// Usage (mothur style, writes mothur's file name next to the input):
//   distseqs fasta=IN [calc=C] [countends=T|F] [cutoff=X] [output=O] [processors=P] [outputdir=D]
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <immintrin.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace {

// ------------------------------------------------------------------ options
struct Opts {
    std::string fasta, out, calc = "onegap", output = "column", outputdir;
    std::string countends = "T", cutoff = "1.0";
    int procs = 1;
    bool mothurStyle = false;
};

std::string upper(std::string s) { for (auto& c : s) c = (char)toupper((unsigned char)c); return s; }

// ------------------------------------------------------------------ istream emulation
// Emulates exactly the std::ifstream operations mothur's Sequence(ifstream&) uses.
struct Stream {
    const char* b; size_t n, p = 0; bool eofb = false, failb = false;
    Stream(const char* buf, size_t len) : b(buf), n(len) {}
    bool ok() const { return !failb; }
    int get() {
        if (failb) { return EOF; }
        if (p < n) return (unsigned char)b[p++];
        eofb = failb = true; return EOF;
    }
    int peek() {
        if (failb) return EOF;
        if (p < n) return (unsigned char)b[p];
        eofb = true; return EOF;
    }
    void putback() { eofb = false; if (!failb && p > 0) --p; }
    static bool ws(int c) { return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r'; }
    // operator>>(istream&, string&)
    void token(std::string& s) {
        if (failb) { failb = true; return; }
        while (p < n && ws((unsigned char)b[p])) ++p;
        if (p >= n) { eofb = true; failb = true; return; }  // sentry failed: s untouched
        s.clear();
        while (p < n && !ws((unsigned char)b[p])) s.push_back(b[p++]);
        if (p >= n) eofb = true;
    }
    // std::getline(istream&, string&)
    void getline(std::string& s) {
        if (failb) { failb = true; return; }
        s.clear();
        size_t q = p;
        while (q < n && b[q] != '\n') ++q;
        s.assign(b + p, q - p);
        if (q < n) p = q + 1;
        else { p = q; eofb = true; if (s.empty()) failb = true; }
    }
};

// mothur gobble(): while(isspace(d=f.get())); if(!f.eof()) f.putback(d);
void gobble(Stream& f) {
    int d;
    for (;;) {
        d = f.get();
        char dc = (char)d;
        if (isspace((int)(unsigned char)dc) && d != EOF) continue;
        break;
    }
    if (!f.eofb) f.putback();
}

std::string getCommentString(Stream& f) {
    std::string temp;
    while (f.ok()) {
        char letter = (char)f.get();
        if (letter == '\r' || letter == '\n' || letter == (char)-1) { gobble(f); break; }
        temp += letter;
    }
    return temp;
}


// byte -> mothur sequence character: toupper, U->T, anything not in .-ACGTN -> N
static unsigned char g_map[256];
static void initMap() {
    for (int c = 0; c < 256; c++) {
        int u = (c < 128) ? toupper(c) : c;
        if (u == 'U') u = 'T';
        if (u != '.' && u != '-' && u != 'A' && u != 'T' && u != 'G' && u != 'C' && u != 'N') u = 'N';
        g_map[c] = (unsigned char)u;
    }
}

// getSequenceString: while (peek != '>' && peek != EOF) { getline; gobble; rtrim; map; append }
// Works on the buffer directly; same stream-state transitions as the istream version.
void getSequenceString(Stream& f, std::string& seq) {
    seq.clear();
    for (;;) {
        int pk = f.peek();
        if (pk == '>' || pk == EOF) break;
        const char* base = f.b; size_t p = f.p, n = f.n;
        const char* nl = (const char*)memchr(base + p, '\n', n - p);
        size_t q = nl ? (size_t)(nl - base) : n;
        size_t ls = p, le = q;
        if (nl) f.p = q + 1; else { f.p = n; f.eofb = true; if (q == p) f.failb = true; }
        gobble(f);
        while (le > ls && isspace((unsigned char)base[le - 1])) --le;
        size_t old = seq.size(); seq.resize(old + (le - ls));
        char* d = &seq[old];
        for (size_t k = ls; k < le; k++) *d++ = (char)g_map[(unsigned char)base[k]];
    }
}

struct Rec { std::string name, aligned; };

// returns false on mothur's "Blank name" error
bool readSequence(Stream& f, Rec& r) {
    r.name.clear(); r.aligned.clear();
    std::string name;
    f.token(name);
    if (!name.empty()) {
        name = name.substr(1);
        for (auto& c : name) if (c == ':') c = '_';
    } else {
        if (!f.eofb) return false;
    }
    // comment records
    while ((name.empty() ? '\0' : name[0]) == '#' && f.ok()) {
        while (!f.eofb) { char c = (char)f.get(); if (c == 10 || c == 13) break; }
        getCommentString(f);
        if (f.ok()) { f.token(name); name = name.size() >= 1 ? name.substr(1) : std::string(); }
        else { name.clear(); break; }
    }
    getCommentString(f);
    std::string& al = r.aligned;
    getSequenceString(f, al);
    // setAligned (characters are already upper case)
    if (!al.empty() && al[0] == '-') {
        for (size_t i = 0; i < al.size(); i++) { if (al[i] == '-') al[i] = '.'; else break; }
        for (size_t i = al.size(); i-- > 0;) { if (al[i] == '-') al[i] = '.'; else break; }
    }
    r.name = std::move(name);
    return true;
}

// ------------------------------------------------------------------ exact scalar calculators
// Direct port of mothur 1.48.6 eachGapIgnoreTermGapDist (same control flow, same arithmetic); used for the
// rare countends=F eachgap pairs whose maxMinLength reaches zero before the last processed column.
inline int setStartIT(const char* A, const char* B, int L, bool& ov) {
    for (int i = 0; i < L; i++) if (A[i] != '.' && B[i] != '.' && A[i] != '-' && B[i] != '-') { ov = true; return i; }
    return 0;
}
inline int setEndIT(const char* A, const char* B, int L, bool& ov) {
    for (int i = L - 1; i >= 0; i--) if (A[i] != '.' && B[i] != '.' && A[i] != '-' && B[i] != '-') { ov = true; return i; }
    return 0;
}




double eachGapIT(const char* A, const char* B, int L, double cutoff) {
    bool overlap = false;
    int start = setStartIT(A, B, L, overlap), end = setEndIT(A, B, L, overlap);
    if (!overlap) return 1.0;
    int mml = end - start + 1, diff = 0; double dist;
    for (int i = start; i < L; i++) {
        char a = A[i], b = B[i];
        if (a == '.' || b == '.') break;
        else if ((a == '-' && b == '-') || (a == '-' && b == '.') || (a == '.' && b == '-')) { mml--; }
        else { if (a != b) diff++; }
        dist = (double)diff / mml;
        if (dist > cutoff) return 1.0;
    }
    if (mml == 0) dist = 1.0; else dist = ((double)diff / (double)mml);
    return dist;
}



// ------------------------------------------------------------------ output
struct Out {
    FILE* f; char* b; size_t len = 0, cap = 1 << 22;
    explicit Out(FILE* fp) : f(fp) { b = (char*)malloc(cap); }
    ~Out() { free(b); }
    Out(const Out&) = delete; Out& operator=(const Out&) = delete;
    inline void reserve(size_t n) {
        if (len + n <= cap) return;
        if (f && len) { flush(); if (n <= cap) return; }
        while (len + n > cap) cap *= 2;
        b = (char*)realloc(b, cap);
    }
    inline void put(const char* s, size_t n) { reserve(n); memcpy(b + len, s, n); len += n; }
    inline void put(const std::string& s) { put(s.data(), s.size()); }
    inline void putc_(char c) { reserve(1); b[len++] = c; }
    // a cached number: slot[0] = length (<= 15), slot[1..15] = characters; copies 16 bytes, advances by length
    inline void putSlot(const char* slot) { reserve(16); memcpy(b + len, slot + 1, 15); len += (unsigned char)slot[0]; }
    void flush() { if (len) { fwrite(b, 1, len, f); len = 0; } }
};

inline int fmtg(char* dst, double d, int prec) { return snprintf(dst, 64, "%.*g", prec, d); }

// ------------------------------------------------------------------ bit-parallel engine
// Every calculator reduces to two integers per pair: diff and maxMinLength (mml) over a
// column range; dist = diff/mml. mothur's running early exit ("return 1.0 once diff/mml
// exceeds cutoff") is equivalent to testing the final ratio, because diff never decreases
// and mml never increases, as long as the final mml is positive. The one case where mml can
// reach zero or below mid-loop with other columns still to come (eachgap, countends=F, with
// '-' columns after the last shared base) is sent to the exact scalar port.
enum Calc { ONEGAP_T = 0, EACHGAP_T = 1, NOGAPS = 2, EACHGAP_F = 3, ONEGAP_F = 4 };

struct Seqs {
    long long n = 0; int L = 0, W = 0;
    std::vector<uint64_t> pl;               // per sequence: P0,P1,P2 (symbol code bits), D ('.'), H ('-'), G (D|H), each W words
    std::vector<int> fnd, lnd, fb, lb;      // first/last non-'.', first/last base (not '.' or '-'); -1 if none
    std::vector<uint8_t> idot;              // has a '.' strictly between fnd and lnd
    std::vector<std::string>* al = nullptr; // aligned strings (for the scalar fallback)
    std::vector<uint64_t> tl;               // lane-transposed planes: [block of 8 sequences][plane][word][lane]
    const uint64_t* tb(long long blk) const { return tl.data() + (size_t)blk * 6 * W * 8; }
    const uint64_t* p(long long i) const { return pl.data() + (size_t)i * 6 * W; }
};

static inline uint64_t loMask0(int w, int lo) { int b = lo - (w << 6); return b <= 0 ? ~0ULL : (b >= 64 ? 0 : (~0ULL << b)); }
static inline uint64_t hiMask0(int w, int hi) { int b = hi - (w << 6); return b >= 63 ? ~0ULL : (b < 0 ? 0 : ((2ULL << b) - 1)); }

// Bit planes of one aligned sequence (codes: A1 C2 G3 T4 N5 -6 .7; planes P0,P1,P2,D,H,G).
static void planesScalar(const char* s, int L, uint64_t* q, int W) {
    for (int c = 0; c < L; c++) {
        unsigned code;
        switch (s[c]) { case 'A': code = 1; break; case 'C': code = 2; break; case 'G': code = 3; break; case 'T': code = 4; break;
                        case 'N': code = 5; break; case '-': code = 6; break; default: code = 7; break; }
        uint64_t bit = 1ULL << (c & 63); int w = c >> 6;
        if (code & 1) q[w] |= bit;
        if (code & 2) q[W + w] |= bit;
        if (code & 4) q[2 * W + w] |= bit;
        if (code == 7) q[3 * W + w] |= bit;
        if (code == 6) q[4 * W + w] |= bit;
        if (code >= 6) q[5 * W + w] |= bit;
    }
}
#pragma GCC push_options
#pragma GCC target("avx2")
static void planesAVX2(const char* s, int L, uint64_t* q, int W) {
    alignas(32) char buf[64];
    for (int w = 0; w < W; w++) {
        int len = L - 64 * w; const char* src = s + 64 * w;
        if (len < 64) { memset(buf, 0, 64); memcpy(buf, src, len); src = buf; }
        __m256i lo = _mm256_loadu_si256((const __m256i*)src), hi = _mm256_loadu_si256((const __m256i*)(src + 32));
        auto mm = [&](char ch) -> uint64_t {
            __m256i v = _mm256_set1_epi8(ch);
            return (uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(lo, v)) | ((uint64_t)(uint32_t)_mm256_movemask_epi8(_mm256_cmpeq_epi8(hi, v)) << 32);
        };
        uint64_t A = mm('A'), C = mm('C'), G = mm('G'), T = mm('T'), N = mm('N'), H = mm('-'), D = mm('.');
        q[w] = A | G | N | D; q[W + w] = C | G | H | D; q[2 * W + w] = T | N | H | D;
        q[3 * W + w] = D; q[4 * W + w] = H; q[5 * W + w] = D | H;
    }
}
#pragma GCC pop_options

void buildSeqs(Seqs& S, std::vector<std::string>& al, int L) {
    S.n = (long long)al.size(); S.L = L; S.W = (L + 63) / 64; S.al = &al;
    int W = S.W;
    S.pl.assign((size_t)S.n * 6 * W, 0);
    S.fnd.assign(S.n, -1); S.lnd.assign(S.n, -1); S.fb.assign(S.n, -1); S.lb.assign(S.n, -1); S.idot.assign(S.n, 0);
    __builtin_cpu_init();
    const char* f = getenv("SWIFTDISTSEQS_ISA");
    bool avx2 = __builtin_cpu_supports("avx2") && !(f && (!strcmp(f, "base") || !strcmp(f, "popcnt")));
    for (long long i = 0; i < S.n; i++) {
        uint64_t* q = S.pl.data() + (size_t)i * 6 * W;
        if (avx2) planesAVX2(al[i].data(), L, q, W); else planesScalar(al[i].data(), L, q, W);
        int fnd = -1, lnd = -1, fb = -1, lb = -1;
        for (int w = 0; w < W; w++) {
            uint64_t valid = hiMask0(w, L - 1);
            uint64_t nd = ~q[3 * W + w] & valid, bs = ~q[5 * W + w] & valid;
            if (nd) { if (fnd < 0) fnd = (w << 6) + __builtin_ctzll(nd); lnd = (w << 6) + 63 - __builtin_clzll(nd); }
            if (bs) { if (fb < 0) fb = (w << 6) + __builtin_ctzll(bs); lb = (w << 6) + 63 - __builtin_clzll(bs); }
        }
        S.fnd[i] = fnd; S.lnd[i] = lnd; S.fb[i] = fb; S.lb[i] = lb;
        if (fnd >= 0) {
            for (int w = (fnd + 1) >> 6; w <= ((lnd - 1) >> 6) && lnd - 1 >= fnd + 1; w++) {
                uint64_t x = q[3 * W + w] & loMask0(w, fnd + 1) & hiMask0(w, lnd - 1);
                if (x) { S.idot[i] = 1; break; }
            }
        }
    }
    // padding so that 8-lane loads past the last sequence stay inside the arrays
    S.fnd.resize(S.n + 8, -1); S.lnd.resize(S.n + 8, -1); S.fb.resize(S.n + 8, -1); S.lb.resize(S.n + 8, -1); S.idot.resize(S.n + 8, 1);
    long long nb = (S.n + 7) / 8;
    S.tl.assign((size_t)nb * 6 * W * 8, 0);
    for (long long j = 0; j < S.n; j++) {
        const uint64_t* q = S.pl.data() + (size_t)j * 6 * W;
        uint64_t* t = S.tl.data() + (size_t)(j / 8) * 6 * W * 8 + (j % 8);
        for (int k = 0; k < 6; k++) for (int w = 0; w < W; w++) t[((size_t)k * W + w) * 8] = q[(size_t)k * W + w];
    }
}

static inline uint64_t loMask(int w, int lo) { int b = lo - (w << 6); return b <= 0 ? ~0ULL : (b >= 64 ? 0 : (~0ULL << b)); }
static inline uint64_t hiMask(int w, int hi) { int b = hi - (w << 6); return b >= 63 ? ~0ULL : (b < 0 ? 0 : ((2ULL << b) - 1)); }

struct KR { int diff, mml; double v; int kind; };  // kind 0: ratio diff/mml (1/1 encodes 1.0); 1: excluded (> cutoff); 2: value v

// formatted-number cache keyed by (diff, mml)
struct FmtCache {
    int prec; std::vector<std::vector<std::array<char, 16>>> rows;
    explicit FmtCache(int p, int L) : prec(p), rows((size_t)L + 2) {}
    // slot for diff/mml (formatted on first use), or nullptr if not cacheable
    inline const char* slot(int diff, int mml, double d) {
        if (mml >= 1 && (size_t)mml < rows.size() && diff >= 0 && diff <= mml) {
            auto& r = rows[mml];
            if (r.empty()) { r.resize((size_t)mml + 1); for (auto& x : r) x[0] = 0; }
            char* e = r[diff].data();
            if (e[0] == 0) {
                char tmp[64]; int k = fmtg(tmp, d, prec);
                if (k > 15) return nullptr;
                memset(e, 0, 16); e[0] = (char)k; memcpy(e + 1, tmp, k);
            }
            return e;
        }
        return nullptr;
    }
};

std::string rootName(const std::string& path) {
    size_t sl = path.find_last_of('/');
    std::string simple = (sl == std::string::npos) ? path : path.substr(sl + 1);
    size_t dot = simple.find_last_of('.');
    if (dot == std::string::npos) return simple + ".";
    return simple.substr(0, dot + 1);
}

static double nowSec() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static bool g_timing = false;
static void tmark(const char* what) { static double t0 = nowSec(), last = t0; if (!g_timing) return; double t = nowSec(); fprintf(stderr, "[time] %-8s %.4f s (total %.4f)\n", what, t - last, t - t0); last = t; }

int fail(const char* msg) { fprintf(stderr, "[ERROR]: %s\n", msg); return 1; }

// ------------------------------------------------------------------ drivers
struct Job {
    const Seqs* S; const std::vector<Rec>* recs; double cutoff; int output;  // 0 column, 1 lt, 2 square
    Out* out;
};

typedef void (*RowsFn)(const Job&, long long, long long);

namespace eng_base {
#define ENG_ISA 0
#include "engine.inc"
#undef ENG_ISA
}
#pragma GCC push_options
#pragma GCC target("popcnt")
namespace eng_popcnt {
#define ENG_ISA 1
#include "engine.inc"
#undef ENG_ISA
}
#pragma GCC pop_options
#pragma GCC push_options
#pragma GCC target("popcnt,avx,avx2,bmi,bmi2")
namespace eng_avx2 {
#define ENG_ISA 2
#include "engine.inc"
#undef ENG_ISA
}
#pragma GCC pop_options
#pragma GCC push_options
#pragma GCC target("popcnt,avx,avx2,bmi,bmi2,avx512f,avx512bw,avx512vl,avx512dq,avx512vpopcntdq")
namespace eng_avx512 {
#define ENG_ISA 3
#include "engine.inc"
#undef ENG_ISA
}
#pragma GCC pop_options

// Run-time choice of the widest supported variant; SWIFTDISTSEQS_ISA=base|popcnt|avx2|avx512 forces one (testing).
const char* g_isa = "base";
RowsFn pickRows(int calc, bool early) {
    __builtin_cpu_init();
    int best = 0;
    if (__builtin_cpu_supports("popcnt")) best = 1;
    if (best >= 1 && __builtin_cpu_supports("avx2") && __builtin_cpu_supports("bmi2")) best = 2;
    if (best >= 2 && __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512dq") && __builtin_cpu_supports("avx512vpopcntdq")) best = 3;
    const char* f = getenv("SWIFTDISTSEQS_ISA");
    if (f) {
        int want = !strcmp(f, "base") ? 0 : !strcmp(f, "popcnt") ? 1 : !strcmp(f, "avx2") ? 2 : !strcmp(f, "avx512") ? 3 : best;
        if (want < best) best = want;
    }
    static const char* names[] = {"base", "popcnt", "avx2", "avx512"};
    g_isa = names[best];
    switch (best) {
        case 3: return eng_avx512::pickRows(calc, early);
        case 2: return eng_avx2::pickRows(calc, early);
        case 1: return eng_popcnt::pickRows(calc, early);
        default: return eng_base::pickRows(calc, early);
    }
}

}  // namespace

int main(int argc, char** argv) {
    Opts o;
    if (argc == 1 || (argc == 2 && (!strcmp(argv[1], "--help") || !strcmp(argv[1], "-h")))) {
        fprintf(stderr,
            "SwiftDistSeqs: byte-identical, fast mothur 1.48.6 dist.seqs for aligned DNA.\n"
            "  distseqs --fasta IN --out OUT [--calc onegap|eachgap|nogaps] [--countends T|F]\n"
            "           [--cutoff X] [--output column|lt|square] [--processors P]\n"
            "  distseqs fasta=IN [calc=..] [countends=..] [cutoff=..] [output=..] [processors=..] [outputdir=..]\n"
            "           (mothur-style; writes IN-root.dist / .phylip.dist / .square.dist like mothur)\n");
        return argc == 1 ? 64 : 0;
    }
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        auto next = [&](std::string& dst) { if (i + 1 < argc) dst = argv[++i]; };
        if (a == "--fasta") next(o.fasta);
        else if (a == "--out") next(o.out);
        else if (a == "--calc") next(o.calc);
        else if (a == "--countends") next(o.countends);
        else if (a == "--cutoff") next(o.cutoff);
        else if (a == "--output") next(o.output);
        else if (a == "--processors") { std::string p; next(p); o.procs = atoi(p.c_str()); }
        else if (a.find('=') != std::string::npos && a.compare(0, 2, "--") != 0) {
            o.mothurStyle = true;
            std::string k = a.substr(0, a.find('=')), v = a.substr(a.find('=') + 1);
            if (k == "fasta") o.fasta = v; else if (k == "calc") o.calc = v; else if (k == "countends") o.countends = v;
            else if (k == "cutoff") o.cutoff = v; else if (k == "output") o.output = v; else if (k == "processors") o.procs = atoi(v.c_str());
            else if (k == "outputdir") o.outputdir = v;
            else { fprintf(stderr, "unknown parameter %s\n", k.c_str()); return 64; }
        } else { fprintf(stderr, "unknown option %s\n", a.c_str()); return 64; }
    }
    if (o.fasta.empty()) return fail("fasta is required");
    if (o.procs < 1) o.procs = 1;
    // option semantics as in DistanceCommand
    std::string calc = o.calc;
    if (calc == "default") calc = "onegap";
    if (calc != "onegap" && calc != "eachgap" && calc != "nogaps") calc = "onegap";
    std::string cu = upper(o.countends);
    bool countends = (cu == "TRUE" || cu == "T");
    // mothur stores cutoff as float: istringstream >> float for strings of [0-9.-], else (float)atof
    float cutoff_f;
    if (!o.cutoff.empty() && o.cutoff.find_first_not_of("0123456789.-") == std::string::npos) cutoff_f = strtof(o.cutoff.c_str(), nullptr);
    else cutoff_f = (float)atof(o.cutoff.c_str());
    double cutoff = (double)cutoff_f;
    std::string output = o.output;
    if (output == "phylip") output = "lt";
    if (output != "column" && output != "lt" && output != "square") output = "column";
    if (o.out.empty()) {
        std::string dir = o.outputdir;
        if (dir.empty()) { size_t sl = o.fasta.find_last_of('/'); if (sl != std::string::npos) dir = o.fasta.substr(0, sl + 1); }
        else if (dir.back() != '/') dir += '/';
        std::string root = dir + rootName(o.fasta);
        o.out = root + (output == "lt" ? "phylip.dist" : output == "square" ? "square.dist" : "dist");
    }

    g_timing = getenv("SWIFTDISTSEQS_TIMING") != nullptr; tmark("start");
    initMap();
    // read input
    FILE* in = fopen(o.fasta.c_str(), "rb");
    if (!in) return fail("cannot open fasta file");
    std::vector<char> data;
    {
        long sz = -1;
        if (fseek(in, 0, SEEK_END) == 0) { sz = ftell(in); if (fseek(in, 0, SEEK_SET) != 0) sz = -1; }
        if (sz >= 0) { data.resize((size_t)sz); size_t got = sz ? fread(data.data(), 1, (size_t)sz, in) : 0; data.resize(got); }
        char tmp[1 << 16]; size_t r;  // anything beyond the measured size (or a non-seekable input)
        while ((r = fread(tmp, 1, sizeof tmp, in)) > 0) data.insert(data.end(), tmp, tmp + r);
    }
    fclose(in);
    tmark("read");
    Stream st(data.data(), data.size());
    std::vector<Rec> recs;
    size_t length = 0; bool same = true;
    while (!st.eofb) {
        Rec r;
        if (!readSequence(st, r)) return fail("Blank name in fasta file");
        if (!r.name.empty()) {
            if (length == 0) length = r.aligned.size();
            if (length != r.aligned.size()) same = false;
            recs.push_back(std::move(r));
        }
        gobble(st);
    }
    if (!same) return fail("your sequences are not the same length, aborting.");
    long long n = (long long)recs.size();
    if (n < 2) return fail("you must have at least 2 sequences to calculate the distances, aborting.");
    int L = (int)length;

    int ci;
    if (calc == "nogaps") ci = NOGAPS;
    else if (calc == "eachgap") ci = countends ? EACHGAP_T : EACHGAP_F;
    else ci = countends ? ONEGAP_T : ONEGAP_F;

    std::vector<std::string> al(n);
    for (long long i = 0; i < n; i++) al[i] = std::move(recs[i].aligned);
    tmark("parse");
    Seqs S;
    buildSeqs(S, al, L);
    tmark("planes");

    FILE* fo = fopen(o.out.c_str(), "wb");
    if (!fo) return fail("cannot open output file");
    Out out(fo);
    int oi = output == "column" ? 0 : (output == "lt" ? 1 : 2);
    Job J{&S, &recs, cutoff, oi, &out};
    if (oi != 0) { char num[64]; int k = snprintf(num, sizeof num, "%lld\n", n); out.put(num, k); }
    RowsFn rows = pickRows(ci, oi == 0);
    const long long r0 = (oi == 0) ? 1 : 0;
    int T = o.procs;
    if (T > 1 && n - r0 >= 2) {
        // Rows are cut into blocks of similar work; workers compute blocks into private buffers
        // and the main thread writes them strictly in row order (same bytes as one thread).
        out.flush();
        std::vector<std::pair<long long, long long>> blocks;
        double total = (oi == 2) ? (double)n * (double)n : (double)n * (double)(n - 1) / 2.0;
        double target = std::max(total / (T * 32.0), 20000.0);
        long long s0 = r0; double acc = 0;
        for (long long i = r0; i < n; i++) {
            acc += (oi == 2) ? (double)n : (double)i;
            if (acc >= target || i == n - 1) { blocks.push_back({s0, i + 1}); s0 = i + 1; acc = 0; }
        }
        const long long nb = (long long)blocks.size();
        std::vector<std::unique_ptr<Out>> bufs(nb);
        std::vector<char> ready(nb, 0);
        std::mutex mu; std::condition_variable cv;
        long long written = 0; std::atomic<long long> next{0};
        const long long WINDOW = 4LL * T + 4;
        auto worker = [&]() {
            for (;;) {
                long long b = next++;
                if (b >= nb) return;
                { std::unique_lock<std::mutex> lk(mu); cv.wait(lk, [&] { return b < written + WINDOW; }); }
                std::unique_ptr<Out> ob(new Out(nullptr));
                Job Jb = J; Jb.out = ob.get();
                rows(Jb, blocks[b].first, blocks[b].second);
                { std::lock_guard<std::mutex> lk(mu); bufs[b] = std::move(ob); ready[b] = 1; }
                cv.notify_all();
            }
        };
        std::vector<std::thread> th;
        for (int t = 0; t < T; t++) th.emplace_back(worker);
        for (long long b = 0; b < nb; b++) {
            std::unique_ptr<Out> ob;
            { std::unique_lock<std::mutex> lk(mu); cv.wait(lk, [&] { return ready[b] != 0; }); ob = std::move(bufs[b]); }
            if (ob->len) fwrite(ob->b, 1, ob->len, fo);
            { std::lock_guard<std::mutex> lk(mu); written = b + 1; }
            cv.notify_all();
        }
        for (auto& t : th) t.join();
    } else {
        rows(J, r0, n);
    }
    out.flush();
    if (fclose(fo) != 0) return fail("write error");
    tmark("compute+write"); if (g_timing) fprintf(stderr, "[isa] %s\n", g_isa);
    return 0;
}
