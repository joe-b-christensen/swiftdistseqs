// Engine body, included once per instruction-set variant inside a "#pragma GCC target" region
// and its own namespace (see swiftdistseqs.cpp). ENG_ISA: 0 baseline x86-64, 1 +popcnt, 2 +AVX2, 3 +AVX-512 VPOPCNTDQ.

// Lower bound on the number of differences: columns where both sequences are non-gap and differ.
// Returns true once the bound reaches thr. Columns outside [lo,hi] never count (callers guarantee
// that no column outside can be both non-gap), so whole words are used.
static inline bool lbReaches(const uint64_t* a, const uint64_t* b, int W, int w0, int w1, int thr) {
#if ENG_ISA == 3
    __m512i acc = _mm512_setzero_si512();
    for (int w = w0; w <= w1; w += 8) {
        int nw = w1 - w + 1; __mmask8 k = nw >= 8 ? (__mmask8)0xff : (__mmask8)((1u << nw) - 1);
        __m512i x = _mm512_xor_si512(_mm512_maskz_loadu_epi64(k, a + w), _mm512_maskz_loadu_epi64(k, b + w));
        x = _mm512_ternarylogic_epi64(x, _mm512_maskz_loadu_epi64(k, a + W + w), _mm512_maskz_loadu_epi64(k, b + W + w), 0xF6);      // x | (y ^ z)
        x = _mm512_ternarylogic_epi64(x, _mm512_maskz_loadu_epi64(k, a + 2 * W + w), _mm512_maskz_loadu_epi64(k, b + 2 * W + w), 0xF6);
        x = _mm512_ternarylogic_epi64(x, _mm512_maskz_loadu_epi64(k, a + 5 * W + w), _mm512_maskz_loadu_epi64(k, b + 5 * W + w), 0x10);  // x & ~(y | z)
        acc = _mm512_add_epi64(acc, _mm512_popcnt_epi64(x));
        if (_mm512_reduce_add_epi64(acc) >= thr) return true;
    }
    return false;
#elif ENG_ISA == 2
    const __m256i lut = _mm256_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4, 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
    const __m256i low = _mm256_set1_epi8(0x0f);
    int lb = 0, w = w0;
    for (; w + 3 <= w1; w += 4) {
        __m256i x = _mm256_or_si256(_mm256_or_si256(
            _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(a + w)), _mm256_loadu_si256((const __m256i*)(b + w))),
            _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(a + W + w)), _mm256_loadu_si256((const __m256i*)(b + W + w)))),
            _mm256_xor_si256(_mm256_loadu_si256((const __m256i*)(a + 2 * W + w)), _mm256_loadu_si256((const __m256i*)(b + 2 * W + w))));
        x = _mm256_andnot_si256(_mm256_or_si256(_mm256_loadu_si256((const __m256i*)(a + 5 * W + w)), _mm256_loadu_si256((const __m256i*)(b + 5 * W + w))), x);
        __m256i c = _mm256_add_epi8(_mm256_shuffle_epi8(lut, _mm256_and_si256(x, low)), _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi16(x, 4), low)));
        __m256i s = _mm256_sad_epu8(c, _mm256_setzero_si256());
        lb += (int)(_mm256_extract_epi64(s, 0) + _mm256_extract_epi64(s, 1) + _mm256_extract_epi64(s, 2) + _mm256_extract_epi64(s, 3));
        if (lb >= thr) return true;
    }
    for (; w <= w1; w++) {
        lb += __builtin_popcountll(~(a[5 * W + w] | b[5 * W + w]) & ((a[w] ^ b[w]) | (a[W + w] ^ b[W + w]) | (a[2 * W + w] ^ b[2 * W + w])));
        if (lb >= thr) return true;
    }
    return false;
#else
    int lb = 0;
    for (int w = w0; w <= w1; w++) {
        lb += __builtin_popcountll(~(a[5 * W + w] | b[5 * W + w]) & ((a[w] ^ b[w]) | (a[W + w] ^ b[W + w]) | (a[2 * W + w] ^ b[2 * W + w])));
        if (lb >= thr) return true;
    }
    return false;
#endif
}

// first column >= from (and <= to) where both sequences have a base; -1 if none
static inline int firstBoth(const Seqs& S, const uint64_t* a, const uint64_t* b, int from, int to) {
    int W = S.W;
    for (int w = from >> 6; w <= (to >> 6); w++) {
        uint64_t x = ~(a[3 * W + w] | a[4 * W + w] | b[3 * W + w] | b[4 * W + w]) & loMask(w, from) & hiMask(w, to);
        if (x) return (w << 6) + __builtin_ctzll(x);
    }
    return -1;
}
static inline int lastBoth(const Seqs& S, const uint64_t* a, const uint64_t* b, int from, int to) {
    int W = S.W;
    for (int w = to >> 6; w >= (from >> 6); w--) {
        uint64_t x = ~(a[3 * W + w] | a[4 * W + w] | b[3 * W + w] | b[4 * W + w]) & loMask(w, from) & hiMask(w, to);
        if (x) return (w << 6) + 63 - __builtin_clzll(x);
    }
    return -1;
}
// first column >= from where either (any=true) or both (any=false) sequences have '.'; L if none
static inline int firstDot(const Seqs& S, const uint64_t* a, const uint64_t* b, int from, bool any) {
    int W = S.W;
    for (int w = from >> 6; w < W; w++) {
        uint64_t x = (any ? (a[3 * W + w] | b[3 * W + w]) : (a[3 * W + w] & b[3 * W + w])) & loMask(w, from);
        if (x) return (w << 6) + __builtin_ctzll(x);
    }
    return S.L;
}

template <int C, bool EARLY>
static inline KR kernel(const Seqs& S, long long ia, long long ib, double cutoff) {
    const int W = S.W;
    const uint64_t* a = S.p(ia); const uint64_t* b = S.p(ib);
    int start, end, last;  // process columns [start, last], truncated at a break column; mml base = end-start+1 (or end-start)
    bool truncBoth = false, truncAny = false, early = EARLY;
    if (C == ONEGAP_T || C == EACHGAP_T) {
        int fa = S.fnd[ia], fb2 = S.fnd[ib];
        if (fa < 0 && fb2 < 0) { start = 0; end = 0; }
        else {
            start = fa < 0 ? fb2 : (fb2 < 0 ? fa : (fa < fb2 ? fa : fb2));
            int la = S.lnd[ia], lb2 = S.lnd[ib];
            end = la > lb2 ? la : lb2;
        }
        last = end; truncBoth = true;
        if (S.L == 0) return KR{0, 1, 0.0, 0};
    } else if (C == NOGAPS) {
        if (!S.idot[ia] && !S.idot[ib]) {
            if (S.fnd[ia] < 0 || S.fnd[ib] < 0) return KR{1, 1, 1.0, 0};
            start = S.fnd[ia] > S.fnd[ib] ? S.fnd[ia] : S.fnd[ib];
            end = S.lnd[ia] < S.lnd[ib] ? S.lnd[ia] : S.lnd[ib];
            if (start > end) return KR{1, 1, 1.0, 0};
        } else {
            // generic: first/last column where neither is '.'
            start = -1; end = -1;
            for (int w = 0; w < W && start < 0; w++) { uint64_t x = ~(a[3 * W + w] | b[3 * W + w]) & hiMask(w, S.L - 1); if (x) start = (w << 6) + __builtin_ctzll(x); }
            if (start < 0) return KR{1, 1, 1.0, 0};
            for (int w = W - 1; w >= 0 && end < 0; w--) { uint64_t x = ~(a[3 * W + w] | b[3 * W + w]) & hiMask(w, S.L - 1); if (x) end = (w << 6) + 63 - __builtin_clzll(x); }
        }
        last = end; truncAny = true;
    } else {  // EACHGAP_F, ONEGAP_F: first/last column where both are bases
        int lo = S.fb[ia] > S.fb[ib] ? S.fb[ia] : S.fb[ib];
        int hi = S.lb[ia] < S.lb[ib] ? S.lb[ia] : S.lb[ib];
        if (S.fb[ia] < 0 || S.fb[ib] < 0 || lo > hi) return KR{1, 1, 1.0, 0};
        start = firstBoth(S, a, b, lo, hi);
        if (start < 0) return KR{1, 1, 1.0, 0};
        end = lastBoth(S, a, b, start, hi);
        if (C == ONEGAP_F) { last = end; }
        else {
            int brk = (!S.idot[ia] && !S.idot[ib]) ? ((S.lnd[ia] < S.lnd[ib] ? S.lnd[ia] : S.lnd[ib]) + 1) : firstDot(S, a, b, start, true);
            if (brk > S.L) brk = S.L;
            last = brk - 1;
            if (brk > end + 1) early = false;  // mml may reach zero: no lower-bound exit
        }
    }
    const int range = (C == ONEGAP_F) ? (end - start) : (end - start + 1);
    int thr = 0x7fffffff;
    if (early && range > 0) {
        double t = cutoff * (double)range;
        if (t < 2147483000.0) thr = (int)__builtin_floor(t) + 2;
    }
    // Cheap exact lower bound: columns where both are non-gap and differ always count as a
    // difference inside the processed range when neither sequence has an internal '.'.
    if (EARLY && early && thr != 0x7fffffff && !S.idot[ia] && !S.idot[ib]) {
        int lo = (C == ONEGAP_F || C == EACHGAP_F) ? start : (S.fnd[ia] > S.fnd[ib] ? S.fnd[ia] : S.fnd[ib]);
        int hi = (C == ONEGAP_F || C == EACHGAP_F) ? end : (S.lnd[ia] < S.lnd[ib] ? S.lnd[ia] : S.lnd[ib]);
        if (lo >= 0 && lo <= hi) {
            if (lbReaches(a, b, W, lo >> 6, hi >> 6, thr)) return KR{0, 0, 0.0, 1};
        }
    }
    int diff = 0, dec = 0;
    uint64_t pX = 0, pY = 0, cX = 0, cY = 0;
    const int ws = start >> 6, we = last >> 6;
    for (int w = ws; w <= we; w++) {
        uint64_t m = loMask(w, start) & hiMask(w, last);
        const uint64_t dA = a[3 * W + w], dB = b[3 * W + w];
        bool stop = false;
        if (truncBoth) { uint64_t bd = dA & dB & m; if (bd) { m &= (bd & (0 - bd)) - 1; stop = true; } }
        if (truncAny) { uint64_t bd = (dA | dB) & m; if (bd) { m &= (bd & (0 - bd)) - 1; stop = true; } }
        const uint64_t neq = (a[w] ^ b[w]) | (a[W + w] ^ b[W + w]) | (a[2 * W + w] ^ b[2 * W + w]);
        if (C == EACHGAP_T) {
            uint64_t gg = (dA | a[4 * W + w]) & (dB | b[4 * W + w]) & m;
            diff += __builtin_popcountll(neq & m & ~gg); dec += __builtin_popcountll(gg);
        } else if (C == NOGAPS) {
            uint64_t dash = (a[4 * W + w] | b[4 * W + w]) & m;
            diff += __builtin_popcountll(neq & m & ~dash); dec += __builtin_popcountll(dash);
        } else if (C == EACHGAP_F) {
            uint64_t gg = a[4 * W + w] & b[4 * W + w] & m;
            diff += __builtin_popcountll(neq & m & ~gg); dec += __builtin_popcountll(gg);
        } else {  // ONEGAP_T, ONEGAP_F
            uint64_t gA, gB;
            if (C == ONEGAP_T) { gA = dA | a[4 * W + w]; gB = dB | b[4 * W + w]; }
            else { gA = a[4 * W + w]; gB = b[4 * W + w]; }
            const uint64_t Wm = gA & gB & m, X = gA & ~gB & m, Y = gB & ~gA & m, Z = ~(gA | gB) & m;
            const uint64_t SX = (X << 1) | pX, SY = (Y << 1) | pY;
            pX = X >> 63; pY = Y >> 63;
            unsigned __int128 sx = (unsigned __int128)Wm + (SX & Wm) + cX;
            unsigned __int128 sy = (unsigned __int128)Wm + (SY & Wm) + cY;
            cX = (uint64_t)(sx >> 64); cY = (uint64_t)(sy >> 64);
            const uint64_t RX = SX | ((uint64_t)sx ^ Wm), RY = SY | ((uint64_t)sy ^ Wm);
            diff += __builtin_popcountll(X & ~RX) + __builtin_popcountll(Y & ~RY) + __builtin_popcountll(Z & neq);
            dec += __builtin_popcountll(Wm) + __builtin_popcountll(X & RX) + __builtin_popcountll(Y & RY);
        }
        if (EARLY && diff >= thr) return KR{0, 0, 0.0, 1};
        if (stop) break;
    }
    int mml = range - dec;
    if (mml > 0) {
        double d = (double)diff / mml;
        if (d > cutoff) return KR{1, 1, 1.0, 0};
        return KR{diff, mml, d, 0};
    }
    if (C != EACHGAP_F) return KR{1, 1, 1.0, 0};
    const std::string& sa = (*S.al)[ia]; const std::string& sb = (*S.al)[ib];
    return KR{0, 0, eachGapIT(sa.data(), sb.data(), S.L, cutoff), 2};
}


#if ENG_ISA == 3
// Eight pairs (row i against j0..j0+7) per vector, for onegap/eachgap with countends=T and nogaps
// when neither sequence has an internal '.'. Lanes it cannot take are left to the scalar kernel.
template <int C, bool EARLY>
static inline void laneBlock(const Seqs& S, long long i, long long j0, int nl, double cutoff, KR* res, unsigned& doneMask, unsigned& exclMask) {
    const int W = S.W;
    __mmask8 act = (__mmask8)(nl >= 8 ? 0xff : ((1u << nl) - 1));
    unsigned idB = 0;
    for (int l = 0; l < 8; l++) idB |= (unsigned)(S.idot[j0 + l] != 0) << l;
    __mmask8 ok = S.idot[i] ? (__mmask8)0 : (__mmask8)(act & ~idB);
    doneMask = ok; exclMask = 0;
    if (!ok) return;
    const __m512i zero = _mm512_setzero_si512(), one = _mm512_set1_epi64(1), ones = _mm512_set1_epi64(-1);
    const __m512i fB = _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i*)(S.fnd.data() + j0)));
    const __m512i lB = _mm512_cvtepi32_epi64(_mm256_loadu_si256((const __m256i*)(S.lnd.data() + j0)));
    const long long fa = S.fnd[i], la = S.lnd[i];
    const __mmask8 eB = _mm512_cmplt_epi64_mask(fB, zero);
    __m512i lo, hi, range; __mmask8 noov = 0;
    if (C == NOGAPS) {
        if (fa < 0) { noov = 0xff; lo = hi = range = zero; }
        else {
            const __m512i fA = _mm512_set1_epi64(fa), lA = _mm512_set1_epi64(la);
            lo = _mm512_max_epi64(fA, fB); hi = _mm512_min_epi64(lA, lB);
            noov = eB | _mm512_cmpgt_epi64_mask(lo, hi);
            range = _mm512_add_epi64(_mm512_sub_epi64(hi, lo), one);
        }
    } else {
        __m512i start, end, brk;
        if (fa < 0) {
            start = _mm512_mask_mov_epi64(fB, eB, zero);
            end = _mm512_mask_mov_epi64(lB, eB, zero);
            brk = _mm512_mask_mov_epi64(_mm512_add_epi64(lB, one), eB, zero);
        } else {
            const __m512i fA = _mm512_set1_epi64(fa), lA = _mm512_set1_epi64(la);
            start = _mm512_mask_mov_epi64(_mm512_min_epi64(fA, fB), eB, fA);
            end = _mm512_mask_mov_epi64(_mm512_max_epi64(lA, lB), eB, lA);
            const __m512i mnl = _mm512_min_epi64(lA, lB), mxl = _mm512_max_epi64(lA, lB);
            __mmask8 adj = _mm512_cmple_epi64_mask(_mm512_max_epi64(fA, fB), _mm512_add_epi64(mnl, one));
            __m512i b2 = _mm512_mask_mov_epi64(_mm512_add_epi64(mnl, one), adj, _mm512_add_epi64(mxl, one));
            brk = _mm512_mask_mov_epi64(b2, eB, _mm512_set1_epi64(la + 1));
        }
        lo = start; hi = _mm512_sub_epi64(brk, one);
        range = _mm512_add_epi64(_mm512_sub_epi64(end, start), one);
    }
    __mmask8 live = ok & ~noov;
    __mmask8 excl = 0;
    if (EARLY && live && fa >= 0) {
        const __m512d thr = _mm512_add_pd(_mm512_roundscale_pd(_mm512_mul_pd(_mm512_set1_pd(cutoff), _mm512_cvtepi64_pd(range)), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC), _mm512_set1_pd(2.0));
        const uint64_t* a = S.p(i); const uint64_t* tb = S.tb(j0 / 8);
        __m512i lb = zero;
        const int w0 = (int)(fa >> 6), w1 = (int)(la >> 6);
        for (int w = w0; w <= w1; w++) {
            __m512i x = _mm512_xor_si512(_mm512_set1_epi64((long long)a[w]), _mm512_loadu_si512(tb + ((size_t)0 * W + w) * 8));
            x = _mm512_ternarylogic_epi64(x, _mm512_set1_epi64((long long)a[W + w]), _mm512_loadu_si512(tb + ((size_t)1 * W + w) * 8), 0xF6);
            x = _mm512_ternarylogic_epi64(x, _mm512_set1_epi64((long long)a[2 * W + w]), _mm512_loadu_si512(tb + ((size_t)2 * W + w) * 8), 0xF6);
            x = _mm512_ternarylogic_epi64(x, _mm512_set1_epi64((long long)a[5 * W + w]), _mm512_loadu_si512(tb + ((size_t)5 * W + w) * 8), 0x10);
            lb = _mm512_add_epi64(lb, _mm512_popcnt_epi64(x));
            if (((w - w0) & 7) == 7 && (_mm512_cmp_pd_mask(_mm512_cvtepi64_pd(lb), thr, _CMP_GE_OQ) & live) == live) break;
        }
        excl = live & _mm512_cmp_pd_mask(_mm512_cvtepi64_pd(lb), thr, _CMP_GE_OQ);
    }
    exclMask = excl;
    __mmask8 todo = live & ~excl;
    alignas(64) long long dv[8], mv[8];
    if (todo) {
        const uint64_t* a = S.p(i); const uint64_t* tb = S.tb(j0 / 8);
        const int w0 = (int)(_mm512_mask_reduce_min_epi64(todo, lo) >> 6);
        long long hmax = _mm512_mask_reduce_max_epi64(todo, hi);
        const int w1 = hmax < 0 ? -1 : (int)(hmax >> 6);
        __m512i diff = zero, dec = zero, pX = zero, pY = zero, cX = zero, cY = zero;
        for (int w = (w0 < 0 ? 0 : w0); w <= w1; w++) {
            const __m512i base = _mm512_set1_epi64(64LL * w);
            __m512i bl = _mm512_max_epi64(_mm512_sub_epi64(lo, base), zero);
            __m512i bh = _mm512_max_epi64(_mm512_sub_epi64(_mm512_set1_epi64(63), _mm512_sub_epi64(hi, base)), zero);
            __m512i m = _mm512_maskz_and_epi64(todo, _mm512_sllv_epi64(ones, bl), _mm512_srlv_epi64(ones, bh));
            __m512i neq = _mm512_xor_si512(_mm512_set1_epi64((long long)a[w]), _mm512_loadu_si512(tb + ((size_t)0 * W + w) * 8));
            neq = _mm512_ternarylogic_epi64(neq, _mm512_set1_epi64((long long)a[W + w]), _mm512_loadu_si512(tb + ((size_t)1 * W + w) * 8), 0xF6);
            neq = _mm512_ternarylogic_epi64(neq, _mm512_set1_epi64((long long)a[2 * W + w]), _mm512_loadu_si512(tb + ((size_t)2 * W + w) * 8), 0xF6);
            if (C == NOGAPS) {
                __m512i dash = _mm512_and_si512(_mm512_or_si512(_mm512_set1_epi64((long long)a[4 * W + w]), _mm512_loadu_si512(tb + ((size_t)4 * W + w) * 8)), m);
                diff = _mm512_add_epi64(diff, _mm512_popcnt_epi64(_mm512_ternarylogic_epi64(neq, m, dash, 0x40)));   // neq & m & ~dash
                dec = _mm512_add_epi64(dec, _mm512_popcnt_epi64(dash));
            } else {
                const __m512i gA = _mm512_set1_epi64((long long)a[5 * W + w]);
                const __m512i gB = _mm512_loadu_si512(tb + ((size_t)5 * W + w) * 8);
                const __m512i Wm = _mm512_ternarylogic_epi64(gA, gB, m, 0x80);                                          // gA & gB & m
                if (C == EACHGAP_T) {
                    diff = _mm512_add_epi64(diff, _mm512_popcnt_epi64(_mm512_ternarylogic_epi64(neq, m, Wm, 0x40)));   // neq & m & ~Wm
                    dec = _mm512_add_epi64(dec, _mm512_popcnt_epi64(Wm));
                } else {
                    const __m512i X = _mm512_ternarylogic_epi64(gA, gB, m, 0x20);   // gA & ~gB & m
                    const __m512i Y = _mm512_ternarylogic_epi64(gA, gB, m, 0x08);   // ~gA & gB & m
                    const __m512i Z = _mm512_ternarylogic_epi64(gA, gB, m, 0x01 << 1); // ~gA & ~gB & m
                    const __m512i SX = _mm512_or_si512(_mm512_slli_epi64(X, 1), pX), SY = _mm512_or_si512(_mm512_slli_epi64(Y, 1), pY);
                    pX = _mm512_srli_epi64(X, 63); pY = _mm512_srli_epi64(Y, 63);
                    const __m512i tX = _mm512_add_epi64(Wm, _mm512_and_si512(SX, Wm)), tY = _mm512_add_epi64(Wm, _mm512_and_si512(SY, Wm));
                    const __mmask8 c1X = _mm512_cmplt_epu64_mask(tX, Wm), c1Y = _mm512_cmplt_epu64_mask(tY, Wm);
                    const __m512i sX = _mm512_add_epi64(tX, cX), sY = _mm512_add_epi64(tY, cY);
                    const __mmask8 c2X = _mm512_cmplt_epu64_mask(sX, tX), c2Y = _mm512_cmplt_epu64_mask(sY, tY);
                    cX = _mm512_maskz_mov_epi64(c1X | c2X, one); cY = _mm512_maskz_mov_epi64(c1Y | c2Y, one);
                    const __m512i RX = _mm512_or_si512(SX, _mm512_xor_si512(sX, Wm)), RY = _mm512_or_si512(SY, _mm512_xor_si512(sY, Wm));
                    __m512i dd = _mm512_add_epi64(_mm512_popcnt_epi64(_mm512_andnot_si512(RX, X)), _mm512_popcnt_epi64(_mm512_andnot_si512(RY, Y)));
                    dd = _mm512_add_epi64(dd, _mm512_popcnt_epi64(_mm512_and_si512(Z, neq)));
                    __m512i de = _mm512_add_epi64(_mm512_popcnt_epi64(Wm), _mm512_popcnt_epi64(_mm512_and_si512(X, RX)));
                    de = _mm512_add_epi64(de, _mm512_popcnt_epi64(_mm512_and_si512(Y, RY)));
                    diff = _mm512_add_epi64(diff, dd); dec = _mm512_add_epi64(dec, de);
                }
            }
        }
        _mm512_store_si512((__m512i*)dv, diff);
        _mm512_store_si512((__m512i*)mv, _mm512_sub_epi64(range, dec));
    }
    for (int l = 0; l < 8; l++) {
        if (!((ok >> l) & 1)) continue;
        if ((noov >> l) & 1) { res[l] = KR{1, 1, 1.0, 0}; continue; }
        if ((excl >> l) & 1) { res[l] = KR{0, 0, 0.0, 1}; continue; }
        long long mm = mv[l], d = dv[l];
        if (mm > 0) {
            double v = (double)d / (double)mm;
            if (v > cutoff) res[l] = KR{1, 1, 1.0, 0}; else res[l] = KR{(int)d, (int)mm, v, 0};
        } else res[l] = KR{1, 1, 1.0, 0};
    }
}
#endif

template <int C, bool EARLY>
static inline __attribute__((always_inline)) void rowsBody(const Job& J, long long i0, long long i1) {
    const Seqs& S = *J.S; const std::vector<Rec>& R = *J.recs; Out& out = *J.out;
    static thread_local FmtCache* fc = nullptr;
    if (!fc || fc->prec != (J.output == 0 ? 6 : 4) || fc->rows.size() != (size_t)S.L + 2) { delete fc; fc = new FmtCache(J.output == 0 ? 6 : 4, S.L); }
    const double cut = (J.output == 0) ? J.cutoff : 1.0;
    for (long long i = i0; i < i1; i++) {
        const std::string& ni = R[i].name;
        long long jend = (J.output == 2) ? S.n : i;
        if (J.output != 0) {
            if (ni.size() < 10) { std::string nm = ni; nm.resize(10, ' '); out.put(nm); } else out.put(ni);
            if (J.output == 2) out.putc_('\t');
        }
        for (long long j0 = 0; j0 < jend; j0 += 8) {
            int nl = (int)((jend - j0) < 8 ? (jend - j0) : 8);
            KR res[8]; unsigned done = 0, excl = 0;
#if ENG_ISA == 3
            if (C == ONEGAP_T || C == EACHGAP_T || C == NOGAPS) {
                if (J.output == 0) laneBlock<C, EARLY>(S, i, j0, nl, cut, res, done, excl);
                else laneBlock<C, false>(S, i, j0, nl, cut, res, done, excl);
            }
#endif
            unsigned todo = ((nl >= 8) ? 0xffu : ((1u << nl) - 1)) & ~excl;
            while (todo) {
                int l = __builtin_ctz(todo); todo &= todo - 1;
                long long j = j0 + l;
                char tmp[64]; const char* sl = nullptr; int len = 0;
                if (J.output == 0) {
                    KR r = ((done >> l) & 1) ? res[l] : kernel<C, EARLY>(S, i, j, cut);
                    if (r.kind == 1) continue;
                    if (!(r.v <= cut)) continue;
                    if (r.kind == 0) sl = fc->slot(r.diff, r.mml, r.v);
                    if (!sl) len = fmtg(tmp, r.v, 6);
                    out.put(ni); out.putc_(' '); out.put(R[j].name); out.putc_(' ');
                    if (sl) out.putSlot(sl); else out.put(tmp, len);
                    out.putc_('\n');
                } else {
                    if (j == i) { len = fmtg(tmp, 0.0, 4); }
                    else {
                        KR r = ((done >> l) & 1) ? res[l] : kernel<C, false>(S, i, j, 1.0);
                        if (r.kind == 0) sl = fc->slot(r.diff, r.mml, r.v);
                        if (!sl) len = fmtg(tmp, r.v, 4);
                    }
                    if (J.output == 1) out.putc_('\t');
                    if (sl) out.putSlot(sl); else out.put(tmp, len);
                    if (J.output == 2) out.putc_('\t');
                }
            }
        }
        if (J.output != 0) out.putc_('\n');
    }
}

template <int C, bool EARLY> void rowsFn(const Job& J, long long a, long long b) { rowsBody<C, EARLY>(J, a, b); }

template <int C> RowsFn pickC(bool early) { return early ? rowsFn<C, true> : rowsFn<C, false>; }
RowsFn pickRows(int calc, bool early) {
    switch (calc) {
        case ONEGAP_T: return pickC<ONEGAP_T>(early);
        case EACHGAP_T: return pickC<EACHGAP_T>(early);
        case NOGAPS: return pickC<NOGAPS>(early);
        case EACHGAP_F: return pickC<EACHGAP_F>(early);
        default: return pickC<ONEGAP_F>(early);
    }
}
