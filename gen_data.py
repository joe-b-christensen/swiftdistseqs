#!/usr/bin/env python3
"""SwiftDistSeqs data generator.

  gen_data.py fuzz --seed S --n N --outdir DIR
      Fresh random adversarial cases for the correctness gate. Each case is one FASTA
      file plus the dist.seqs options to run it with (calc, countends, cutoff, output)
      and a thread count for the candidate. Writes DIR/cases.json.

  gen_data.py sets
      Builds the real development and held-out sets from the processed public data
      (deterministic subsets, fixed seeds) into data/sets/ and writes data/sets/sets.json.
"""
import argparse, json, os, random, sys, hashlib

ROOT = os.path.dirname(os.path.abspath(__file__))

CALCS = ["onegap", "eachgap", "nogaps"]
COUNTENDS = ["T", "F", "t", "f", "true", "false", "TRUE", "FALSE"]
CUTOFFS = ["0", "0.0", "0.01", "0.02", "0.03", "0.04", "0.05", "0.1", "0.10", "0.15", "0.2",
           "0.25", "0.3", "0.333333", "0.5", "0.6", "0.75", "0.9", "0.99", "1", "1.0", "1.5", "2"]
OUTPUTS = ["column", "column", "column", "lt", "lt", "square", "phylip"]
IUPAC = "RYKMSWBDHVN"
ODD = "*~?0123456789xXzZjJ"
NAMECH = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.|-"


def rand_name(r):
    n = r.choice([1, 2, 3, 5, 8, 9, 10, 11, 12, 15, 20, 25])
    s = "".join(r.choice(NAMECH) for _ in range(n))
    if r.random() < 0.08:  # mothur rewrites ':' to '_'
        p = r.randrange(len(s) + 1)
        s = s[:p] + ":" + s[p:]
    if s[0] == "#":
        s = "x" + s
    return s


def mutate(r, anc, sub, gap, ambig, lower, odd, internal_dot):
    s = list(anc)
    L = len(s)
    for i in range(L):
        x = r.random()
        if x < sub:
            s[i] = r.choice("ACGT")
        elif x < sub + gap:
            s[i] = "-"
    # gap runs (indels)
    for _ in range(r.choice([0, 0, 1, 2, 4])):
        if L == 0:
            break
        a = r.randrange(L)
        b = min(L, a + r.choice([1, 2, 3, 5, 9, 20]))
        for i in range(a, b):
            s[i] = "-"
    for i in range(L):
        y = r.random()
        if y < ambig:
            s[i] = r.choice(IUPAC + "U")
        elif y < ambig + odd:
            s[i] = r.choice(ODD)
        elif y < ambig + odd + internal_dot:
            s[i] = "."
    # terminal padding: leading/trailing '.' or '-'
    if L > 0 and r.random() < 0.7:
        lead = r.choice([0, 0, 1, 3, L // 10, L // 4, L // 2, L])
        trail = r.choice([0, 0, 1, 3, L // 10, L // 4, L // 2, L])
        lc = r.choice([".", "-", "-", "."])
        tc = r.choice([".", "-", lc, lc])
        for i in range(min(lead, L)):
            s[i] = lc
        for i in range(max(0, L - trail), L):
            s[i] = tc
    if lower and L:
        s = [c.lower() if r.random() < lower else c for c in s]
    return "".join(s)


def make_case(r):
    kind = r.random()
    if kind < 0.03:
        n = r.choice([0, 1])
    elif kind < 0.10:
        n = r.randint(120, 320)
    else:
        n = r.randint(2, 60)
    L = r.choice([1, 2, 3, 5, 8, 16, 31, 32, 33, 50, 63, 64, 65, 100, 127, 128, 129, 200, 255,
                  256, 257, 300, 383, 384, 385, 450, 511, 512, 513, 700, 1024, 1500])
    if n > 100:
        L = min(L, 400)
    anc = "".join(r.choice("ACGT") for _ in range(L))
    # shared alignment gap columns (as in real alignments)
    shared = set(r.sample(range(L), k=min(L, r.choice([0, 0, L // 20, L // 5, L // 2])))) if L else set()
    sub = r.choice([0.0, 0.002, 0.01, 0.03, 0.08, 0.2, 0.5])
    gap = r.choice([0.0, 0.0, 0.005, 0.02, 0.1, 0.3])
    ambig = r.choice([0.0, 0.0, 0.0, 0.001, 0.01, 0.05])
    odd = r.choice([0.0, 0.0, 0.0, 0.0, 0.002, 0.02])
    idot = r.choice([0.0, 0.0, 0.0, 0.0, 0.003, 0.03])
    lower = r.choice([0.0, 0.0, 0.0, 0.1, 1.0])
    clusters = [anc] + ["".join(c if r.random() > 0.15 else r.choice("ACGT") for c in anc)
                        for _ in range(r.choice([0, 1, 3]))]
    seqs = []
    for k in range(n):
        base = r.choice(clusters)
        if seqs and r.random() < 0.05:
            s = r.choice(seqs)[1]  # exact duplicate sequence
        elif r.random() < 0.02:
            s = r.choice([".", "-", "N"]) * L  # all-gap / all-N rows
        else:
            s = mutate(r, base, sub, gap, ambig, lower, odd, idot)
            if shared:
                s = "".join("-" if (i in shared and s[i] not in ".") else s[i] for i in range(L))
        seqs.append((rand_name(r), s))
    if n >= 2 and r.random() < 0.03:  # unequal lengths -> mothur aborts
        j = r.randrange(n)
        nm, s = seqs[j]
        seqs[j] = (nm, s + r.choice("ACGT-."))
    # serialise with formatting noise
    eol = "\r\n" if r.random() < 0.15 else "\n"
    wrap = r.choice([0, 0, 0, 60, 80, 7])
    parts = []
    for nm, s in seqs:
        hdr = ">" + nm
        if r.random() < 0.2:
            hdr += r.choice([" ", "\t"]) + r.choice(["comment", "x=1 y=2", "len=%d" % len(s), ""])
        parts.append(hdr + eol)
        if wrap and len(s) > wrap:
            lines = [s[i:i + wrap] for i in range(0, len(s), wrap)]
        else:
            lines = [s]
        for ln in lines:
            pre = " " * r.choice([0, 0, 0, 0, 0, 1, 2])
            post = r.choice(["", "", "", "", " ", "\t", "  "])
            parts.append(pre + ln + post + eol)
        if r.random() < 0.1:
            parts.append(eol)
    text = "".join(parts)
    if text and r.random() < 0.1:
        text = text.rstrip("\r\n")  # no final newline
    opts = {
        "calc": r.choice(CALCS),
        "countends": r.choice(COUNTENDS),
        "cutoff": r.choice(CUTOFFS) if r.random() < 0.9 else "%.4f" % r.random(),
        "output": r.choice(OUTPUTS),
        "processors": r.choice([1, 1, 2, 2, 3, 4]),
    }
    return text, opts


def cmd_fuzz(a):
    r = random.Random(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    cases = []
    for i in range(a.n):
        text, opts = make_case(r)
        fn = os.path.join(a.outdir, "case_%04d.fasta" % i)
        with open(fn, "w", newline="") as f:
            f.write(text)
        opts["fasta"] = fn
        opts["id"] = "fuzz_%d_%04d" % (a.seed, i)
        cases.append(opts)
    # make sure every option value of every output-affecting option is present
    combos = [(c, e, o) for c in CALCS for e in ["T", "F"] for o in ["column", "lt", "square"]]
    for j, (c, e, o) in enumerate(combos):
        text, opts = make_case(r)
        fn = os.path.join(a.outdir, "combo_%02d.fasta" % j)
        with open(fn, "w", newline="") as f:
            f.write(text)
        opts.update({"calc": c, "countends": e, "output": o, "fasta": fn, "id": "fuzz_%d_combo_%02d" % (a.seed, j)})
        cases.append(opts)
    json.dump({"seed": a.seed, "cases": cases}, open(os.path.join(a.outdir, "cases.json"), "w"), indent=1)
    print("wrote %d cases to %s" % (len(cases), a.outdir))


# ---------------------------------------------------------------- real sets
def read_fasta(fn):
    recs, name, buf = [], None, []
    with open(fn) as f:
        for line in f:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def write_fasta(fn, recs):
    with open(fn, "w") as f:
        for n, s in recs:
            f.write(">%s\n%s\n" % (n, s))


def filter_vertical_trump(recs):
    """Same column filter as mothur filter.seqs(vertical=T, trump=.): drop columns that are
    a gap ('-' or '.') in every sequence or contain '.' in any sequence."""
    L = len(recs[0][1])
    keep = []
    for i in range(L):
        col = [s[i] for _, s in recs]
        if any(c == "." for c in col):
            continue
        if all(c in "-." for c in col):
            continue
        keep.append(i)
    return [(n, "".join(s[i] for i in keep)) for n, s in recs]


def cmd_sets(a):
    out = os.path.join(ROOT, "data", "sets")
    os.makedirs(out, exist_ok=True)
    sop = os.path.join(ROOT, "data", "sop_run")
    uniq = read_fasta(os.path.join(sop, "stability.trim.contigs.unique.good.filter.unique.fasta"))
    final = read_fasta(os.path.join(sop, "stability.trim.contigs.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.fasta"))
    v4 = read_fasta(os.path.join(sop, "silva.v4.fasta"))
    full = read_fasta(os.path.join(ROOT, "data", "raw", "silva", "silva.seed_v138_2.align"))
    tax = {}
    for line in open(os.path.join(ROOT, "data", "raw", "silva", "silva.seed_v138_2.tax")):
        k, v = line.rstrip("\n").split("\t", 1)
        tax[k] = v
    sets = []

    # SOP unique sequences before pre.cluster: two disjoint random subsets (dev / held-out)
    r = random.Random(20261104)
    idx = list(range(len(uniq)))
    r.shuffle(idx)
    n_sop = a.n_sop
    dev_sop = [uniq[i] for i in sorted(idx[:n_sop])]
    held_sop = [uniq[i] for i in sorted(idx[n_sop:2 * n_sop])]
    write_fasta(os.path.join(out, "dev_sop_uniq.fasta"), dev_sop)
    write_fasta(os.path.join(out, "held_sop_uniq.fasta"), held_sop)
    sets.append({"id": "dev_sop_uniq", "role": "dev", "fasta": "data/sets/dev_sop_uniq.fasta",
                 "calc": "onegap", "countends": "T", "cutoff": "0.03", "output": "column",
                 "desc": "MiSeq SOP V4 unique sequences before pre.cluster, random %d; SOP-style dist.seqs(cutoff=0.03)" % n_sop})
    sets.append({"id": "held_sop_uniq", "role": "held", "fasta": "data/sets/held_sop_uniq.fasta",
                 "calc": "onegap", "countends": "T", "cutoff": "0.03", "output": "column",
                 "desc": "MiSeq SOP V4 unique sequences, disjoint random %d; dist.seqs(cutoff=0.03)" % n_sop})

    # SOP final sequences (after pre.cluster, chimera removal, lineage and Mock removal): lt for clearcut
    write_fasta(os.path.join(out, "dev_sop_final.fasta"), final)
    sets.append({"id": "dev_sop_final_lt", "role": "dev", "fasta": "data/sets/dev_sop_final.fasta",
                 "calc": "onegap", "countends": "T", "cutoff": "1.0", "output": "lt",
                 "desc": "MiSeq SOP final %d sequences; the SOP's dist.seqs(output=lt) for clearcut" % len(final)})

    # SILVA seed full-length bacteria: disjoint subsets, each filtered (vertical=T, trump=.)
    bact = [x for x in full if tax.get(x[0], "").startswith("Bacteria;")]
    r = random.Random(138)
    idx = list(range(len(bact)))
    r.shuffle(idx)
    n_full = a.n_full
    dev_full = filter_vertical_trump([bact[i] for i in sorted(idx[:n_full])])
    held_full = filter_vertical_trump([bact[i] for i in sorted(idx[n_full:2 * n_full])])
    write_fasta(os.path.join(out, "dev_silva_full.fasta"), dev_full)
    write_fasta(os.path.join(out, "held_silva_full.fasta"), held_full)
    sets.append({"id": "dev_silva_full", "role": "dev", "fasta": "data/sets/dev_silva_full.fasta",
                 "calc": "onegap", "countends": "T", "cutoff": "0.03", "output": "column",
                 "desc": "SILVA SEED v138.2 full-length Bacteria, random %d, filtered to %d columns; cutoff 0.03" % (n_full, len(dev_full[0][1]))})
    sets.append({"id": "held_silva_full", "role": "held", "fasta": "data/sets/held_silva_full.fasta",
                 "calc": "nogaps", "countends": "T", "cutoff": "0.10", "output": "column",
                 "desc": "SILVA SEED v138.2 full-length Bacteria, disjoint random %d, filtered to %d columns; calc=nogaps, cutoff 0.10" % (n_full, len(held_full[0][1]))})

    # SILVA seed V4 region (pcr.seqs output), held-out: lt with eachgap, countends=F
    v4b = [x for x in v4 if tax.get(x[0], "").startswith("Bacteria;")]
    r = random.Random(4)
    idx = list(range(len(v4b)))
    r.shuffle(idx)
    held_v4 = filter_vertical_trump([v4b[i] for i in sorted(idx[:a.n_v4])])
    write_fasta(os.path.join(out, "held_silva_v4.fasta"), held_v4)
    sets.append({"id": "held_silva_v4_lt", "role": "held", "fasta": "data/sets/held_silva_v4.fasta",
                 "calc": "eachgap", "countends": "F", "cutoff": "1.0", "output": "lt",
                 "desc": "SILVA SEED v138.2 V4 region (pcr.seqs 11894-25319), Bacteria, random %d, filtered to %d columns; calc=eachgap, countends=F, output=lt" % (a.n_v4, len(held_v4[0][1]))})

    for s in sets:
        p = os.path.join(ROOT, s["fasta"])
        s["bytes"] = os.path.getsize(p)
        s["sha256"] = hashlib.sha256(open(p, "rb").read()).hexdigest()
        s["nseqs"] = sum(1 for line in open(p) if line.startswith(">"))
    json.dump(sets, open(os.path.join(out, "sets.json"), "w"), indent=1)
    for s in sets:
        print(s["id"], s["role"], s["nseqs"], s["bytes"])


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    f = sp.add_parser("fuzz")
    f.add_argument("--seed", type=int, required=True)
    f.add_argument("--n", type=int, default=60)
    f.add_argument("--outdir", required=True)
    s = sp.add_parser("sets")
    s.add_argument("--n-sop", type=int, default=6000)
    s.add_argument("--n-full", type=int, default=2000)
    s.add_argument("--n-v4", type=int, default=2500)
    a = ap.parse_args()
    {"fuzz": cmd_fuzz, "sets": cmd_sets}[a.cmd](a)


if __name__ == "__main__":
    main()
