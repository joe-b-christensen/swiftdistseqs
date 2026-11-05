#!/usr/bin/env python3
"""SwiftDistSeqs correctness gate and timing harness.

A build is a directory whose bin/distseqs takes
    --fasta IN --out OUT --calc C --countends E --cutoff X --output O --processors P
and must write exactly the bytes that mothur 1.48.6 dist.seqs writes (same exit status,
same presence or absence of the output file). The live reference is ref/bin/distseqs
(mothur 1.48.6, processors=1).

  evaluate.py gate BUILD [--fuzz N] [--seed S] [--held]   correctness gate
  evaluate.py time BUILD [--sets dev|held|all] [--threads T] [--reps K] [--window W]
                                                        speed-up of BUILD over the reference

Timing: every sample is one single run (one process). The output file is deleted before
each repetition, outside the timer. Repetitions stop at K samples or when the fixed window
of W seconds for that (build, set) is used up, whichever comes first (at least one sample).
"""
import argparse, hashlib, json, math, os, statistics, subprocess, sys, tempfile, time

ROOT = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(ROOT, "ref")
LOGDIR = os.path.join(ROOT, "logs", "eval")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()



def load_sets(role):
    sets = json.load(open(os.path.join(ROOT, "data", "sets", "sets.json")))
    if role == "all":
        return sets
    return [s for s in sets if s["role"] == role]


def run(build_bin, case, out, processors):
    if os.path.exists(out):
        os.remove(out)
    cmd = [build_bin, "--fasta", case["fasta"], "--out", out, "--calc", case["calc"],
           "--countends", case["countends"], "--cutoff", case["cutoff"], "--output", case["output"],
           "--processors", str(processors)]
    t0 = time.perf_counter()
    p = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    t1 = time.perf_counter()
    return p.returncode, t1 - t0


def bin_of(build):
    b = os.path.join(os.path.abspath(build), "bin", "distseqs")
    if not os.access(b, os.X_OK):
        sys.exit("no executable %s" % b)
    return b


def compare(case, cand_bin, tmp, cand_procs):
    ro, co = os.path.join(tmp, "ref.out"), os.path.join(tmp, "cand.out")
    rrc, rt = run(os.path.join(REF, "bin", "distseqs"), case, ro, 1)
    crc, ct = run(cand_bin, case, co, cand_procs)
    rex, cex = os.path.exists(ro), os.path.exists(co)
    ok = (rrc == crc) and (rex == cex)
    detail = ""
    if ok and rex:
        ok = os.path.getsize(ro) == os.path.getsize(co) and sha256(ro) == sha256(co)
        if not ok:
            detail = "bytes differ (ref %d B, cand %d B)" % (os.path.getsize(ro), os.path.getsize(co))
    elif not ok:
        detail = "exit ref=%d cand=%d; file ref=%s cand=%s" % (rrc, crc, rex, cex)
    size = os.path.getsize(ro) if rex else -1
    for f in (ro, co):
        if os.path.exists(f):
            os.remove(f)
    return ok, detail, rrc, size, rt, ct


def cmd_gate(a):
    cand = bin_of(a.build)
    seed = a.seed if a.seed is not None else int(time.time_ns() % 2_000_000_000)
    rec = {"build": os.path.abspath(a.build), "seed": seed, "start": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "fails": []}
    with tempfile.TemporaryDirectory(prefix="gate.") as tmp:
        fz = os.path.join(tmp, "fuzz")
        subprocess.run([sys.executable, "-I", os.path.join(ROOT, "gen_data.py"), "fuzz", "--seed", str(seed), "--n", str(a.fuzz), "--outdir", fz], check=True, stdout=subprocess.DEVNULL)
        cases = json.load(open(os.path.join(fz, "cases.json")))["cases"]
        fixed = os.path.join(ROOT, "data", "adversarial", "cases.json")
        if os.path.exists(fixed):  # fixed adversarial suite, run in addition to the fresh cases
            for c in json.load(open(fixed))["cases"]:
                c = dict(c)
                c["fasta"] = os.path.join(ROOT, c["fasta"])
                cases.append(c)
        nf, rcs = 0, {}
        for c in cases:
            ok, detail, rrc, size, rt, ct = compare(c, cand, tmp, c["processors"])
            rcs[rrc] = rcs.get(rrc, 0) + 1
            if not ok:
                nf += 1
                rec["fails"].append({"case": c["id"], "opts": {k: c[k] for k in ("calc", "countends", "cutoff", "output", "processors")}, "detail": detail})
                if a.keep_failures:
                    os.makedirs(a.keep_failures, exist_ok=True)
                    import shutil
                    shutil.copy(c["fasta"], os.path.join(a.keep_failures, c["id"] + ".fasta"))
                    json.dump(c, open(os.path.join(a.keep_failures, c["id"] + ".json"), "w"))
        rec["fuzz_cases"] = len(cases)
        rec["fuzz_fail"] = nf
        rec["ref_exit_codes"] = rcs
        print("fuzz: %d cases (seed %d), %d failed; reference exit codes %s" % (len(cases), seed, nf, rcs))
        roles = ["dev"] + (["held"] if a.held else [])
        for role in roles:
            for s in load_sets(role):
                case = dict(s)
                case["fasta"] = os.path.join(ROOT, s["fasta"])
                ro = os.path.join(tmp, "ref_set.out")
                rrc, rt = run(os.path.join(REF, "bin", "distseqs"), case, ro, 1)  # live reference, once per set
                rex = os.path.exists(ro)
                rsize = os.path.getsize(ro) if rex else -1
                rhash = sha256(ro) if rex else None
                for p in sorted(set([1, a.threads])):
                    co = os.path.join(tmp, "cand_set.out")
                    crc, ct = run(cand, case, co, p)
                    cex = os.path.exists(co)
                    ok = (rrc == crc) and (rex == cex) and (not rex or (os.path.getsize(co) == rsize and sha256(co) == rhash))
                    detail = "" if ok else "exit ref=%d cand=%d; file ref=%s cand=%s; size ref=%d cand=%d" % (
                        rrc, crc, rex, cex, rsize, os.path.getsize(co) if cex else -1)
                    if cex:
                        os.remove(co)
                    print("%-5s %-18s procs=%d %s  ref %.2fs cand %.3fs  (%d bytes) %s" % (role, s["id"], p, "PASS" if ok else "FAIL", rt, ct, rsize, detail))
                    if not ok:
                        rec["fails"].append({"case": s["id"], "procs": p, "detail": detail})
                if rex:
                    os.remove(ro)
    rec["pass"] = len(rec["fails"]) == 0
    os.makedirs(LOGDIR, exist_ok=True)
    fn = os.path.join(LOGDIR, "gate_%s.json" % time.strftime("%Y%m%dT%H%M%S", time.gmtime()))
    json.dump(rec, open(fn, "w"), indent=1)
    print("GATE %s  (%s)" % ("PASS" if rec["pass"] else "FAIL", fn))
    sys.exit(0 if rec["pass"] else 1)


def measure(binpath, case, procs, reps, window, tmp):
    out = os.path.join(tmp, "t.out")
    times, used = [], 0.0
    while len(times) < reps and (not times or used < window):
        rc, t = run(binpath, case, out, procs)
        if rc != 0:
            sys.exit("timing run failed (exit %d) for %s" % (rc, binpath))
        times.append(t)
        used += t
    if os.path.exists(out):
        os.remove(out)
    return times


def gmean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def cmd_time(a):
    cand = bin_of(a.build)
    refb = os.path.join(REF, "bin", "distseqs") if not a.vs else bin_of(a.vs)
    rows = []
    with tempfile.TemporaryDirectory(prefix="time.") as tmp:
        for s in load_sets(a.sets):
            case = dict(s)
            case["fasta"] = os.path.join(ROOT, s["fasta"])
            rt = measure(refb, case, a.threads, a.reps, a.window, tmp)
            ct = measure(cand, case, a.threads, a.reps, a.window, tmp)
            row = {"set": s["id"], "role": s["role"], "threads": a.threads, "ref_times": rt, "cand_times": ct,
                   "ref_best": min(rt), "cand_best": min(ct), "ref_median": statistics.median(rt), "cand_median": statistics.median(ct)}
            row["speedup_best"] = row["ref_best"] / row["cand_best"]
            row["speedup_median"] = row["ref_median"] / row["cand_median"]
            rows.append(row)
            print("%-18s T=%d ref best %.3fs med %.3fs | cand best %.4fs med %.4fs | speed-up best %.1fx median %.1fx (n=%d/%d)" % (
                s["id"], a.threads, row["ref_best"], row["ref_median"], row["cand_best"], row["cand_median"], row["speedup_best"], row["speedup_median"], len(rt), len(ct)))
    res = {"build": os.path.abspath(a.build), "reference": refb, "sets": a.sets, "threads": a.threads,
           "reps": a.reps, "window": a.window, "rows": rows,
           "gmean_speedup_best": gmean([r["speedup_best"] for r in rows]),
           "gmean_speedup_median": gmean([r["speedup_median"] for r in rows]),
           "when": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    print("geometric mean speed-up over %s sets (T=%d): best-of-%d %.2fx, median %.2fx" % (a.sets, a.threads, a.reps, res["gmean_speedup_best"], res["gmean_speedup_median"]))
    os.makedirs(LOGDIR, exist_ok=True)
    fn = a.json or os.path.join(LOGDIR, "time_%s_%s_T%d.json" % (a.sets, time.strftime("%Y%m%dT%H%M%S", time.gmtime()), a.threads))
    json.dump(res, open(fn, "w"), indent=1)
    print("wrote", fn)


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    g = sp.add_parser("gate")
    g.add_argument("build")
    g.add_argument("--fuzz", type=int, default=60)
    g.add_argument("--seed", type=int)
    g.add_argument("--held", action="store_true")
    g.add_argument("--threads", type=int, default=2)
    g.add_argument("--keep-failures")
    t = sp.add_parser("time")
    t.add_argument("build")
    t.add_argument("--vs", help="compare against this build instead of the reference")
    t.add_argument("--sets", default="dev", choices=["dev", "held", "all"])
    t.add_argument("--threads", type=int, default=1)
    t.add_argument("--reps", type=int, default=3)
    t.add_argument("--window", type=float, default=300.0)
    t.add_argument("--json")
    a = ap.parse_args()
    {"gate": cmd_gate, "time": cmd_time}[a.cmd](a)


if __name__ == "__main__":
    main()
