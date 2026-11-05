#!/usr/bin/env python3
"""Application: pairwise distances and OTU counts versus the
distance cutoff in the mothur MiSeq SOP data.

Runs, on data/sets/dev_sop_final.fasta (2,929 sequences) and its count table:
  1. SwiftDistSeqs at cutoff 1.0 (all pairs) and at cutoffs 0.01..0.10 (onegap, countends=T, column)
  2. mothur 1.48.6 cluster(method=opti) at each cutoff, after set.seed(seed=19760620)
  3. the same eleven dist.seqs runs with mothur 1.48.6 (timing, and a byte check against SwiftDistSeqs)
Writes app/results/*.tsv, app/results/summary.json. All runs use 1 thread (processors=1).
"""
import hashlib, json, os, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK = os.path.join(ROOT, "application", "work")
RES = os.path.join(ROOT, "application", "results")
SWIFT = os.path.join(ROOT, sys.argv[1] if len(sys.argv) > 1 else "solver", "bin", "distseqs")
REFW = os.path.join(ROOT, "ref", "bin", "distseqs")
MOTHUR = os.path.join(ROOT, "ref", "bin", "mothur")
FASTA = os.path.join(ROOT, "data", "sets", "dev_sop_final.fasta")
COUNT = os.path.join(ROOT, "data", "sop_run", "stability.trim.contigs.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table")
CUTS = ["0.01", "0.02", "0.03", "0.04", "0.05", "0.06", "0.07", "0.08", "0.09", "0.10"]


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run_dist(binary, cutoff, out):
    if os.path.exists(out):
        os.remove(out)
    t0 = time.perf_counter()
    subprocess.run([binary, "--fasta", FASTA, "--out", out, "--calc", "onegap", "--countends", "T",
                    "--cutoff", cutoff, "--output", "column", "--processors", "1"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return time.perf_counter() - t0


def main():
    shutil.rmtree(WORK, ignore_errors=True)
    os.makedirs(WORK)
    os.makedirs(RES, exist_ok=True)
    os.symlink(COUNT, os.path.join(WORK, "final.count_table"))
    rows = []
    # 1 + 3: distance files with both tools
    cut_list = ["1.0"] + CUTS
    for c in cut_list:
        so = os.path.join(WORK, "swift_%s.dist" % c)
        ro = os.path.join(WORK, "mothur_%s.dist" % c)
        ts = run_dist(SWIFT, c, so)
        tr = run_dist(REFW, c, ro)
        same = sha(so) == sha(ro)
        npairs = sum(1 for _ in open(so))
        rows.append({"cutoff": c, "swift_s": ts, "mothur_s": tr, "identical": same, "pairs_written": npairs, "bytes": os.path.getsize(so)})
        print("dist cutoff %s: swift %.3f s, mothur %.2f s, identical=%s, %d pairs" % (c, ts, tr, same, npairs), flush=True)
        os.remove(ro)
    # 2: opti clustering at each cutoff (mothur 1.48.6)
    clus = []
    for c in CUTS:
        src = os.path.join(WORK, "swift_%s.dist" % c)
        col = os.path.join(WORK, "c%s.dist" % c)
        os.replace(src, col)
        t0 = time.perf_counter()
        p = subprocess.run([MOTHUR, "#set.seed(seed=19760620); cluster(column=c%s.dist, count=final.count_table, method=opti, cutoff=%s)" % (c, c)],
                           cwd=WORK, capture_output=True, text=True)
        tc = time.perf_counter() - t0
        if p.returncode != 0:
            print(p.stdout[-2000:])
            sys.exit("cluster failed at %s" % c)
        base = "c%s.opti_mcc" % c
        lst = [f for f in os.listdir(WORK) if f.startswith(base) and f.endswith(".list")][0]
        ss = [f for f in os.listdir(WORK) if f.startswith(base) and f.endswith(".sensspec")][0]
        clus.append({"cutoff": c, "list": lst, "sensspec": ss, "cluster_s": tc})
        print("cluster %s: %.1f s (%s)" % (c, tc, lst), flush=True)
    # summaries
    counts = {}
    with open(COUNT) as f:  # mothur count table (compressed format: '#' comment lines, then a header)
        for line in f:
            if line.startswith("#") or line.startswith("Representative_Sequence"):
                continue
            x = line.rstrip("\n").split("\t")
            counts[x[0]] = int(x[1])
    total_reads = sum(counts.values())
    otu_rows = []
    for cl in clus:
        with open(os.path.join(WORK, cl["list"])) as f:
            f.readline()
            parts = f.readline().rstrip("\n").split("\t")
        label, notus, otus = parts[0], int(parts[1]), parts[2:]
        sizes = [sum(counts[n] for n in o.split(",")) for o in otus]
        nonsingle = sum(1 for s in sizes if s >= 2)
        single_seq = sum(1 for o in otus if "," not in o)
        with open(os.path.join(WORK, cl["sensspec"])) as f:
            h = f.readline().rstrip("\n").split("\t"); v = f.readline().rstrip("\n").split("\t")
        ssd = dict(zip(h, v))
        otu_rows.append({"cutoff": cl["cutoff"], "label": label, "otus": notus, "otus_ge2_reads": nonsingle,
                         "otus_one_sequence": single_seq, "mcc": float(ssd["mcc"]), "sensitivity": float(ssd["sensitivity"]),
                         "specificity": float(ssd["specificity"]), "cluster_s": cl["cluster_s"]})
    # distance distribution from the full file
    import numpy as np
    names = [l[1:].split()[0] for l in open(FASTA) if l.startswith(">")]
    idx = {n: k for k, n in enumerate(names)}
    nn = np.full(len(names), np.inf)
    import array
    dl = array.array("d")
    with open(os.path.join(WORK, "swift_1.0.dist")) as f:
        for line in f:
            a, b, d = line.split()
            d = float(d)
            dl.append(d)
            ia, ib = idx[a], idx[b]
            if d < nn[ia]: nn[ia] = d
            if d < nn[ib]: nn[ib] = d
    D = np.frombuffer(dl, dtype=np.float64)
    n_pairs = len(D)
    hist = np.histogram(D, bins=200, range=(0.0, 1.0))[0]
    le = {c: int(np.count_nonzero(D <= float(c))) for c in CUTS}
    with open(os.path.join(RES, "dist_hist.tsv"), "w") as f:
        f.write("bin_lo\tbin_hi\tpairs\n")
        for k in range(200):
            f.write("%.3f\t%.3f\t%d\n" % (k * 0.005, (k + 1) * 0.005, hist[k]))
    np.savetxt(os.path.join(RES, "nearest_neighbour.tsv"), nn, fmt="%.6g", header="nearest_neighbour_distance", comments="")
    with open(os.path.join(RES, "dist_runs.tsv"), "w") as f:
        f.write("cutoff\tswift_s\tmothur_s\tidentical\tpairs_written\tbytes\n")
        for r in rows:
            f.write("%s\t%.4f\t%.3f\t%s\t%d\t%d\n" % (r["cutoff"], r["swift_s"], r["mothur_s"], r["identical"], r["pairs_written"], r["bytes"]))
    with open(os.path.join(RES, "otus_by_cutoff.tsv"), "w") as f:
        f.write("cutoff\totus\totus_ge2_reads\totus_one_sequence\tmcc\tsensitivity\tspecificity\tpairs_le_cutoff\tfrac_pairs_le_cutoff\tfrac_seqs_nn_le_cutoff\tcluster_s\n")
        for o in otu_rows:
            c = o["cutoff"]
            f.write("%s\t%d\t%d\t%d\t%.6f\t%.6f\t%.6f\t%d\t%.6g\t%.6g\t%.2f\n" % (c, o["otus"], o["otus_ge2_reads"], o["otus_one_sequence"], o["mcc"],
                    o["sensitivity"], o["specificity"], le[c], le[c] / n_pairs, float(np.mean(nn <= float(c))), o["cluster_s"]))
    summary = {
        "n_sequences": len(names), "n_pairs": n_pairs, "total_reads": total_reads,
        "swift_total_s": sum(r["swift_s"] for r in rows), "mothur_total_s": sum(r["mothur_s"] for r in rows),
        "all_identical": all(r["identical"] for r in rows),
        "cluster_total_s": sum(o["cluster_s"] for o in otu_rows),
        "median_pair_distance": float(np.median(D)), "mean_pair_distance": float(np.mean(D)),
        "frac_pairs_le_0.03": le["0.03"] / n_pairs,
        "nn_median": float(np.median(nn)),
        "when": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    json.dump(summary, open(os.path.join(RES, "summary.json"), "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
