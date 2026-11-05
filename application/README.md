# Application

Pairwise distances and OTU counts against the distance cutoff in the mothur MiSeq SOP data: 2,929 final sequences from
19 fecal samples of one mouse (`data/sets/dev_sop_final.fasta` and the count table in `data/sop_run/`).

`scripts/run_application.py` runs SwiftDistSeqs and mothur at a cutoff of 1.0 and at 0.01 to 0.10 (onegap, countends=T,
column, one thread), checks that both wrote the same bytes, and clusters each file with mothur's `cluster(method=opti)`
after `set.seed(seed=19760620)`. `reproduce.sh` runs it after the comparison and timings. It uses about 1 GB of disk for
the distance files in `application/work/`.

`results/` holds the published tables:

- `otus_by_cutoff.tsv`: OTUs, OTUs with at least two reads, OTUs of one sequence, opti MCC, and pairs and sequences within
  each cutoff.
- `dist_runs.tsv`: run times, byte check, pairs written and file size for each cutoff.
- `dist_hist.tsv`: histogram of all 4,288,056 distances in bins of 0.005.
- `nearest_neighbour.tsv`: each sequence's nearest-neighbour distance.
- `summary.json`: totals.
