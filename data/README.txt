Data sets used for the comparison and timings.

sets/        Six alignments, with options, in sets.json. Three development sets and three held-out sets, fixed before any
             timing of SwiftDistSeqs.
adversarial/ 58 fixed random and edge-case inputs, with options in cases.json. gen_data.py writes fresh ones for every run of
             evaluate.py gate.
sop_run/     sop.batch, the mothur 1.48.6 commands that processed the MiSeq SOP reads, and the count table of the final
             2,929 sequences.

Sources, downloaded on 4 November 2026:
- MiSeq SOP data (Kozich et al. 2013, Schloss et al. 2012):
  https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
  36,879,679 bytes, sha256 728f7544824f8ca0223618b18bb60fce6a5dc20ae3d476267dbd5cf5594a35d4
- SILVA SEED v138.2, as packaged by mothur (Quast et al. 2013):
  https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.seed_v138_2.tgz
  15,358,615 bytes, sha256 2f62d070ae67a4121c4597299c76cc43534b050f8f0de77c5cdb2632991fdf1f

The SILVA-derived files (sets/dev_silva_full.fasta, sets/held_silva_full.fasta, sets/held_silva_v4.fasta) are subsets of
SILVA release 138.2, (c) the SILVA authors, used under CC BY 4.0: https://creativecommons.org/licenses/by/4.0/
They were sampled at random and had their empty alignment columns removed; see https://www.arb-silva.de/silva-license-information

"gen_data.py sets" rebuilds sets/ from the downloads after sop.batch has been run in data/sop_run/ with the SILVA files
alongside. That path was not rerun for this release.
