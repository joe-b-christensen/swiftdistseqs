#!/bin/bash
# Compares SwiftDistSeqs with mothur byte for byte, times both, and reruns the application.
# Needs python3; the application step also needs numpy.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p results
python3 -I evaluate.py gate solver --held
for s in dev held; do for t in 1 2; do
  python3 -I evaluate.py time solver --sets "$s" --threads "$t" --reps 3 --json "results/time_${s}_T${t}.json"
done; done
python3 application/scripts/run_application.py solver
