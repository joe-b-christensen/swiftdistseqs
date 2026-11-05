#!/bin/bash
# Builds the reference (mothur 1.48.6, from source) and SwiftDistSeqs.
# Needs Linux on x86-64, g++ 9 or newer, make, git and zlib headers. JOBS sets the number of compile jobs (default 2).
set -euo pipefail
cd "$(dirname "$0")"
J=${JOBS:-2}
mkdir -p logs
if [ ! -x ref/bin/mothur ]; then
  [ -d ref/mothur_src ] || git clone --quiet https://github.com/mothur/mothur.git ref/mothur_src
  (cd ref/mothur_src && git checkout --quiet v1.48.6 && git apply ../patches/build_makefile.patch \
     && make -j"$J" USEREADLINE=no CPPFLAGS="-include memory" > ../../logs/mothur_build.log 2>&1)
  cp ref/mothur_src/mothur ref/bin/
fi
solver/build.sh
echo "prepared: ref/bin/mothur, solver/bin/distseqs"
