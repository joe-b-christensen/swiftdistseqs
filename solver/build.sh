#!/bin/bash
# Build SwiftDistSeqs into solver/bin/distseqs (portable x86-64 baseline; CPU-specific kernels are dispatched at run time).
set -e
cd "$(dirname "$0")"
mkdir -p bin
${CXX:-g++} -O3 -std=c++17 -pthread -Wall -Wextra -o bin/distseqs src/swiftdistseqs.cpp
