# swiftdistseqs

swiftdistseqs is a drop-in replacement for the [mothur](https://github.com/mothur/mothur) 1.48.6 command `dist.seqs` on
aligned DNA. It writes mothur's distance files, byte for byte, at any thread count, and does so much faster.

On three held-out data sets on a 2-vCPU AVX-512 machine it was 92-fold faster than mothur at one thread (57- to 172-fold per
set) and 93-fold faster at two. Without AVX-512, the AVX2 code path on the same machine gave 54-fold. Details are in the preprint (see
Citation).

## Build

Requirements: Linux on x86-64 and g++ with C++17.

```
solver/build.sh        # produces solver/bin/distseqs
```

The program picks an AVX-512, AVX2, POPCNT or baseline x86-64 code path at run time. All four write the same bytes.
`SWIFTDISTSEQS_ISA=avx2|popcnt|base` forces a slower one.

## Usage

```
solver/bin/distseqs fasta=final.fasta cutoff=0.03 processors=2     # writes final.dist, as mothur does
solver/bin/distseqs fasta=final.fasta output=lt                    # writes final.phylip.dist
solver/bin/distseqs --fasta in.fasta --out out.dist --calc eachgap --countends F --cutoff 0.1
```

- `calc` is onegap, eachgap or nogaps; `output` is column, lt or square; `countends` and `cutoff` work as in mothur.
- The output doesn't depend on the thread count and equals mothur's file with `processors=1`. mothur's own column output with
  two or more processors has the same lines in a different order.
- Not covered: protein calculators, `oldfasta`, `column` and `fitcalc` appending, `compress`, sequence names starting with `#`,
  and processors other than x86-64.

## Check it yourself

```
./prepare.sh           # builds mothur 1.48.6 from source and swiftdistseqs
./reproduce.sh         # byte comparison, timings and the application
```

`prepare.sh` needs make, git, g++ and zlib headers and took 14 minutes on one core of a laptop (`JOBS` sets the number of compile jobs).
`reproduce.sh` needs python3 with numpy for the application and takes about 20 minutes. `evaluate.py gate` prints `GATE PASS`
and exits 0 when every output matched mothur's, on fresh random inputs, on the 58 fixed inputs in `data/adversarial/` and on the
six real sets in `data/sets/`.

`data/README.txt` gives the sources of the data, and `application/` holds the analysis of OTU counts against the cutoff and its
results.

## Licence

GPL-3.0 (`LICENSE`), as mothur. swiftdistseqs reimplements the behaviour of `dist.seqs`. The SILVA-derived sets in `data/sets/`
are CC BY 4.0; see `data/README.txt`.

## Citation

See `CITATION.cff`. Please also cite mothur (Schloss et al., Appl Environ Microbiol 75:7537-7541, 2009).
