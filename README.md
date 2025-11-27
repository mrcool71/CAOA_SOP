# CAOA Toolkit — Multithreaded Strength-3 and Strength-4

## Overview
We search and verify binary Cyclically-Arranged Orthogonal Arrays (CAOA) with two approaches: a De Bruijn/Eulerian constructor for strength-3 and a residue-statistics search for strength-3 and strength-4. This README describes the multithreaded C++ implementations and how to run them.

## Files
- `strength3_mt.cpp` — Multithreaded strength-3 residue search (OpenMP). Writes `generation_results_mt.txt`.
- `strength4_mt.cpp` — Multithreaded strength-4 search with bandwidth tolerance (tol=1) and Theorem-2 screen. Writes `results_strength4_bw1_mt.txt`.
- Reference files in this folder:
  - `generation-verification-version0.py` (Python constructor + verifier for strength-3)
  - `generation-version1.cpp`, `generation-version2.cpp` (earlier strength-3 residue searches)
  - `strength4-tentative.cpp` (exact strength-4, single n)
  - `strength4-v2.cpp` (serial strength-4 with bandwidth tolerance)
- Dataset: `CAOA_vectors.xlsx` (strength-3/4 vectors and timings).

## Build
Requires a compiler with OpenMP (e.g., GCC/Clang).
```bash
g++ -O2 -fopenmp strength3_mt.cpp -o strength3_mt
g++ -O2 -fopenmp strength4_mt.cpp -o strength4_mt
```

## Run
```bash
./strength3_mt
# outputs generation_results_mt.txt for n = 1..10

./strength4_mt
# outputs results_strength4_bw1_mt.txt for n = 8..18
```

## What the MT versions do
- Enumerate balanced partitions with element 1 fixed in V0.
- Parallelize per-partition evaluation with OpenMP (`schedule(dynamic,64)`).
- Track the best strength k and one best partition.
- For strength-4, also record the first partition that passes the Theorem-2 modular screen under bandwidth tolerance.

## Interpreting outputs
- `generation_results_mt.txt`: columns `n`, `k`, generating vector (0 marks V0), runtime in ms.
- `results_strength4_bw1_mt.txt`: for each n, runtime, best k by ID logic, best partition (V0/V1), and the first theorem-valid partition if found.

## Notes and options
- Bandwidth tolerance in strength-4 is set to 1; adjust `BANDWIDTH_TOL` in `strength4_mt.cpp` if needed.
- For exact equality checks, use `strength4-tentative.cpp` (single n) or the serial `strength4-v2.cpp` with tolerance set to 0.
- For larger n, increase threads via `OMP_NUM_THREADS`.

## Suggested workflow
1. For quick strength-3 baselines: run `strength3_mt`.
2. For strength-4 exploration with tolerance: run `strength4_mt` and inspect `results_strength4_bw1_mt.txt`.
3. Cross-check any promising vector with the Python verifier (`generation-verification-version0.py`) or the exact strength-4 checker when you need strict equality.

## Future work
We plan to push to larger n by keeping bandwidth/modular prefilters, adding deeper parallelism (per-residue counting and per-partition batching), and possibly adding a SAT/SMT front-end with multithreaded exact verification.
