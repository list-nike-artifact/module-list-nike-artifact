# `CORAL`, `log p = 509`

The row of the comparison table.  CORAL is benchmarked through its own harness,
with only the medians extracted by `median.py`, so nothing about the measurement
is ours.  `coral-c/` is upstream's reference implementation, vendored unmodified,
and `sh fetch.sh` diffs it against CORAL-nike/CORAL at the revision we measured.
The repository's 368 MB Sage specification is not vendored, since it is not what
we measured.

    cd coral-c
    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    build/src/nike/ref/p_500/test/test_nike_p_500 --single --iterations=1000 \
        | tr '\r' '\n' | python3 ../median.py 4.05

`p_500` is the 509-bit prime (`Nbits 509` in `src/gf/ref/p_500/fp_p_500_64.c`),
which is CORAL's parameter set at the security level of ours; the larger primes
target higher levels and are not comparable.  `tr` is needed because the harness
overwrites one progress line with carriage returns.

On the machine of the paper this gives 32.30 million cycles for key generation
and 15.11 million for a derivation.
