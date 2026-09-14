# `ML-KEM-768`, reference C

The row of the comparison table.  `bench.c` is the harness; `ref/` is upstream's
reference implementation, vendored unmodified, and `sh fetch.sh` diffs it against
pq-crystals/kyber at the revision we measured.

    cc -O3 -fomit-frame-pointer -DKYBER_K=3 -Iref -o bench \
       ref/kem.c ref/indcpa.c ref/polyvec.c ref/poly.c ref/ntt.c ref/cbd.c \
       ref/reduce.c ref/verify.c ref/fips202.c ref/symmetric-shake.c \
       ref/randombytes.c bench.c
    ./bench 4.05

The harness is compiled from here rather than copied into `ref/`, so that `ref/`
stays byte-identical to upstream and `sh fetch.sh` keeps passing.

The compiler flags are upstream's own, minus `-z noexecstack`, which the Apple
linker does not accept.  Before trusting the build, check it against the test
vectors that ship with it:

    cd ref && cc -O3 -fomit-frame-pointer -DKYBER_K=3 kem.c indcpa.c polyvec.c \
       poly.c ntt.c cbd.c reduce.c verify.c fips202.c symmetric-shake.c \
       test/test_vectors.c -o /tmp/test_vectors768
    /tmp/test_vectors768 | shasum -a 256   # must equal tvecs768 in SHA256SUMS

On the machine of the paper this reports 101 000, 110 000 and 139 000 cycles for
key generation, encapsulation and decapsulation, and confirms
`|ct| = 1088`, `|pk| = 1184`, `|sk| = 2400`.

Every row of the paper's comparison table is the implementation its authors
published, which for this scheme is the reference C above and not an optimised
build.  A vector implementation of the same KEM is faster, so the table is not a
statement about how fast `ML-KEM-768` can be made to run.
