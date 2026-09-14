# pswoosh-native

Passive-SWOOSH [USENIX Security 2024], vendored from the authors' reference
Rust implementation so that it can be benchmarked on aarch64 (Apple silicon)
next to the list-NIKE of this repository.

Upstream is the `rust/ref0` tree of the reference implementation published with
that paper, at revision `e3cc35360a30ae2ea3c449c58e4958994d506907`; its
repository is linked from the paper.  We do not vendor it here, so that the
diff below is the whole of what we changed.

Upstream is x86-64 only in exactly two places, and only those two were changed:

1. `src/arithmetic/fq.rs` — upstream links `fq.s`, 2718 lines of Jasmin-emitted
   x86-64 assembly, but the rest of the crate calls only five of its symbols
   (`fp_add`, `fp_sub`, `fp_mul`, `fp_toM`, `fp_fromM`).  They are reimplemented
   here for q = 2^214 - 255 and R = 2^256: 4-limb add/sub with a conditional
   subtraction, and CIOS Montgomery multiplication.  Since upstream's arithmetic
   is hand-written assembly, replacing it with portable Rust would have made
   Passive-SWOOSH look slower than it is, and so flattered the comparison, so
   `mul` is hand-written `aarch64` assembly too (`mul_aarch64`): the carries run
   through the condition flags rather than being materialised into registers, and
   the reduction uses the shape of the modulus, m*Q = (m << 214) - 255*m, in
   place of four limb multiplications per round.  `mul_portable` is the
   straightforward version with `u128` intermediates, kept both as the fallback
   on other architectures and as the oracle that `test_mul_agrees_with_portable`
   checks the assembly against.
2. `src/aes256ctr.rs` — replaces `aesenc-int.c`, which expands the public matrix
   with AES-NI intrinsics, by a table-based AES-256 in counter mode.  Two known
   answer tests (FIPS-197 C.3 and SP 800-38A F.5.5) are included.

Three further mechanical changes, each marked `PORT:` in `src/lib.rs`:
`asm!("popcnt")` becomes `u64::count_ones`; the 37 MB `sysA.rs` table of the
public matrix is dropped and the matrix is built instead by upstream's own
`genmatrix` on a fixed seed, outside every timed region; and `with_stack` runs
the binaries on a 256 MB thread, since one matrix is 8 MB and upstream held it
in a `static`.

    cargo test --release      # 19 tests, including a two-party key agreement
    cargo run --release --bin bench
    cargo run --release --bin kat > kat.txt

    # the same tree with portable Rust arithmetic instead of the assembly,
    # which is how the sensitivity of the figures to the arithmetic is measured
    cargo test --release --features portable-fq
    cargo run --release --features portable-fq --bin bench

On an Apple M3 performance core the two builds differ by 11% on key generation
(44.5 M cycles against 49.4 M) and 9% on a derivation (3.36 M against 3.66 M).
Both are about three times below the 146.9 M and 10.6 M the authors publish for
a Skylake core, so neither choice of arithmetic disadvantages Passive-SWOOSH.

`bin/kat` prints a deterministic 324-line test vector: 289 add/sub/mul triples
over range-spanning probes, Montgomery round-trips, two entries of the expanded
matrix, a serialised ring element, the head and tail of the AES keystream, the
noise polynomial, and a seeded two-party key generation and derivation.  Building
the same `kat.rs` against upstream's `fq.s` and `aesenc-int.c` for
`x86_64-apple-darwin` and running it under Rosetta 2 gives byte-identical output,
which is how this port is checked.  Both builds here, with the assembly and with
`--features portable-fq`, match it.

`bin/bench` pins itself to a performance core, warms up, and reports the median
of 1000 calls.  Apple silicon exposes no cycle counter to unprivileged code, so
it times with the monotonic raw clock and converts at 4.05 GHz, printing an
in-process calibration of the same core as a cross-check.  **Turn Low Power Mode
off and run on AC**: under Low Power Mode the core is capped near 2.2 GHz and
every figure is about twice too large, which the printed calibration will show.
