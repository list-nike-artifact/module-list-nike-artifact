# `X25519`

The `X25519` row of the comparison table, from the reference implementation.

    cc -O3 -Iglue -Iref10 -o bench bench.c ref10/base.c ref10/fe_*.c ref10/scalarmult.c
    ./bench 4.05

The argument is the clock in GHz, used only to turn the measured times into cycle
counts.  `bench.c` checks itself against the test vectors of RFC 7748 Sec. 6.1,
both directions, and refuses to print a figure unless they pass.

`ref10/` is `crypto_scalarmult/curve25519/ref10` from SUPERCOP, by D. J.
Bernstein, in the public domain; `sh fetch.sh` diffs it against upstream and is
the check that it is unmodified.  `glue/` is five headers of ours that SUPERCOP
generates at build time and that are therefore not part of the directory: four
integer typedefs and one declaration.  Each says in a comment what it is for.

Key generation and shared-key derivation cost the same here, and that is not a
mistake: `ref10` computes a public key with the same Montgomery ladder it uses
for a shared key, over the fixed base point, without the precomputed tables that
an optimised implementation would use.

On the machine of the paper, idle and on AC power, this reports

    X25519       keygen 217 000   derive 217 000

Nothing here is pinned to a core the way a kernel-level benchmark would be, so an
unquiet machine inflates every line, by 20% or more, and unevenly between them.
Check that the load average is near zero before believing a figure; the frequency
calibration that `../module-list-nike/` prints is the quickest way to tell
whether the core is where you think it is.
