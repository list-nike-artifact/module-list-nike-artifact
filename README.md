# Artifact: a list-NIKE from module lattices

Everything needed to reproduce the implementation section of the submission: the
scheme, the script that verifies its parameter set, and the four schemes it is
compared against.  Each directory has a README with its own build and run lines.

| directory | what it is |
| --- | --- |
| `coral/` | the `CORAL` reference C implementation, and the script that reduces its own harness's output to medians |
| `mlkem/` | harness for the pq-crystals reference C implementation of `ML-KEM-768` |
| `module-list-nike/` | our implementation of the module list-NIKE (Rust, no dependencies) |
| `params/` | one Sage script, which re-derives the parameter table of the paper |
| `pswoosh/` | `Passive-SWOOSH` [USENIX Security 2024], ported to aarch64 so it can be measured on the same machine |
| `x25519/` | harness for `X25519`, over the `ref10` reference implementation |

Only `module-list-nike/`, `params/` and `pswoosh/` contain code of ours.  In
`coral/`, `mlkem/` and `x25519/` the implementation is upstream's, vendored
unmodified so that the artifact builds without network access, and only a benchmark
harness and a README are ours.  Each of those three carries a `fetch.sh` that
clones upstream at the revision we measured and diffs the vendored tree against it,
so that nothing about those schemes has to be taken on trust:

    sh coral/fetch.sh && sh mlkem/fetch.sh && sh x25519/fetch.sh

`x25519/` names no revision, because SUPERCOP is distributed as dated tarballs
rather than as a git repository; its `fetch.sh` diffs against the current upstream
files instead.

Start with `module-list-nike/README.md`, which maps every line of the scheme's two
figures to a place in the code and lists the tests that pin the fast paths to
literal transcriptions of the paper's definitions.

## Reproducing the figures

Every figure in the submission was measured on one performance core of an Apple M3
at 4.05 GHz, as the median over 1000 independent calls, and every harness here
follows the same three conventions so that the figures can be compared with one
another:

1. the process asks the scheduler for a performance core
   (`pthread_set_qos_class_self_np(0x21, 0)`);
2. it busy-waits for about 300 ms first, so the core is already at its top
   frequency when timing starts;
3. it times with the monotonic raw clock and converts at a frequency given on the
   command line, which is why each harness takes the clock in GHz as an argument.

Our two Rust harnesses need nothing but `cargo`:

    cargo run --release --manifest-path module-list-nike/Cargo.toml --bin bench
    cargo run --release --manifest-path pswoosh/Cargo.toml --bin bench

The first is Table 4 of the paper, the second the `Passive-SWOOSH` row of Table 5.
Every command in this file is written to be run from this directory, the root of
the artifact.

The other three build a vendored upstream tree; see `mlkem/README.md`,
`x25519/README.md` and `coral/README.md`.

**Run the machine on AC power with Low Power Mode off, and otherwise idle.**  On
battery the core settles about 12% below its rated frequency and in Low Power Mode
roughly half of it, and a machine that is merely busy is worse still: every figure
then comes out too large, by 20% or more and by different factors for different
schemes.  Both Rust harnesses print an in-process calibration of the core they are
running on next to the 4.05 GHz they convert at, so

    calibrated clock: 3.73 GHz (using 4.05)

means the core was throttled and the run should be thrown away.  On a machine with
a different clock, pass its own frequency instead of 4.05; cycle counts should
carry over, wall-clock times will not.

## Parameters

    cd params
    sh fetch-estimator.sh
    sage verify-params.sage lattice-estimator

`verify-params.sage` is the only script here, and it re-derives the parameter table
of the submission rather than reading it: the sizes from the parameters, the
correctness bound, the hypothesis relating the noise to the reconciliation radius,
and the MLWE security estimate.  Only the last needs the lattice estimator; run it
with no argument and it reports the other three and says it is skipping that one.
