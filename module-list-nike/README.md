# `module-list-nike` — the module list-NIKE of the paper, in Rust

The construction of Figure `fig:mod-list-nike` at the recommended parameter set
of Table `tab:mod-params`:

    d = 128   rho = 2   c = 128   beta = 16   n = 9   m = 100   q = 2^23 - 21
    eta = 2   ell = 2^24      MLWE 133 bits, correctness 2^-103
    |pk| = 40112 B   |sk| = 12912 B   |shk| = 32 B

Under 1900 lines, no dependencies, no vector intrinsics and no assembly, and it
follows the pseudocode of the figure line by line so that it can be read against
the paper.  Nothing needs the network.

    cargo test    --release            # 23 tests, see "How to check it" below
    cargo run     --release --bin check    # one page of checks against the paper
    cargo run     --release --bin bench    # the cycle counts of Table 4
    cargo run     --release --bin listpos  # where the peer's key sits in the list

## Figures

Median over 1000 independent calls on one performance core of an Apple M3 at
4.05 GHz.  These are Table 4 of the paper.

| | cycles | time |
| --- | --- | --- |
| `Stp`: sample `A` in `R_q^{9 x 100}` | 442 500 | 109 us |
| `Gen`: `(sk, pk)`, `2nm = 1800` products in `R_q` | 4 980 000 | 1.23 ms |
| `Sdk`, initiator | 291 100 | 71.9 us |
| &nbsp;&nbsp;of which `ListRec` | 2 361 | 0.6 us |
| `Sdk`, responder, honest peer | 291 300 | 71.9 us |
| &nbsp;&nbsp;of which `ListRec` | 2 527 | 0.6 us |
| `Sdk`, responder, worst case | 145 370 000 | 35.9 ms |
| &nbsp;&nbsp;of which `ListRec` | 145 080 000 | 35.8 ms |

Run to run these move by about 2%, so the digits past the third are not
meaningful; the calibration line `bench` prints first is what says whether a run
is usable at all.

`Sdk` returns a whole list of `ell` candidate keys, and correctness only asks
that the first entry of one party's list occur somewhere in the list of its peer.
The initiator, the party that names the key, therefore stops at the first entry:
`ListRec` is step 5 of `Sdk` and returns the list in a compact form that already
determines that entry.  The responder has to recognise which of its own `ell`
candidates the named key is, so it expands its list and compares entry by entry,
at 8.6 cycles an entry.

Enumerating all `2^24` entries is a bound and not the expected cost, which is
what `bin/listpos` measures.  Of 1000 honest pairs, 915 needed no deviation at
all on the 24 hedged coefficients and the other 85 needed exactly one, and none
needed more than three.  Enumerating lightest first, as the code does, the peer's
key is found after 1 entry at the median, 22 at the 99th percentile and 25 at the
worst, 2.1 on average; in index order the same sample needs up to `2^23` entries
and 53 000 on average.

Inside `Sdk`: `k_stat = s^T x` is `n = 9` ring products at 23 500 cycles and
`k_comp = b^T u` is `m = 100` at 258 900, so `k_comp` is 89% of what the
initiator pays.  Inside `Gen`: the 1800 ring products at about 2750 cycles each
are 4.95 M, so the sampling of `s`, `e` and `u` together is under 1%.

## Reading the code against the paper

Each line of the two figures has one place in the code.

| paper | code |
| --- | --- |
| `R_q = Z_q[X]/(X^d+1)` | `poly.rs` (`mac`, `mul_ref`) over `fq.rs` |
| `NIKEStp`: `A <- R_q^{n x m}` | `nike::stp` |
| `NIKEGen`: `s <- R_q^n` | `sample::uniform` |
| `e <- CbdVec(eta)` | `sample::cbd` |
| `u <- ([-beta,beta]^d)^m` | `sample::short` |
| `pk_L := s^T A + e^T` | `nike::gen`, first loop |
| `pk_R := A u` | `nike::gen`, second loop |
| `NIKESdk`: `k_stat := s^T x` | `nike::k_stat` |
| `k_comp := b^T u` | `nike::k_comp` |
| `ListRec(k_stat + k_comp)`, step 5 | `nike::raw_key` then `listrec::plan` |
| `Rnd`, `Rnd_rec`, `DistanceToRB` | `listrec::round_coeff` (all three from one division) |
| sort `T` by distance, hedge the first `log ell` | `listrec::plan` |
| `L[k] <- L[k] || L_i[...]` | `Plan::entry` (literal) and `Plan::enumerate_into` (fast) |

## How to check it

Every fast path is pinned to a literal transcription of the paper's definition by
a test, and `bin/check` re-runs the interesting ones with output you can read.

* `listrec::rounding_matches_definitions` — `round_coeff` against literal
  `floor(k 2^rho / q)`, `floor(k 2^(rho+1) / q) mod 2`, and a slow minimum over
  all `2^rho + 1` scaled boundaries.
* `listrec::enumeration_matches_direct_construction` — the `O(ell)` toggling
  enumeration against `Plan::entry`, which builds each entry from the figure's
  indexing directly, for 200 000 entries.
* `listrec::hedged_are_the_closest_to_a_boundary` — the hedged set really is the
  `log ell` nearest ones.
* `poly::mac_matches_schoolbook` — the expanded-form product against textbook
  negacyclic convolution; `negacyclic_relation` checks `X^d = -1`.
* `fq::reduce_matches_rem_euclid` — the pseudo-Mersenne reduction against
  `rem_euclid`, on 200 000 values up to `2^57` plus the edge cases.
* `aes::fips197_vector` — the AES-128 block against the FIPS 197 test vector, and
  `neon_matches_reference` against the portable path.
* `nike::raw_key_gap_is_inside_the_reconciliation_radius` — the gap
  `k_1 - k_2 = e_2^T u_1 - e_1^T u_2` of `eq:mod-gap` stays inside
  `q / 2^(rho+1)`.
* `nike::list_correctness_holds_for_honest_pairs` — `ListRec(k_1)[1]` is in
  `ListRec(k_2)` and conversely, over independent key pairs.
* `nike::serialised_sizes_match_the_table` — serialises a freshly generated key
  pair and *measures* the length.  The table's 40112 / 12912 / 32 bytes are the
  expectation; nothing in the code is built from them.

## Design notes

* **No NTT.**  `f = ord_{2d}(q) = 64` gives `g = 2` factors, so at most one
  splitting layer could exist, and it would need `q = 1 (mod 4)`.  Here
  `q = 2^23 - 21 = 3 (mod 4)`, so it does not.  Products are schoolbook, which
  costs little because `q < 2^23`: a whole inner product accumulates in `i64` and
  is reduced once per output coefficient, not once per product.  This reaches
  about 5.9 coefficient products per cycle.
* **Expanded form.**  One operand of every product is stored as
  `[-b_0, ..., -b_{d-1}, b_0, ..., b_{d-1}]`, which turns the negacyclic product
  into a sliding dot product and folds the wrap-around sign into the data.  The
  upper half is the ordinary coefficient vector, so the same array serves both
  roles.
* **The clock.**  `bench` measures wall-clock nanoseconds and multiplies by a
  frequency it also calibrates itself: it times two loops whose bodies are chains
  of 8 and of 24 dependent `add` instructions, so the loop overhead cancels in
  the difference and what remains is 16 cycles per iteration.  On the machine of
  the paper that returns 4.052 GHz against a published performance-core clock of
  4.05 GHz, which is the sanity check for the whole method.  See `cycles.rs`.
* **Constant time in the secret.**  `poly::mac` has no data-dependent branch (an
  early exit on zero coefficients would leak `s`), and `fq` uses masked
  conditional subtraction rather than `if`.  The enumeration of `ListRec` depends
  on the raw key, which both parties share.  The responder's scan, however, stops
  at the matching entry, so its running time reveals where in the list that entry
  lay; a deployment that cares about this should walk the whole list.
