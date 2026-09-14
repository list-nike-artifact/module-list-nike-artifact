# Parameters

One script.  It re-derives the parameter table of the submission rather than
reading it, and prints each check with the arithmetic that produced it.

    sh fetch-estimator.sh
    sage verify-params.sage lattice-estimator

    # or, if the estimator is already somewhere:
    sage verify-params.sage /path/to/lattice-estimator
    ESTIMATOR_PATH=/path/to/lattice-estimator sage verify-params.sage

    # and to see the neighbouring parameter sets the text refers to:
    sage verify-params.sage lattice-estimator --variants

Four checks, of which only the last needs the estimator.  Without it the script
reports the first three and says it is skipping the fourth.

1. **Direct.** Runs honest exchanges of the scheme, in Sage, and verifies the gap
   identity of the paper against them.
2. **The hypothesis (H)** relating the noise to the reconciliation radius, and the
   correctness bound of the reconciliation lemma that follows from it.
3. **The LoHL condition**, which is what forces `m` to be large and therefore what
   sets the public-key size.
4. **MLWE security**, via the lattice estimator, which reports `2^133` against the
   128-bit target.

`fetch-estimator.sh` pins the estimator at
`f18533a19433f6fb1d9fb396006f462adc6b8ad3`, which is the revision the paper cites.
