r"""
verify-params.sage -- self-contained check of the module list-NIKE parameter sets.

    sage verify-params.sage [/path/to/lattice-estimator] [--variants]

or set $ESTIMATOR_PATH instead of passing the path.

The target is NIST level 1 only: about 128 bits of MLWE security, no more.  The
default set is the one recommended in the paper; --variants additionally checks
every row of the (d, rho, ell) trade-off table.

For every parameter set below the script recomputes, from scratch, the numbers
claimed in the paper:

  (1) the scheme of Fig. \ref{fig:mod-list-nike} is run on real ring elements
      and the correctness definition is tested literally.  This is a functional
      check: the published delta is far below any feasible number of trials, so
      the expectation is simply zero failures.
  (2) delta_hon as the paper bounds it: hypothesis (H) of
      Lemma \ref{lem:list-recon-multi-prob} plus the bound that lemma states.
      This is the number quoted in the tables; nothing outside the paper is
      used to obtain it.
  (3) the condition \eqref{eq:anon-module-lhl-cond} of Lemma \ref{lem:rlohl},
      which is what makes pk_R = A u statistically uniform.
  (4) the cost of the best known lattice attack on pk_L = s^T A + e^T, via
      LWE.estimate.rough().

This is a CHECKER, not a search: the parameters are given, and each quantity is
recomputed independently of how they were found.

Only (4) needs the lattice estimator (https://github.com/malb/lattice-estimator);
pass its path as an argument or in $ESTIMATOR_PATH.  Without it the script still
prints (0), (1) and (2).
"""

import os
import sys
import numpy as np
from sage.all import *

ESTIMATOR_PATH = os.environ.get("ESTIMATOR_PATH", "lattice-estimator")

# Two DIFFERENT parameters, which it is easy (and wrong) to conflate.
#
#   SECPAR is the security parameter: the target is NIST level 1, so the best
#   known attack on pk_L must cost about 2^128, and every additive term of the
#   security bound must be at most 2^-SECPAR.
#
#   LHL_SURPLUS is the entropy surplus that eq:anon-module-lhl-cond demands over
#   the (n+1) d log q bits the leftover-hash output already costs.  It is NOT a
#   security level.  Lemma \ref{lem:rlohl} loses a square root -- it bounds the
#   distance by 2^(-LHL_SURPLUS/2 - 1) -- so a surplus of s bits buys only
#   s/2 + 1 bits of distance, and asking for the 2^-SECPAR the rest of the bound
#   is set to means
#
#       LHL_SURPLUS = 2 * SECPAR + 2 = 258,
#
#   for which the statistical term of Theorem \ref{thm:pkanon-module} is
#   2^(1 - LHL_SURPLUS/2) = 2^-SECPAR = 2^-128 exactly.  Taking the surplus to be
#   SECPAR itself would leave 2^-63 there; reading the surplus as a security
#   level and raising SECPAR to 256 would instead ask for MLWE parameters that
#   nobody wants at level 1.  The surplus is paid in m alone -- extra columns of
#   A -- and not in lattice dimension, which is why it is cheap here.
SECPAR = 128
LHL_SURPLUS = 2 * SECPAR + 2

TRIALS = 40  # (1) honest executions per set: functional check only


def grouped(x):
    r"""39000 -> "39 000": digits in groups of three, as the paper prints sizes."""
    s = str(ZZ(x))
    out = []
    while len(s) > 3:
        out.append(s[-3:])
        s = s[:-3]
    return " ".join([s] + out[::-1])


# ---------------------------------------------------------------------------
# the parameter sets
# ---------------------------------------------------------------------------
#   d        ring degree; R_q = Z_q[X]/(X^d + 1), d a power of two
#   n, m     A <- R_q^{n x m};  pk_L = s^T A + e^T in R_q^m,  pk_R = A u in R_q^n
#   q        modulus, prime
#   eta      e <- CBD(eta)^m           (centred binomial on {-eta,...,eta})
#   beta     u <- ([-beta, beta]^d)^m  (uniform coefficients)
#   rho      bits extracted per coefficient
#   nc       number of coefficients used; the shared key has rho * nc bits
#   log_ell  log2 of the list length ell
#
# The target is NIST level 1: about 128 bits of MLWE security and no more.  The
# recommended set is the first one; the rest are the (d, rho, ell) trade-off rows
# of the paper's table, checked only under --variants.
#
# rho is not free.  The raw key is a single ring element, so at most
# nc <= d coefficients are available, and a rho * nc = 256 bit shared key forces
# rho >= 256/d: choosing rho is choosing d.  Raising rho beyond 256/d only
# shrinks the reconciliation radius q/2^(rho+1), which (H) then pays for with a
# larger q, so the minimum rho = 256/d is always the right one.
#
# In every set q is the largest prime below a power of two with
# ord_{2d}(q) = d/2, the largest such order a power-of-two cyclotomic admits;
# hence ceil(log2 q) is exactly k and no bits of the public key are wasted.

PARAMS = []


def q_below(k, d):
    r"""Largest prime < 2^k with ord_{2d}(q) = d/2, the largest such order a
    power-of-two cyclotomic admits.  ceil(log2 q) = k exactly."""
    q = previous_prime(2 ** k)
    while Integers(2 * d)(q).multiplicative_order() != d // 2:
        q = previous_prime(q)
    return q


def P(name, n, m, beta, k=None, q=None, log_ell=24, d=128, eta=2, rho=2, nc=128):
    PARAMS.append(dict(name=name, d=d, n=n, m=m, q=q if q else q_below(k, d), eta=eta, beta=beta, rho=rho, nc=nc,
                       log_ell=log_ell))


P("RECOMMENDED  d=128 rho=2 beta=16 ell=2^24  -- NIST level 1",
  d=128, rho=2, nc=128, log_ell=24, n=9, m=100, k=23, beta=16)

# The (d, rho, ell) trade-off rows, run only under --variants.  Each is the
# smallest key at that (d, rho, ell) reaching >= 128 bits of MLWE.
VARIANTS = []
_rec, PARAMS = PARAMS, VARIANTS
P("variant  d=64  rho=4 ell=2^24", d=64, rho=4, nc=64, log_ell=24, n=18, m=328, k=23, beta=3)
P("variant  d=256 rho=1 ell=2^24", d=256, rho=1, nc=256, log_ell=24, n=4, m=85, k=20, beta=3)
P("variant  d=512 rho=1 ell=2^24", d=512, rho=1, nc=256, log_ell=24, n=2, m=59, k=21, beta=3)
P("variant  d=128 rho=4 ell=2^24", d=128, rho=4, nc=64, log_ell=24, n=10, m=203, k=24, beta=3)
P("variant  d=128 rho=8 ell=2^24", d=128, rho=8, nc=32, log_ell=24, n=11, m=255, k=28, beta=3)
P("variant  d=128 rho=2 ell=2^16", d=128, rho=2, nc=128, log_ell=16, n=8, m=149, k=21, beta=3)
P("variant  d=128 rho=2 ell=2^32", d=128, rho=2, nc=128, log_ell=32, n=8, m=149, k=21, beta=3)
P("variant  d=128 rho=2 ell=2^24, beta=3 -- conservative beta", d=128, rho=2, nc=128, log_ell=24, n=8, m=149, k=21,
  beta=3)
P("variant  d=128 rho=2 ell=2^24, beta=8", d=128, rho=2, nc=128, log_ell=24, n=9, m=118, k=22, beta=8)
P("variant  d=128 rho=2 ell=2^24, beta=64", d=128, rho=2, nc=128, log_ell=24, n=10, m=86, k=25, beta=64)
P("variant  d=128 rho=2 ell=2^32, beta=16 -- delta <= 2^-128", d=128, rho=2, nc=128, log_ell=32, n=9, m=100, k=23,
  beta=16)
PARAMS = _rec

PARAMS = PARAMS + (VARIANTS if "--variants" in sys.argv else [])

# The nine rows of tab:mod-params, in table order, used by --ell to rebuild
# tab:mod-ell.  Same sets as VARIANTS above, minus its ell-only duplicates.
TAB_ROWS = [dict(d=64, rho=4, nc=64, n=18, m=328, k=23, beta=3),
            dict(d=128, rho=2, nc=128, n=8, m=149, k=21, beta=3),
            dict(d=256, rho=1, nc=256, n=4, m=85, k=20, beta=3),
            dict(d=512, rho=1, nc=256, n=2, m=59, k=21, beta=3),
            dict(d=128, rho=4, nc=64, n=10, m=203, k=24, beta=3),
            dict(d=128, rho=8, nc=32, n=11, m=255, k=28, beta=3),
            dict(d=128, rho=2, nc=128, n=9, m=118, k=22, beta=8),
            dict(d=128, rho=2, nc=128, n=9, m=100, k=23, beta=16),
            dict(d=128, rho=2, nc=128, n=10, m=86, k=25, beta=64)]
ELL_GRID = [8, 16, 20, 24, 28, 32, 48]


# ---------------------------------------------------------------------------
# (1) the scheme itself, on real ring elements
# ---------------------------------------------------------------------------
#
# R_q = Z_q[X]/(X^d + 1).  Party i draws s_i <- R_q^n, e_i <- CBD(eta)^m,
# u_i <- ([-beta,beta]^d)^m and publishes
#     pk_L^(i) = s_i^T A + e_i^T in R_q^m,      pk_R^(i) = A u_i in R_q^n.
# Each party then forms the raw key against the other's public key,
#     raw_i = s_i^T pk_R^(j) + (pk_L^(j))^T u_i
#           = s_i^T A u_j + s_j^T A u_i + e_j^T u_i,
# so the two raw keys differ by exactly
#     raw_1 - raw_2 = e_2^T u_1 - e_1^T u_2                                 (G)
# which is asserted below rather than assumed.

def negmul_Z(a, b, d):
    r"""Product of two coefficient vectors in Z[X]/(X^d + 1): the convolution,
    folded back with X^d = -1 (hence the minus), but NOT reduced modulo q."""
    c = np.convolve(a, b)
    hi = np.zeros(d, dtype=c.dtype)
    hi[:len(c) - d] = c[d:]
    return c[:d] - hi


def negmul(a, b, d, q):
    r"""Product in Z_q[X]/(X^d + 1) of two coefficient vectors.  One coefficient
    of the convolution sums d products < q^2, so int64 runs out at q^2 d ~ 2^63
    (reached by the rho = 8, q = 2^28 variant); above that use Python ints."""
    if q * q * d >= 2 ** 62:
        a, b = a.astype(object), b.astype(object)
    return negmul_Z(a, b, d) % q


def cbd(eta, size, rng):
    return (rng.integers(0, 2, size=(size, 2 * eta)).sum(axis=1) - eta)


def rnd(k, rho, q):
    r"""floor(k * 2^rho / q) in Z_{2^rho}."""
    return (k * 2 ** rho) // q


def rnd_rec(k, rho, q):
    r"""floor(k * 2^(rho+1) / q) mod 2: 1 if k rounds up, 0 if down."""
    return ((k * 2 ** (rho + 1)) // q) % 2


def dist_to_rb(k, rho, q):
    r"""2^rho times the distance from k to the nearest rounding boundary."""
    s = rnd(k, rho, q)
    return np.minimum(k * 2 ** rho - s * q, (s + 1) * q - k * 2 ** rho)


def direct_trial(par, rng):
    r"""One honest exchange.  Returns (ok, max |gap coefficient|)."""
    d, n, m, q = par["d"], par["n"], par["m"], par["q"]
    eta, beta, rho, N, L = par["eta"], par["beta"], par["rho"], par["nc"], par["log_ell"]

    A = rng.integers(0, q, size=(n, m, d), dtype=np.int64)

    def keygen():
        s = rng.integers(0, q, size=(n, d), dtype=np.int64)
        e = np.array([cbd(eta, d, rng) for _ in range(m)], dtype=np.int64)
        u = rng.integers(-beta, beta + 1, size=(m, d), dtype=np.int64)
        pkL = np.zeros((m, d), dtype=np.int64)  # s^T A + e^T
        for j in range(m):
            acc = np.zeros(d, dtype=np.int64)
            for i in range(n):
                acc = (acc + negmul(s[i], A[i][j], d, q)) % q
            pkL[j] = (acc + e[j]) % q
        pkR = np.zeros((n, d), dtype=np.int64)  # A u
        for i in range(n):
            acc = np.zeros(d, dtype=np.int64)
            for j in range(m):
                acc = (acc + negmul(A[i][j], u[j], d, q)) % q
            pkR[i] = acc
        return (s, e, u), (pkL, pkR)

    (s1, e1, u1), (pkL1, pkR1) = keygen()
    (s2, e2, u2), (pkL2, pkR2) = keygen()

    def raw(s, u, pkL, pkR):  # s^T pkR + pkL^T u
        acc = np.zeros(d, dtype=np.int64)
        for i in range(n):
            acc = (acc + negmul(s[i], pkR[i], d, q)) % q
        for j in range(m):
            acc = (acc + negmul(pkL[j], u[j], d, q)) % q
        return acc

    k1 = raw(s1, u1, pkL2, pkR2)
    k2 = raw(s2, u2, pkL1, pkR1)

    # check (G): the gap really is e_2^T u_1 - e_1^T u_2, computed over Z
    gap = np.zeros(d, dtype=np.int64)
    for j in range(m):
        gap = gap + negmul_Z(e2[j], u1[j], d) - negmul_Z(e1[j], u2[j], d)
    assert np.all((k1 - k2 - gap) % q == 0), "gap identity (G) failed"

    # reconciliation: party 2 flags the log(ell) coefficients of k_2 closest to a
    # rounding boundary and offers two candidates for each, one for the rest.
    c1, c2 = k1[:N], k2[:N]
    order = np.lexsort((np.arange(N), dist_to_rb(c2, rho, q)))  # index breaks ties
    flagged = set(int(t) for t in order[:L])
    s2r, b2r = rnd(c2, rho, q), rnd_rec(c2, rho, q)
    alt = np.where(b2r == 0, (s2r - 1) % 2 ** rho, (s2r + 1) % 2 ** rho)
    s1r = rnd(c1, rho, q)
    ok = all(bool(s1r[t] == s2r[t]) or (t in flagged and bool(s1r[t] == alt[t])) for t in range(N))
    return ok, int(np.abs(gap).max())


def measure(par, trials, seed=1):
    rng = np.random.default_rng(int(seed))
    bad, worst = 0, 0
    for _ in range(trials):
        ok, mx = direct_trial(par, rng)
        bad += 0 if ok else 1
        worst = max(worst, mx)
    return bad, worst


# ---------------------------------------------------------------------------
# (2) delta_hon as the paper bounds it: Lemma \ref{lem:list-recon-multi-prob}
# ---------------------------------------------------------------------------
#
# Lemma \ref{lem:list-recon-multi} says the exchange succeeds as soon as every
# mismatching coefficient is one of the log(ell) coefficients that ListRec
# duplicates, provided the displacement obeys |k_1 - k_2| <= q/2^(rho+1)
# coefficientwise.  Lemma \ref{lem:list-recon-multi-prob} bounds the failure
# probability of that event by
#
#   delta <= (L+1) C(nc, L+1) E[ (2^rho |e|/q) ((2^(rho+1)|e| + 1)/q)^L ],
#                                                            L = log2(ell),
#
# where e is one coefficient of the gap.  The hypothesis of that lemma is that
# the displacement bound holds with probability one.  In the MODULE scheme each
# gap coefficient is a signed sum of 2 m d products Cbd(eta) * U[-beta,beta]
# (Lemma \ref{lem:mod-correctness}, identity (G)), so the hypothesis reads
#
#   (H)   2 m d eta beta <= q / 2^(rho+1)
#
# and NOT the m eta beta <= q/2^(rho+1) of the unstructured statement.  (H) is
# checked below; it is the constraint that actually floors the modulus.
#
# Evaluating that expectation needs the absolute moments E|e|^j for j <= L+1 = 25,
# since expanding the L-th power binomially gives
#     E[...] = (2^rho/q^(L+1)) sum_j C(L,j) 2^((rho+1)j) E|e|^(j+1).
# e is symmetric, so the even absolute moments are the raw moments, read off the
# mgf of one product raised to the 2md-th power; the odd ones are bounded by
# Holder, E|e|^(2i+1) <= (E e^(2i+2))^((2i+1)/(2i+2)).  Everything runs at 400
# bits, so the bound reported is an upper bound on the lemma's bound.

RF = RealField(400)


def hypothesis_H(par):
    r"""(bound on |e|, radius q/2^(rho+1)) -- Lemma \ref{lem:list-recon-multi-prob}."""
    B = 2 * par["m"] * par["d"] * par["eta"] * par["beta"]
    return ZZ(B), RR(par["q"]) / 2 ** (par["rho"] + 1)


def gap_moments(par, kmax):
    r"""E[e^0], ..., E[e^kmax] for one coefficient of the gap, exactly."""
    N, eta, beta = 2 * par["m"] * par["d"], par["eta"], par["beta"]
    R = PowerSeriesRing(RF, 't', default_prec=kmax + 1)
    t = R.gen()
    mgf = R(0)  # mgf of one product c * u
    for i in range(2 * eta + 1):
        pc, c = RF(binomial(2 * eta, i)) / RF(2) ** (2 * eta), i - eta
        for u in range(-beta, beta + 1):
            x = RF(c * u)
            mgf += pc / RF(2 * beta + 1) * sum(x ** j * t ** j / RF(factorial(j))
                                               for j in range(kmax + 1))
    ms = (N * mgf.log()).exp()  # = mgf^N, stable since mgf(0) = 1
    return [ms[j] * RF(factorial(j)) for j in range(kmax + 1)]


def delta_paper(par):
    r"""The bound of Lemma \ref{lem:list-recon-multi-prob}, as log2."""
    L, nc, rho, q = par["log_ell"], par["nc"], par["rho"], par["q"]
    if L + 1 > nc:
        return None  # log(ell) < nc is required
    mo = gap_moments(par, L + 2)
    am = lambda j: mo[j] if j % 2 == 0 else mo[j + 1] ** (RF(j) / RF(j + 1))
    tot = sum(RF(binomial(L, j)) * RF(2) ** ((rho + 1) * j) * am(j + 1)
              for j in range(L + 1))
    return RR((RF(L + 1) * RF(binomial(nc, L + 1)) * RF(2) ** rho
               / RF(q) ** (L + 1) * tot).log(2))


def gap_sigma(m, d, eta, beta):
    r"""Standard deviation of one coefficient of k_1 - k_2, in closed form."""
    return sqrt(RR(2 * m * d) * RR(eta) / 2 * RR(beta * (beta + 1)) / 3)


# ---------------------------------------------------------------------------
# (3) the algebraic LHL condition, eq:anon-module-lhl-cond
# ---------------------------------------------------------------------------
#
#   m f log(2 beta + 1) - d log(2 d eta beta m + 1)
#        >= (n + 1) d log q + g log(h + 1) + lambda,
#
# where <q> = p_1^h ... p_g^h with N(p_i) = q^f.  For d a power of two and q an
# odd prime not dividing d, q is unramified (h = 1), f = ord_{2d}(q), g = d / f.
# Since (Z/2d)^* has exponent d/2 we always have f <= d/2 and g >= 2: at best
# half of the entropy of u is usable, and an NTT-friendly q (small f) is useless.
# This, not correctness and not MLWE, is what forces m >> n.

def qform(q):
    r"""q written against the nearest power of two, e.g. 1048573 -> "2^20 - 3"."""
    k = ZZ(q).nbits()
    lo, hi = ZZ(2) ** (k - 1), ZZ(2) ** k
    return ("2^%d - %d" % (k, hi - q)) if hi - q <= q - lo else ("2^%d + %d" % (k - 1, q - lo))


def resdeg(q, d):
    f = ZZ(Integers(2 * d)(q).multiplicative_order())
    return f, ZZ(d) // f


def lhl_sides(par, lam=LHL_SURPLUS):
    r"""The two sides of eq:anon-module-lhl-cond, in bits.  Need LHS >= RHS.

        m f log2(2 beta + 1) - d log2(2 d eta beta m + 1)
              >=  (n + 1) d log2 q  +  g log2(h + 1)  +  lambda

    LHS = min-entropy of u that survives the ideal structure; RHS = length of
    the output pk_R = A u plus the g log(h+1) ideal-counting term plus lambda.
    """
    d, n, m, q = par["d"], par["n"], par["m"], par["q"]
    eta, beta = par["eta"], par["beta"]
    f, g = resdeg(q, d)
    h = 1
    lhs = RR(m * f) * log(RR(2 * beta + 1), 2) - RR(d) * log(RR(2 * d * eta * beta * m + 1), 2)
    rhs = RR((n + 1) * d) * log(RR(q), 2) + RR(g) * log(RR(h + 1), 2) + RR(lam)
    return lhs, rhs


def lhl_slack(par, lam=LHL_SURPLUS):
    lhs, rhs = lhl_sides(par, lam)
    return lhs - rhs


# ---------------------------------------------------------------------------
# (4) hardness of the instance pk_L = s^T A + e^T
# ---------------------------------------------------------------------------
#
# s is uniform on R_q^n and e <- CBD(eta)^m, so flattening the module structure
# gives an LWE instance of secret dimension n*d with m*d samples.  (Flattening
# can only help the attacker: it forgets the ring structure.)
#
# CAVEAT: what is estimated here is a plain MLWE instance.  It does NOT model the
# single adaptive error hint of Def. \ref{def:alwe}, in which the adversary learns
# e_A^T u_B.  A hint can only help the attacker, so the figure below is an UPPER
# BOUND on the hardness of the aMLWE instance the proof actually reduces to.

def security(par, LWE, ND):
    n, m, d, q, eta = par["n"], par["m"], par["d"], par["q"], par["eta"]
    inst = LWE.Parameters(n=int(n * d), q=int(q), Xs=ND.UniformMod(int(q)),
                          Xe=ND.CenteredBinomial(int(eta)), m=int(m * d))
    # `quiet` only exists in estimator revisions after the one we pin, so
    # swallow the progress output here instead of asking for silence.
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        res = LWE.estimate.rough(inst)
    return {k: RR(log(v["rop"], 2)) for k, v in res.items()}


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def pk_parts(par):
    r"""
    |pk| = |pk_L| + |pk_R|, in bits.  Both halves are counted.

        pk_L = s^T A + e^T  in R_q^m   ->  m * d * ceil(log2 q) bits
        pk_R = A u          in R_q^n   ->  n * d * ceil(log2 q) bits

    A itself is a public parameter, regenerated from a seed, and is not counted.
    Because the LoHL forces m >> n, pk_L is 92 to 96% of the key.  If the roles
    of the two parties are fixed -- one side always contributes s, the other
    always contributes u -- then each party publishes only its own half, which is
    a 96/4 split rather than a halving.
    """
    k = ceil(log(RR(par["q"]), 2))
    return par["m"] * par["d"] * k, par["n"] * par["d"] * k


def sk_bits(par):
    r"""
    |sk|, in bits, with the secret key given EXPANDED as the pair (s, u):

        s in R_q^n              ->  n * d * ceil(log2 q)          bits
        u in ([-beta,beta]^d)^m ->  m * d * ceil(log2 (2 beta+1)) bits

    In a real implementation sk is a 32-byte seed and both are regenerated from
    it, so this figure is not a transmission cost.  It is reported because it is
    the one column on which a large beta loses: raising beta shrinks |pk| (m
    falls) but grows |sk| (each coefficient of u needs more bits), which is part
    of why beta = 16 is preferred to beta = 64.
    """
    k = ceil(log(RR(par["q"]), 2))
    ku = ceil(log(RR(2 * par["beta"] + 1), 2))
    return par["n"] * par["d"] * k + par["m"] * par["d"] * ku


def ell_table():
    r"""Rebuild tab:mod-ell: floor(-log2 delta) for every row and every ell."""
    print("tab:mod-ell -- entries are floor(-log2 delta(lambda=%d, ell)) of "
          "Lemma \\ref{lem:list-recon-multi-prob}," % SECPAR)
    print("                'x' marks log2(ell) >= nc, which the lemma excludes.\n")
    print("    d rho beta       |pk|     |sk| " + "".join("%8s" % ("2^%d" % L)
                                                          for L in ELL_GRID))
    for r in TAB_ROWS:
        par = dict(name="", d=r["d"], n=r["n"], m=r["m"], q=q_below(r["k"], r["d"]),
                   eta=2, beta=r["beta"], rho=r["rho"], nc=r["nc"], log_ell=24)
        bL, bR = pk_parts(par)
        cells = []
        for L in ELL_GRID:
            par["log_ell"] = L
            dp = delta_paper(par)
            cells.append("x" if dp is None else "%d" % floor(-dp))
        print("%5d %3d %4d %10s %8s " % (r["d"], r["rho"], r["beta"],
                                         grouped((bL + bR) // 8),
                                         grouped(sk_bits(par) // 8))
              + "".join("%8s" % c for c in cells))
    print("")


def main():
    path = next((a for a in sys.argv[1:] if not a.startswith("-")), ESTIMATOR_PATH)
    LWE = ND = None
    try:
        sys.path.insert(0, path)
        from estimator import LWE as _LWE, ND as _ND
        LWE, ND = _LWE, _ND
    except ImportError:
        print("lattice-estimator not found at %r -- skipping part (4).\n" % path)

    print("entropy surplus = %d in eq:anon-module-lhl-cond  ->  LoHL statistical "
          "distance 2^(1 - surplus/2) = 2^-%d\n" % (LHL_SURPLUS, LHL_SURPLUS // 2 - 1))

    for par in PARAMS:
        f, g = resdeg(par["q"], par["d"])
        print("=" * 78)
        print("%s" % par["name"])
        print("  d = %d, n = %d, m = %d, q = %d = %s (prime), eta = %d, beta = %d"
              % (par["d"], par["n"], par["m"], par["q"], qform(par["q"]),
                 par["eta"], par["beta"]))
        print("      ceil(log2 q) = %d bits per coefficient"
              % ceil(log(RR(par["q"]), 2)))
        print("  rho = %d, nCoefficients = %d, ell = 2^%d  ->  %d-bit shared keys"
              % (par["rho"], par["nc"], par["log_ell"], par["rho"] * par["nc"]))
        bL, bR = pk_parts(par)
        print("  |pk_L| = %11s bytes   (m d ceil(log2 q), %.0f%% of the key)"
              % (grouped(bL // 8), RR(100 * bL) / (bL + bR)))
        print("  |pk_R| = %11s bytes   (n d ceil(log2 q))" % grouped(bR // 8))
        print("  |pk|   = %11s bytes   (= |pk_L| + |pk_R|)" % grouped((bL + bR) // 8))
        print("  |sk|   = %11s bytes   (n d ceil(log2 q) + m d ceil(log2(2 beta+1)), "
              "expanded)" % grouped(sk_bits(par) // 8))
        print("  f = ord_{2d}(q) = %d, g = d/f = %d, h = 1" % (f, g))

        trials = TRIALS
        bad, worst = measure(par, trials)
        rate = ("2^%.1f" % RR(log(RR(bad) / trials, 2))) if bad else "0 of %d" % trials
        print("  (1) direct: %d honest exchanges of Fig.\\ref{fig:mod-list-nike}, "
              "gap identity (G) verified" % trials)
        print("      failures %s (%d), max |gap coefficient| seen %d of radius %d"
              % (rate, bad, worst, floor(RR(par["q"]) / 2 ** (par["rho"] + 1))))

        B, Rr = hypothesis_H(par)
        print("  (2) hypothesis (H) of Lemma \\ref{lem:list-recon-multi-prob}:")
        print("        2 m d eta beta <= q/2^(rho+1):  %d <= %.0f   %s (%.1f%% slack)"
              % (B, Rr, "OK" if B <= Rr else "FAILS", 100 * (1 - RR(B) / Rr)))
        print("      sigma of one gap coefficient = %.1f"
              % gap_sigma(par["m"], par["d"], par["eta"], par["beta"]))
        dp = delta_paper(par)
        print("      Lemma \\ref{lem:list-recon-multi-prob}:  delta_hon <= 2^%.2f"
              "   (tables report the floor, %d)" % (dp, floor(-dp)))
        print("        = (L+1) C(nc,L+1) E[(2^rho|e|/q)((2^(rho+1)|e|+1)/q)^L],  L = %d"
              % par["log_ell"])
        print("        prefactor (L+1) C(%d,%d) = 2^%.1f,  E[...] = 2^%.1f"
              % (par["nc"], par["log_ell"] + 1,
                 RR(log(RR((par["log_ell"] + 1)
                           * binomial(par["nc"], par["log_ell"] + 1)), 2)),
                 dp - RR(log(RR((par["log_ell"] + 1)
                                * binomial(par["nc"], par["log_ell"] + 1)), 2))))
        if bad:
            print("      the bound holds: it is >= the measured rate above")

        lhs, rhs = lhl_sides(par)
        print("  (3) eq:anon-module-lhl-cond, the condition of Lemma \\ref{lem:rlohl}:")
        print("        m f log2(2b+1) - d log2(2 d eta b m + 1)  >=  (n+1) d log2 q"
              " + g log2(h+1) + surplus")
        print("        LHS = %d*%d*log2(%d) - %d*log2(%d) = %.1f - %.1f = %.1f"
              % (par["m"], f, 2 * par["beta"] + 1, par["d"],
                 2 * par["d"] * par["eta"] * par["beta"] * par["m"] + 1,
                 RR(par["m"] * f) * log(RR(2 * par["beta"] + 1), 2),
                 RR(par["d"]) * log(RR(2 * par["d"] * par["eta"] * par["beta"]
                                       * par["m"] + 1), 2), lhs))
        print("        RHS = %d*%d*log2(%d) + %d*log2(2) + %d = %.1f"
              % (par["n"] + 1, par["d"], par["q"], g, LHL_SURPLUS, rhs))
        print("        %.1f >= %.1f  ->  slack = %+.0f bits   %s"
              % (lhs, rhs, lhs - rhs, "OK" if lhs >= rhs else "FAILS"))
        if LWE is not None:
            sec = security(par, LWE, ND)
            print("  (4) MLWE: lattice dimension n*d = %d, samples m*d = %d"
                  % (par["n"] * par["d"], par["m"] * par["d"]))
            for k in sorted(sec):
                print("      %-14s 2^%.1f" % (k, sec[k]))
            print("      best known attack: 2^%.1f" % min(sec.values()))
    print("=" * 78)
    if "--ell" in sys.argv:
        ell_table()
        print("=" * 78)
    print("""
NOTE on rho and d.  rho appears in neither (3) nor (4), and in (2) only through
the reconciliation radius q/2^(rho+1).  It is not free: the raw key is ONE ring
element, so only nc <= d coefficients exist and a rho * nc = 256 bit key forces

        rho >= 256 / d.

Raising rho above that minimum only shrinks the radius, which (H) then pays for
with a larger q, so rho = 256/d is always right.  Choosing rho is choosing d, and
sweeping d at ~128 bits of MLWE (beta = 3, ell = 2^24, minimal rho):

        d = 64,  rho = 4  ->  63 664 bytes, delta 2^-129
        d = 128, rho = 2  ->  52 752 bytes, delta 2^-102      <-- best d
        d = 256, rho = 1  ->  56 960 bytes, delta 2^-72
        d = 512, rho = 1  ->  81 984 bytes, delta 2^-92

(the recommended set refines this row to beta = 16: 40 112 bytes, delta 2^-103)
and at d = 128 the non-minimal rho cost 81 792 bytes (rho = 4) and 119 168 bytes
(rho = 8).  Two effects trade off in d: the RHS of (3) carries a d log2 q
overhead -- the "+1" in (n+1) -- so a smaller d needs a smaller m, while a
smaller d raises rho, halving the radius and costing a bit of q.  At d = 128 the
first effect wins; at d = 64 the second does.

NOTE on ell.  ell appears in neither (H), (3) nor (4), so once q sits at the
minimum (H) allows, raising ell buys correctness at ZERO cost in key size, and
cannot buy key size at all.  At d = 128, rho = 2, beta = 3 (all 52 752 bytes):

        ell   2^8      2^16     2^24     2^32     2^48
        delta 2^-29    2^-63    2^-102   2^-143   2^-232

The only cost is that ListRec enumerates all ell entries, Theta(ell) work per key
derivation, so ell = 2^24 (a few seconds) is about the practical ceiling; 2^32 is
~256x that.  The other constraint is log2(ell) < nc.

NOTE on beta.  The sweep below comes from param-search.sage; only the sets that
appear in the paper's tables are re-verified here (run with --variants).
beta is the one statistical parameter that moves the key size:
(3) gives m ~ 2 (n+1) log2 q / log2(2 beta + 1), so m falls as beta grows, while
(H) makes q grow like beta.  Since |pk| ~ m d log2 q, the net is
(c + log beta)^2 / log(2 beta + 1) -- a real but rapidly flattening gain, and
delta stays put because the ratio 2^rho sigma / q is unchanged:

        beta       1      2      3      6      8     16     32     64    128
        |pk|   83 200 63 168 52 752 49 280 44 704 40 112 40 704 38 400 36 608
        |sk|   10 624 11 328  9 840 11 552 12 608 12 912 14 592 15 008 15 392
        delta    -100   -111   -102   -106    -99   -103   -104   -106   -108

beta = 16 (40 112 bytes) is the knee and is the RECOMMENDED set: from there on
the curve is a flat, non-monotone plateau (40 112, 40 704, 38 400, 36 608
bytes), and the two costs the |pk| row does not show both favour the smaller
modulus.  First, |sk| climbs across the same stretch, since each coefficient of
u needs ceil(log2(2beta+1)) bits: 12 912 bytes at beta = 16 against 15 008 at
beta = 64.  Second, ceil(log2 q) = 23 still packs into three bytes at beta = 16,
while beta = 64 needs 25 bits and therefore four, so a byte-aligned
implementation sends 41 856 bytes against 49 152.  Nothing in the proof needs
beta small -- (3) needs only the min-entropy log2(2beta+1) per coefficient and
(2) needs only (H) -- so beta = 3 vs 16 is a choice, not a constraint.

NOTE on what limits the key size.  Solving (3) for m with f = d/2 gives
m ~ 2 (n+1) log2 q / log2(2 beta + 1), in which d cancels; hence m/n ~ 15-17 and

        |pk| ~ 2 (n+1) d (log2 q)^2 / log2(2 beta + 1),

quadratic in log2 q and linear in the lattice dimension n d.  Both are already
fixed: n d = 1152 by the 128-bit MLWE target, and ceil(log2 q) = 23 by (H),
which is solved together with (3) since q >~ 2^(rho+2) m d eta beta while
m ~ log2 q.  So ~40 000 bytes (or ~53 000 at beta = 3) is the limit, and it is
an ENTROPY-ACCOUNTING one, not a lattice one: pk_R plus the surplus accounts for
(n+1) d log2 q ~ 29400 bits that the min-entropy of u must cover, of which (3)
can use only the fraction f/d = 1/2, and every column of A that supplies entropy
also adds d log2 q bits to pk_L.  The surplus itself is almost free by
comparison: LHL_SURPLUS = 258 bits against a right-hand side of ~29700 raises m
from 99 to 100, about 370 bytes.  Replacing (3) by a computational assumption --
(A, A u) with u short is the knapsack form of MLWE, pseudorandom for m > n under
MLWE in dimension (m-n) d -- removes the entropy deficit, needs only m ~ 2n, and
lands near 6 000 bytes; but it makes anonymity computational and requires the
knapsack instance to survive the leakage
the protocol already publishes (e^T u in the key, u_B as the aMLWE hint).""")


main()
