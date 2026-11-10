#!/usr/bin/env python3
"""
Rung 3, step 1 (method P): DeWitt's 1965 covariant recursion for the heat-kernel coefficient a_2 of a second-order operator on a
BUNDLE -- connection N_mu and endomorphism E as matrices -- in Riemann normal coordinates on random algebraic curvature tensors, the
machinery of rung1_dewitt.py (functions copied verbatim with credit; refactoring the closed rung's script is a housekeeping item).
Method M is hv1974_master.py ('t Hooft-Veltman (3.40)); rule 12 is satisfied for the vector when the two agree.
Summary and plan: ../notes/2026-11-04_rung3_thooft_veltman_summary.md.

Operator D = -(Box_N + E) on sections of a bundle of dimension d:  nabla_mu f = d_mu f + N_mu f,
    Box_N f = g^{mu nu} nabla_mu nabla_nu f = (1/sqrt g) d_mu (sqrt g g^{mu nu} V_nu) + g^{mu nu} N_mu V_nu,   V_nu = nabla_nu f.
Ansatz K = (4 pi t)^{-n/2} Delta^{1/2} exp(-sigma/2t) sum_k t^k a_k, a_k(x, 0) a d x d matrix (bundle index at x on the left, at the
base point on the right); transport equations
    sigma^{;mu} nabla_mu a_k + k a_k = Delta^{-1/2} (Box_N + E)(Delta^{1/2} a_{k-1}),      sigma^{;mu} nabla_mu a_0 = 0, a_0(0) = 1.
In normal coordinates sigma^{;mu} = x^mu, so x.nabla = x.d + (x.N): the connection is NOT in radial gauge (for Levi-Civita on
vectors x^mu Gamma^a_{mu b} = -(1/3) R^a_{mu b nu} x^mu x^nu + ... != 0), a_0 is the parallel propagator along the geodesic, and
each a_k is solved degree by degree: (d + k) P_d = RHS_d - sum_{m>=2} (x.N)_m P_{d-m}.  Degrees: a_0 to x^4, a_1 to x^2, a_2 at 0.
The endomorphism is a FIELD E(x) (rung 1 kept it constant): for the vector E = -Ric(x) from the ansatz metric to x^2.
Input metric (rung 1's, with the x^4 term that makes Box R(0) = 0 on every background, so no Box R or Box E structure is fitted):
    g = delta - (1/3) R x x + (2/45) R.R x x x x.

Bundles:  scalar(xi): d = 1, N = 0, E = -xi R(x).   vector: d = 4, N_mu = Gamma^a_{mu b}(x), E^a_b = -R^a_b(x) (Feynman gauge).
          maxwell = vector - 2 scalar(0).   Standard curvature sign (hv1974_master.py header).

PREDICTIONS (rule 25; written 4 Nov 19:40 before the first run):
  P1 [rung 1]: Tr a_2(0) for scalar(xi) = (xi^2/2 - xi/6 + 1/72) R^2 - Ric^2/180 + Riem^2/180  (rung 1's [a_2] with E = -xi R).
  P2 [memory: Gilkey's a_2 = (1/360) Tr(180 E^2 + 60 R E + 5 R^2 - 2 Ric^2 + 2 Riem^2 + 30 Omega^2); marked]:
      vector: Tr a_2(0) = -R^2/9 + (43/90) Ric^2 - (11/180) Riem^2;  maxwell: (1/180)(-25 R^2 + 88 Ric^2 - 13 Riem^2).
  P3 [the rule-12 check, exact if both methods are right]: (1/2) Tr a_2 - (3.40) of hv1974_master.py = +(d/360)(R^2 - 4 Ric^2 + Riem^2)
      as fitted coefficients, i.e. (1/360)(1, -4, 1) for scalar(any xi) and (1/90)(1, -4, 1) for vector; after Gauss-Bonnet the two
      methods agree exactly: scalar(xi) -> (xi^2/4 - xi/12 + 1/240) R^2 + Ric^2/120; vector -> -R^2/40 + (7/60) Ric^2; maxwell ->
      -R^2/30 + Ric^2/10.
  P4 [exact]: the parallel propagator a_0 is not 1 for the vector (its x^2 part is -(1/6) R^a_{mu b nu} x^mu x^nu, the
      Levi-Civita holonomy to second order [derived by hand from x.N a_0 = -x.d a_0 with Gamma to first order; the run decides]).
  Input checks as rung 1's (normal coordinates, Riemann at 0, Delta^{1/2}) on every background; FIT via exact solve on three backgrounds,
  the rest as residuals (which, as in hv1974_master.py, test bookkeeping: the content is P1-P3).
Status PASS if the input checks, the zero residuals, P1 and P3 hold; FAIL otherwise.  (Review, 4 Nov: for the vector P3 is
numerically equivalent to P2 given the hard-coded (3.40) values, so a P2 miss fails through P3 -- "memory" marks its origin, not
its weight.)  WHAT THESE BACKGROUNDS CANNOT SEE (review, 4 Nov): the x^4 term makes Box R(0) = 0 on every background and here
Tr Box E = -xi Box R (scalar), -Box R (vector), so the two total-derivative coefficients of a_2 (Gilkey's 60 Box E, 12 Box R) and
the relation between them are invisible to this run; and for the vector the fit tests only the sums (-40, 172, -22)/360 of the
Omega^2, E^2, R E and pure-geometry pieces -- the scalar's xi-dependence separates E^2 from R E, the Omega^2 coefficient is
confirmed only in combination (a bundle with Omega independent of Riemann, a U(1) field strength say, would separate it).

SMOKE RECORD (rule 11; sandbox 4 Nov 20:10, one background, --samples 1 --bundles scalar,vector, 935 s on one core):
  input checks PASS; Tr a_1 = R/6 - xi R (scalar), -R/3 (vector); Tr a_2 EXACTLY the predicted values P1 and P2 on the background
  (scalar: 7225 xi^2/288 - 7225 xi/864 + 14967413/466560; vector: 38414491/233280; diff 0 for both) -- the Gilkey memory holds.
  P4 FAILED as written: the x^2 part of a_0 has magnitude 1/6 but the hand derivation's sign was wrong; from x.N a_0 = -x.d a_0 with
  Gamma^a_{mu b} = -(1/3)(R^a_{mu b nu} + R^a_{b mu nu}) x^nu one gets 2 P_2 = -(x.Gamma) P_0 = +(1/3) R^a_{mu b nu} x^mu x^nu, i.e.
  c = +1/6 -- the sign error was in the hand transcription, not the recursion.  Re-scored as c = +1/6, reported per background.
  Timing: 935 s per background single-core (sandbox); the author's full run: --samples 4 --jobs 4, wall ~ one background plus
  overhead, ~16 min at the sandbox's speed, under 10 min if the author's cores keep the 2 Nov ratio.  The tail of the script (the fit,
  P3, GB, maxwell, P4's re-scored path, status) has NOT run on real data: exercised only on mock backgrounds carrying the predicted
  values (review, 4 Nov).  The 'a0_x2_vector': False in the smoke line is the pre-rewrite -1/6 assertion.

Run:  python3 rung3_dewitt_bundle.py --samples 4 --jobs 4 --log results/rung3_dewitt_bundle.log      [the author, 4 Nov: 2118 s with --jobs 1]
Checkpoint/resume (9 Nov): each finished background is written to <log>.ckpt.json; a stopped run continues with --resume (same
command plus the flag; the log is appended to); without --resume a non-empty log is refused, never overwritten; the checkpoint is
ignored if seed or bundles differ.  Delete the log and the .ckpt.json to start over.
Step 2 (8 Nov): the graviton bundle (10-dimensional symmetric tensors, the induced Levi-Civita connection, the endomorphism derived in
graviton_loop/graviton_operator_derive.py as a field E(x)) and 't Hooft-Veltman's ghost (vector, E = +Ric):
      python3 rung3_dewitt_bundle.py --samples 4 --jobs 4 --bundles graviton,ghost --log results/rung3_dewitt_graviton.log
  Predictions (rule 25, 8 Nov 16:40, before the first run): the pre-GB Tr a_2 equals 2 x method M's pre-GB coefficients plus
  (d/180)(1, -4, 1) -- the identity that held for the scalar and the vector -- i.e. graviton (11/36, 56/45, -17/18), ghost
  (2/9, 43/90, -11/180); after Gauss-Bonnet half of it equals method M: graviton (7/24, 7/12), ghost (17/60, 7/30) [the ghost here
  is the plain vector-with-E=+Ric a_2; the FP factor -2 is applied in the sum]; and graviton - 2 ghost after GB = (1/120, 7/20), (5.24).
  RECORD (rule 11).  The author's run of 8-9 Nov (results/rung3_dewitt_graviton.log, 4680 s for four backgrounds with --jobs 4):
  the ghost matched exactly; the graviton's Tr a_2 = (59/36, -55/18, 19/18) pre-GB, (7/12, 7/6) after, so (1/2) Tr a_2 after GB =
  (7/24, 7/12) = method M's h-loop EXACTLY, and graviton - 2 ghost after GB = (1/120, 7/20) = (5.24) -- the physics passed.  The run
  printed status=FAIL because the prediction table and the per-sample "predicted" values for the graviton were derived from a WRONG
  hand copy of method M's pre-GB coefficients ((1/8, 11/15, -1/2) typed; (19/24, -17/12, 1/2) in the M log), the third hand-copied
  number of the week to be wrong.  Fixed 9 Nov: method M's (3.40) is computed pointwise inside this script on the same backgrounds
  and every prediction is derived from it; nothing is typed from a log any more.  With the right M, (1/2) Tr a_2 - M =
  (10/360)(1, -4, 1) exactly, the identity predicted.  Step 2 is closed under rule 12 by that log; its FAIL line stays.
  Timing: ~78 min for four backgrounds in parallel on the author's machine (single-core sandbox: > 25 min, unfinished).
Backend (9 Nov evening): --backend flint does the same exact-rational polynomial arithmetic with python-flint (fmpq_mpoly) instead of
sympy; sympy stays the default and the reference.  Needs  pip install -U python-flint  (0.9.0 here).  The recursion, the bundles and
the fits are the same lines for both backends.
Step 3a, method P (9 Nov): gravity + scalar, the 11-dimensional bundle (h_{ab}, phi) and the ghost with its phi-tilde term, on
backgrounds that carry a scalar field phi-tilde(x) = a polynomial of degree 4 with a random rational coefficient on every monomial of degree 1-4 (seed + 1) besides the curvature:
      python3 rung3_dewitt_bundle.py --samples 4 --bundles gravscalar,ghostscalar --backend flint --log results/rung3_dewitt_gravscalar.log
  E for (h, phi), covariant form of 't Hooft-Veltman (4.31)-(4.33) after the P^-1 replacement (v_a = d_a phi-tilde, W_ab = D_a D_b
  phi-tilde as FIELDS; h lower; all contractions with g^-1(x)):
      (E h)_ab = [the graviton's] - (v_a (hv)_b + (hv)_a v_b) + (1/2)(g_ab vhv + v_a v_b tr h) - (1/4) g_ab v^2 tr h + (1/2) v^2 h_ab - 2 W_ab phi,
      (E phi)  = (1/2) tr W tr h - W^ab h_ab - v^2 phi;          ghost: E^a_b = Ric^a_b - v^a v_b.
  PREDICTIONS (rule 25; written 9 Nov 18:50, before the first run):
  Q1 [exact]: E(0) and Tr Omega^2 of both bundles equal, as matrices, what graviton_loop/hv1974_gravscalar_M.py builds from the printed
      (4.31)-(4.33) at the same (Riemann, v, W(0)) -- its build_E is imported and evaluated in place, nothing typed.  (The covariant
      form above is my hand reduction of that script's X; Q1 is the test of the hand reduction.)
  Q2 [derived in code 9 Nov on the scalar bundle with a symbolic E(x), flat and curved: coefficient +1/6; literature: Gilkey's
      Thm 4.18.16 as quoted in arXiv 1411.7851 has +60 tr(nabla^2 E)/360, arXiv 1507.06308 (85) prints the opposite sign for both
      box terms]:  Tr a_2(0) = 2 x (3.40) + (d/180) GB + (1/6) Box(Tr E)(0)  on every background, both bundles.  Box(Tr E)(0) is
      computed from the same E(x) (normal coordinates: the flat Laplacian of Tr E at 0).  This is the first run in which E is not
      built from curvature alone, so the first that sees the Box E term (zero in steps 1-2 by the ansatz's Box R(0) = 0).
  Q3 [follows if Q2 holds]: the cubic and quartic jets of phi-tilde enter Tr a_2(0) only through Box(Tr E)(0) -- (3.40) knows only v
      and W(0).
  Q4: Tr a_1(0) = (d/6) R + Tr E(0).
  [As first written, 18:50: "Closing 3a under rule 12 = Q1 and Q2 on every background".  WITHDRAWN after the fresh review, below.]
  SMOKE RECORD (rule 11; sandbox 9 Nov 19:00, four backgrounds, 26 s, flint): Q1, Q2, Q4 hold on every background, both bundles;
  without the Box E term the Q2 diff is exactly -(1/6) Box(Tr E)(0).  The ghost bundle agrees between the sympy and flint backends
  (sympy 14 min for one background); gravscalar under sympy was not run (too slow).
  FRESH REVIEW (9 Nov, separate agent; mutation tests on scratch copies) -- what this run does and does NOT establish:
   * It does establish: on an 11-dimensional bundle with mixing blocks and an E not built from curvature, DeWitt's recursion gives
     2 x (3.40) + (d/180) GB + (1/6) Box Tr E.  Two methods for the MAP  E -> divergence.
   * It does NOT test E: Q2 holds for ANY matrix field E(x) (both sides are built from the same Ef; a mutated coefficient in endo_phi
     still gives diff 0), and Q1 compares E(0) with method M's own input, the print.  So E -- the phi-tilde terms of X (4.31), Y (4.32),
     Z (4.33), the ghost's v v term -- is single-source until step 3b derives it from the action.  Rule 12 is NOT satisfied for 3a's
     numbers by this run; "closed" was my overstatement.
   * Nothing here tests the x-dependence of E(x): Tr a_2(0) feels it only through Box(Tr E)(0), taken from the same field on both
     sides (mutations of g(x) -> delta, v^a -> v_a, extra x-dependent factors: status stayed PASS), and the x-dependence of the
     off-diagonal blocks drops out of the trace altogether.  The physical (non-total-derivative) part needs only E(0).  Guard added
     after the review: Tr E(x) = -6 R(x) - v^2(x) for (h, phi) and R(x) - v^2(x) for the ghost, to x^2 (catches the traced mutations;
     the off-diagonal blocks' x-dependence stays untestable by any trace).  The reviewer checked Box(Tr E)(0) = -Box(v^2)(0) from the
     raw jets independently on all four backgrounds, and the +1/6 from the Mehler kernel (my "derived in code" used this same recursion).
   * Q3 is vacuous for the quartic jets: they enter only W's x^2 part, which is off-diagonal -- dropping them changes no printed number.
     The cubic jets enter through Box(v^2) only.
(8 Nov: the body is under a main guard -- on Windows multiprocessing re-imports the file and the unguarded 4 Nov version forked
workers forever, each truncating the log; the log is now refused if non-empty, and a line is printed per finished background.)
"""
import time, argparse, itertools, random, sys
import sympy as sp

out = None


def version_line():
    """interpreter and library versions, for the log (rung 0's scripts had it; dropped since; restored 9 Nov as a default)"""
    import platform, sympy
    try:
        import mpmath; mpv = mpmath.__version__
    except Exception:
        mpv = '-'
    try:
        import flint; flv = flint.__version__
    except Exception:
        flv = '-'
    return f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sympy.__version__}, mpmath {mpv}, python-flint {flv}; {sys.executable}"


def say(*a, flush_ckpt=None):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


N = 4
X = sp.symbols('x1:5', real=True)
xi = sp.Symbol('xi')
DEG = 4
idx = range(N)
QUARTIC = sp.Rational(2, 45)


# ----------------------------------------------------------------- polynomial backends (9 Nov evening)
# 'sympy' is the reference: the functions of rung1_dewitt.py, verbatim.  'flint' does the same exact-rational arithmetic with
# python-flint's fmpq_mpoly (FLINT, C) in the variables x1..x4, xi: a wrapper class FP with the usual operators, so the recursion,
# the bundles and the fits below are the same lines for both.  Only multiply / truncate / differentiate / evaluate-at-0 differ.
SX = X                      # the sympy symbols, always
BACKEND = 'sympy'
_ctx = None


def use_backend(name):
    """set the backend in this process (called at the top of run_background, so pool workers on Windows get it too)"""
    global BACKEND, X, _ctx
    BACKEND = name
    if name == 'flint':
        import flint
        _ctx = flint.fmpq_mpoly_ctx.get(('x1', 'x2', 'x3', 'x4', 'xi'), 'degrevlex')
        X = tuple(FP(gen) for gen in _ctx.gens()[:4])
    elif name == 'sympy':
        X = SX
    else:
        raise ValueError(name)


class FP:
    """an exact polynomial in x1..x4, xi over Q (flint fmpq_mpoly); mixes with ints and sympy numbers / polynomials in xi"""
    __slots__ = ('p',)
    _op_priority = 20.0

    def __init__(self, p):
        self.p = p

    @staticmethod
    def co(o):
        import flint
        if isinstance(o, FP):
            return o.p
        if isinstance(o, int):
            return _ctx.constant(flint.fmpq(o))
        o = sp.sympify(o)
        if o.is_Rational:
            return _ctx.constant(flint.fmpq(int(o.p), int(o.q)))
        P = sp.Poly(o, *SX, xi)
        return _ctx.from_dict({mon: flint.fmpq(int(c.p), int(c.q)) for mon, c in P.terms()})

    def __add__(self, o): return FP(self.p + FP.co(o))
    __radd__ = __add__
    def __sub__(self, o): return FP(self.p - FP.co(o))
    def __rsub__(self, o): return FP(FP.co(o) - self.p)
    def __mul__(self, o): return FP(self.p * FP.co(o))
    __rmul__ = __mul__
    def __neg__(self): return FP(-self.p)
    def __pow__(self, n): return FP(self.p ** int(n))

    def __truediv__(self, o):
        o = sp.sympify(o)
        assert o.is_Rational and o != 0, "FP: division by a rational number only"
        return self * (1 / o)

    def __eq__(self, o): return (self.p - FP.co(o)) == 0
    def __ne__(self, o): return not self.__eq__(o)
    __hash__ = None

    def to_sympy(self):
        return sp.Add(*[sp.Rational(int(c.p), int(c.q)) * sp.Mul(*[v**e for v, e in zip(SX + (xi,), mon)]) for mon, c in self.p.to_dict().items()])

    def filt(self, keep):
        return FP(_ctx.from_dict({m: c for m, c in self.p.to_dict().items() if keep(sum(m[:4]))}))


def trunc(expr, deg=DEG):
    if isinstance(expr, FP):
        return expr.filt(lambda t: t <= deg)
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(expr, *SX)
    return sp.Add(*[c * sp.Mul(*[x**e for x, e in zip(SX, mon)]) for mon, c in P.terms() if sum(mon) <= deg])


def homog(expr, d):
    if isinstance(expr, FP):
        return expr.filt(lambda t: t == d)
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(expr, *SX)
    return sp.Add(*[c * sp.Mul(*[x**e for x, e in zip(SX, mon)]) for mon, c in P.terms() if sum(mon) == d])


def at0(expr):
    """the value at x = 0, always returned as a sympy expression (a rational, or a polynomial in xi)"""
    if isinstance(expr, FP):
        return expr.filt(lambda t: t == 0).to_sympy()
    return sp.expand(expr).subs({x: 0 for x in SX})


def pd(expr, m):
    """d/dx^m"""
    if isinstance(expr, FP):
        return FP(expr.p.derivative(m))
    return sp.diff(expr, SX[m])


def tosp(expr):
    return expr.to_sympy() if isinstance(expr, FP) else sp.expand(expr)


def riemann_sample(rng, nmat=8):
    R = {k: sp.Integer(0) for k in itertools.product(idx, repeat=4)}
    for _ in range(nmat):
        A = [[None] * N for _ in idx]
        for a in idx:
            for b in idx:
                if a <= b:
                    A[a][b] = A[b][a] = sp.Rational(rng.randint(-3, 3), rng.randint(1, 3))
        w = sp.Rational(rng.randint(-2, 2) or 1, rng.randint(1, 2))
        for a, b, c, d in itertools.product(idx, repeat=4):
            R[(a, b, c, d)] += w * (A[a][c] * A[b][d] - A[a][d] * A[b][c])
    return R


def invariants(R):
    Ric = {(b, d): sum(R[(a, b, a, d)] for a in idx) for b in idx for d in idx}
    Rs = sum(Ric[(a, a)] for a in idx)
    Ric2 = sum(Ric[(a, b)]**2 for a in idx for b in idx)
    Riem2 = sum(v**2 for v in R.values())
    return Ric, sp.expand(Rs), sp.expand(Ric2), sp.expand(Riem2)


def geometry(R):
    g = [[None] * N for _ in idx]
    for mu in idx:
        for nu in idx:
            q2 = -sp.Rational(1, 3) * sum(R[(mu, a, nu, b)] * X[a] * X[b] for a in idx for b in idx)
            q4 = QUARTIC * sum(R[(a, mu, b, l)] * R[(l, c, nu, d)] * X[a] * X[b] * X[c] * X[d]
                               for a in idx for b in idx for c in idx for d in idx for l in idx)
            g[mu][nu] = trunc(sp.Integer(mu == nu) + q2 + q4)
    A = [[trunc(g[m][n] - sp.Integer(m == n)) for n in idx] for m in idx]
    A2 = [[homog(A[m][n], 2) for n in idx] for m in idx]
    A2sq = [[trunc(sum(A2[m][k] * A2[k][n] for k in idx)) for n in idx] for m in idx]
    ginv = [[trunc(sp.Integer(m == n) - A[m][n] + A2sq[m][n]) for n in idx] for m in idx]
    L = trunc(sum(A[m][m] for m in idx) - sp.Rational(1, 2) * sum(A2sq[m][m] for m in idx))
    sqrtg = trunc(1 + L / 2 + L**2 / 8)
    isqrtg = trunc(1 - L / 2 + L**2 / 8)
    Dhalf = trunc(1 - L / 4 + L**2 / 32)
    iDhalf = trunc(1 + L / 4 + L**2 / 32)
    return g, ginv, sqrtg, isqrtg, Dhalf, iDhalf


def christoffel(g, ginv):
    dg = [[[trunc(pd(g[a][b], c)) for c in idx] for b in idx] for a in idx]
    G = [[[None] * N for _ in idx] for _ in idx]
    for mu in idx:
        for a in idx:
            for b in idx:
                G[mu][a][b] = trunc(sp.Rational(1, 2) * sum(ginv[mu][l] * (dg[l][b][a] + dg[l][a][b] - dg[a][b][l]) for l in idx))
    return G


def riemann_field(G, g):
    Rup = {}
    for r, s, m, n in itertools.product(idx, repeat=4):
        Rup[(r, s, m, n)] = trunc(pd(G[r][n][s], m) - pd(G[r][m][s], n)
                                  + sum(G[r][m][l] * G[l][n][s] - G[r][n][l] * G[l][m][s] for l in idx), 2)
    Rlow = {}
    for a, s, m, n in itertools.product(idx, repeat=4):
        Rlow[(a, s, m, n)] = trunc(sum(g[a][r] * Rup[(r, s, m, n)] for r in idx), 2)
    return Rup, Rlow


# ----------------------------------------------------------------- bundle pieces (new)
def mzero(d):
    return [[sp.Integer(0)] * d for _ in range(d)]


def meye(d):
    return [[sp.Integer(i == j) for j in range(d)] for i in range(d)]


def madd(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in range(len(A))] for i in range(len(A))]


def mscale(A, c):
    return [[c * A[i][j] for j in range(len(A))] for i in range(len(A))]


def mmul(A, B, deg=DEG):
    d = len(A)
    return [[trunc(sum(A[i][k] * B[k][j] for k in range(d)), deg) for j in range(d)] for i in range(d)]


def mmap(f, A):
    return [[f(A[i][j]) for j in range(len(A))] for i in range(len(A))]


def mtr(A):
    return tosp(sum(A[i][i] for i in range(len(A))))


def box_bundle(f, Nc, ginv, sqrtg, isqrtg, deg):
    """Box_N f = (1/sqrt g) d_mu (sqrt g g^{mu nu} V_nu) + g^{mu nu} N_mu V_nu,  V_nu = d_nu f + N_nu f; f, N_mu d x d matrices."""
    V = [madd(mmap(lambda e, nu=nu: pd(e, nu), f), mmul(Nc[nu], f, deg + 1)) for nu in idx]
    d = len(f)
    tot = mzero(d)
    for mu in idx:
        inner = mzero(d)
        for nu in idx:
            inner = madd(inner, mmap(lambda e, mu=mu, nu=nu: trunc(sqrtg * ginv[mu][nu] * e, deg + 1), V[nu]))
        tot = madd(tot, mmap(lambda e, mu=mu: pd(e, mu), inner))
        tot = madd(tot, mmap(lambda e: trunc(e, deg), mmul(Nc[mu], sum_mat([mmap(lambda e, nu=nu: ginv[mu][nu] * e, V[nu]) for nu in idx]), deg)))
    return mmap(lambda e: trunc(isqrtg * e, deg), tot)


def sum_mat(ms):
    tot = mzero(len(ms[0]))
    for m in ms:
        tot = madd(tot, m)
    return tot


def transport(k, RHS, Nc, d, maxdeg, init=None):
    """solve (x.d + k) a + (x.N) a = RHS degree by degree; for k = 0 the degree-0 part is init."""
    xN = sum_mat([mmap(lambda e, mu=mu: trunc(X[mu] * e, maxdeg + 1), Nc[mu]) for mu in idx])   # starts at degree 2
    assert all(homog(e, 1) == 0 for row in xN for e in row), "x.N has a degree-1 part: N(0) != 0, not handled (review, 4 Nov)"
    xN_h = {m: mmap(lambda e, m=m: homog(e, m), xN) for m in range(2, maxdeg + 1)}
    P = {}
    for dd in range(maxdeg + 1):
        rhs = mmap(lambda e, dd=dd: homog(e, dd), RHS) if RHS is not None else mzero(d)
        for m, Qm in xN_h.items():
            if dd - m in P:
                rhs = madd(rhs, mmul(Qm, P[dd - m], maxdeg), -1)
        if dd + k == 0:
            P[dd] = init
        else:
            P[dd] = mscale(rhs, sp.Rational(1, dd + k))
    return sum_mat([P[dd] for dd in P])


def run_background(R, label, bundles, backend='sympy', phijets=None):
    use_backend(backend)
    rep = {'label': label}
    Ric, Rs, Ric2, Riem2 = invariants(R)
    rep['inv'] = (Rs, Ric2, Riem2)
    g, ginv, sqrtg, isqrtg, Dhalf, iDhalf = geometry(R)
    G = christoffel(g, ginv)
    nc = [trunc(sum(G[m][a][b] * X[a] * X[b] for a in idx for b in idx), 3) for m in idx]
    rep['normal_coords'] = all(v == 0 for v in nc)
    Rup, Rf = riemann_field(G, g)
    rep['riemann_at_0'] = all(sp.expand(at0(Rf[k]) - R[k]) == 0 for k in Rf)
    Rsf = trunc(sum(ginv[a][c] * ginv[b][d] * Rf[(a, b, c, d)] for a in idx for b in idx for c in idx for d in idx), 2)
    rep['boxR_at_0'] = sp.expand(at0(sum(pd(pd(Rsf, m), m) for m in idx)))
    rep['Dhalf_x2'] = tosp(homog(Dhalf, 2) - sp.Rational(1, 12) * sum(Ric[(a, b)] * X[a] * X[b] for a in idx for b in idx)) == 0
    # Ricci field with one index up: R^a_b(x) = g^{ac} R_{cb}, R_{cb} = R^r_{c r b}
    Ricf = {(c, b): trunc(sum(Rup[(r, c, r, b)] for r in idx), 2) for c in idx for b in idx}
    Ricup = [[trunc(sum(ginv[a][c] * Ricf[(c, b)] for c in idx), 2) for b in idx] for a in idx]
    # fields needed by the graviton bundles: R(x), Ric^a_b(x), R^a_{m b n}(x) as polynomials to x^2 (all from the ansatz)
    Rupf = {k: trunc(v, 2) for k, v in Rup.items()}                      # R^a_{s m n}(x)
    if phijets is not None:
        # the background scalar phi-tilde(x) and its covariant jets as fields: v_a = d_a phi (lower), W_ab = d_a d_b phi - Gamma^c_ab v_c
        phit = sum(c * X[0]**m[0] * X[1]**m[1] * X[2]**m[2] * X[3]**m[3] for m, c in phijets)
        vlow = [trunc(pd(phit, a), 3) for a in idx]
        Wlow = [[trunc(pd(vlow[a], b) - sum(G[c][a][b] * vlow[c] for c in idx), 2) for b in idx] for a in idx]
        vup = [trunc(sum(ginv[a][c] * vlow[c] for c in idx), 2) for a in idx]
        v2f = trunc(sum(vup[a] * vlow[a] for a in idx), 2)
        trWf = trunc(sum(ginv[a][b] * Wlow[a][b] for a in idx for b in idx), 2)
        v0 = sp.Matrix([at0(vlow[a]) for a in idx]); W0 = sp.Matrix(4, 4, lambda a, b: at0(Wlow[a][b]))
        rep['scalar_bg'] = (sp.expand((v0.T * v0)[0]), sp.expand(W0.trace()))
        # method M's own construction at the same point, from the sibling script (computed in place, never typed)
        import os
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'graviton_loop'))
        import hv1974_gravscalar_M as GM
        RicM = sp.Matrix(4, 4, lambda b, d: Ric[(b, d)])
        E_M, Om2_M = GM.build_E(R, RicM, Rs, v0, W0)
        Eg_M = RicM - v0 * v0.T
    for name in bundles:
        if name == 'scalar':
            d = 1; Nc = [mzero(1) for _ in idx]; Ef = [[-xi * Rsf]]
        elif name == 'vector':
            d = 4; Nc = [[[G[a][mu][b] for b in idx] for a in idx] for mu in idx]          # (N_mu)^a_b = Gamma^a_{mu b}
            Ef = [[-Ricup[a][b] for b in idx] for a in idx]
        elif name == 'ghost':                                                              # the FP ghost of 't Hooft-Veltman (5.19): E = +Ric
            d = 4; Nc = [[[G[a][mu][b] for b in idx] for a in idx] for mu in idx]
            Ef = [[Ricup[a][b] for b in idx] for a in idx]
        elif name == 'ghostscalar':                                                        # (5.19) with its phi-tilde term: E = Ric - v v
            d = 4; Nc = [[[G[a][mu][b] for b in idx] for a in idx] for mu in idx]
            Ef = [[trunc(Ricup[a][b] - vup[a] * vlow[b], 2) for b in idx] for a in idx]
            rep['trE_field_ghostscalar'] = tosp(trunc(sum(Ef[i][i] for i in idx) - Rsf + v2f, 2)) == 0      # guard (review, 9 Nov)
        elif name in ('graviton', 'gravscalar'):
            # the 10-dimensional bundle of symmetric tensors h_{ab} (lower indices) in the basis E_aa, E_ab + E_ba; the induced
            # connection (N_mu h)_{ab} = -Gamma^c_{mu a} h_{cb} - Gamma^c_{mu b} h_{ac}; the endomorphism of graviton_operator_derive.py
            # as a field, with the index positions of the lower-index h:  (E h)_{ab} = -R h_{ab} + (1/2) R g_{ab} tr h + (Ric h + h Ric)_{ab}
            # - g_{ab} tr(Ric h) - Ric_{ab} tr h + 2 R_{a m b n} h^{m n},  traces with g^{-1}(x).
            d = 10
            SP = [(a, b) for a in idx for b in idx if a <= b]
            def sym_of_coords(c):
                h = [[sp.Integer(0)] * 4 for _ in idx]
                for (a, b), v in zip(SP, c):
                    h[a][b] = h[a][b] + v
                    if a != b:
                        h[b][a] = h[b][a] + v
                return h
            def coords_of_sym(h):
                return [h[a][b] for (a, b) in SP]
            def op_mat(f):
                cols = []
                for p in SP:
                    e = [[sp.Integer(0)] * 4 for _ in idx]
                    e[p[0]][p[1]] += 1
                    if p[0] != p[1]:
                        e[p[1]][p[0]] += 1
                    cols.append(coords_of_sym(f(e)))
                return [[trunc(cols[j][i], 2) for j in range(10)] for i in range(10)]
            def conn(mu):
                def f(h):
                    return [[-sum(G[c][mu][a] * h[c][b] + G[c][mu][b] * h[a][c] for c in idx) for b in idx] for a in idx]
                return f
            Nc = [op_mat(conn(mu)) for mu in idx]
            Rlowf = {}
            for a, m, b, n in itertools.product(idx, repeat=4):
                Rlowf[(a, m, b, n)] = trunc(sum(g[a][r] * Rupf[(r, m, b, n)] for r in idx), 2)
            Riclowf = {(a, b): trunc(sum(Rupf[(r, a, r, b)] for r in idx), 2) for a in idx for b in idx}
            def endo(h):
                trh = sum(ginv[m][n] * h[m][n] for m in idx for n in idx)
                hup = [[sum(ginv[m][p] * h[p][q] * ginv[q][n] for p in idx for q in idx) for n in idx] for m in idx]     # h^{mn}
                trRich = sum(Riclowf[(m, n)] * hup[m][n] for m in idx for n in idx)
                Rich = [[sum(Riclowf[(a, c)] * ginv[c][p] * h[p][b] for c in idx for p in idx) for b in idx] for a in idx]       # (Ric h)_{ab} = Ric_{ac} g^{cp} h_{pb}
                out_ = [[None] * 4 for _ in idx]
                for a in idx:
                    for b in idx:
                        out_[a][b] = trunc(-Rsf * h[a][b] + Rsf / 2 * g[a][b] * trh + Rich[a][b] + Rich[b][a] - g[a][b] * trRich
                                           - Riclowf[(a, b)] * trh + 2 * sum(Rlowf[(a, m, b, n)] * hup[m][n] for m in idx for n in idx), 2)
                return out_
            Ef = op_mat(endo)
            if name == 'gravscalar':
                d = 11
                def endo_phi(h):
                    trh = sum(ginv[m][n] * h[m][n] for m in idx for n in idx)
                    hv = [sum(h[a][b] * vup[b] for b in idx) for a in idx]                      # (h v)_a = h_ab v^b
                    vhv = sum(vup[a] * hv[a] for a in idx)
                    return [[trunc(-(vlow[a] * hv[b] + hv[a] * vlow[b]) + sp.Rational(1, 2) * (g[a][b] * vhv + vlow[a] * vlow[b] * trh)
                                   - sp.Rational(1, 4) * g[a][b] * v2f * trh + sp.Rational(1, 2) * v2f * h[a][b], 2) for b in idx] for a in idx]
                Ephi = op_mat(endo_phi)
                Wup = [[sum(ginv[a][p] * Wlow[p][q] * ginv[q][b] for p in idx for q in idx) for b in idx] for a in idx]
                def bottom(h):                                                                   # <Y, h> = (1/2) tr W tr h - W^ab h_ab
                    return sp.Rational(1, 2) * trWf * sum(ginv[m][n] * h[m][n] for m in idx for n in idx) - sum(Wup[a][b] * h[a][b] for a in idx for b in idx)
                brow = []
                for p_ in SP:
                    e = [[sp.Integer(0)] * 4 for _ in idx]
                    e[p_[0]][p_[1]] += 1
                    if p_[0] != p_[1]:
                        e[p_[1]][p_[0]] += 1
                    brow.append(trunc(bottom(e), 2))
                E11 = [[trunc(Ef[i][j] + Ephi[i][j], 2) for j in range(10)] + [trunc(-2 * Wlow[SP[i][0]][SP[i][1]], 2)] for i in range(10)]
                E11.append(brow + [-v2f])
                Ef = E11
                rep['trE_field_gravscalar'] = tosp(trunc(sum(Ef[i][i] for i in range(11)) + 6 * Rsf + v2f, 2)) == 0   # guard (review, 9 Nov)
                Nc = [[Nc[mu][i] + [sp.Integer(0)] for i in range(10)] + [[sp.Integer(0)] * 11] for mu in idx]   # the scalar's line is flat
        else:
            raise ValueError(name)
        a0 = transport(0, None, Nc, d, 4, init=meye(d))
        if name == 'vector':
            # a_0's x^2 part as c * R^a_{mu b nu} x^mu x^nu: fit c on the first nonzero component, verify on all (P4; see the smoke record)
            basis = [[trunc(sum(R[(a, m, b, n)] * X[m] * X[n] for m in idx for n in idx), 2) for b in idx] for a in idx]
            c = None
            for a in idx:
                for b in idx:
                    if basis[a][b] != 0 and c is None:
                        c = sp.cancel(tosp(homog(a0[a][b], 2)) / tosp(basis[a][b]))
            rep['a0_x2_coeff'] = c
            rep['a0_x2_vector'] = c is not None and c.is_Rational and all(tosp(homog(a0[a][b], 2) - c * basis[a][b]) == 0 for a in idx for b in idx)
        f0 = mmap(lambda e: trunc(Dhalf * e, 4), a0)
        rhs1 = madd(box_bundle(f0, Nc, ginv, sqrtg, isqrtg, 2), mmul(Ef, f0, 2))
        rhs1 = mmap(lambda e: trunc(iDhalf * e, 2), rhs1)
        a1 = transport(1, rhs1, Nc, d, 2)
        f1 = mmap(lambda e: trunc(Dhalf * e, 2), a1)
        rhs2 = madd(box_bundle(f1, Nc, ginv, sqrtg, isqrtg, 0), mmul(Ef, f1, 0))
        a2_0 = mmap(at0, mscale(rhs2, sp.Rational(1, 2)))
        rep[f'a1_{name}'] = mtr(mmap(at0, a1))
        rep[f'a2_{name}'] = mtr(a2_0)
        # method M at the same point: (3.40) = (1/24) Tr Omega^2 + 1/4 Tr(E(0) + R/6)^2 + (d/120)(Ric^2 - R^2/3), with
        # Omega_{mn} = [D_m, D_n] from the connection's curvature at 0: Omega = dN_n/dx_m - dN_m/dx_n + N_m N_n - N_n N_m
        E0 = mmap(at0, Ef)
        rep[f'boxTrE_{name}'] = sp.expand(at0(sum(pd(pd(sum(Ef[i][i] for i in range(d)), m), m) for m in idx)))   # Box(Tr E)(0)
        rep[f'trE_{name}'] = mtr(E0)
        Om2 = sp.Integer(0)
        for m in idx:
            for n in idx:
                Omn = madd(madd(mmap(lambda e, m=m: at0(pd(e, m)), Nc[n]), mmap(lambda e, n=n: at0(pd(e, n)), Nc[m]), -1),
                           madd(mmul(mmap(at0, Nc[m]), mmap(at0, Nc[n]), 0), mmul(mmap(at0, Nc[n]), mmap(at0, Nc[m]), 0), -1))
                Om2 += mtr(mmul(Omn, Omn, 0))
        Xm = madd(E0, mscale(meye(d), Rs / 6))
        rep[f'hv340_{name}'] = sp.expand(sp.Rational(1, 24) * Om2 + sp.Rational(1, 4) * mtr(mmul(Xm, Xm, 0)) + sp.Rational(d, 120) * (Ric2 - Rs**2 / 3))
        if name == 'gravscalar':        # Q1
            rep['Q1_E0_equals_M_gravscalar'] = bool(sp.Matrix(E0).applyfunc(sp.expand) == E_M.applyfunc(sp.expand)) and sp.expand(Om2 - Om2_M) == 0
            rep['hvM_gravscalar'] = GM.hv340(11, Om2_M, E_M, Rs, Ric2)
        if name == 'ghostscalar':
            rep['Q1_E0_equals_M_ghostscalar'] = bool(sp.Matrix(E0).applyfunc(sp.expand) == Eg_M.applyfunc(sp.expand)) and sp.expand(Om2 + Riem2) == 0
            rep['hvM_ghostscalar'] = GM.hv340(4, -Riem2, Eg_M, Rs, Ric2)
    return rep


def fit(reps, key):
    c1, c2, c3 = sp.symbols('c_R2 c_Ric2 c_Riem2')
    eqs = [sp.Eq(c1 * r['inv'][0]**2 + c2 * r['inv'][1] + c3 * r['inv'][2], r[key]) for r in reps]
    sol = sp.solve(eqs[:3], [c1, c2, c3], dict=True)
    if not sol:
        say("FIT-FAILED: the first three backgrounds are degenerate in (R^2, Ric^2, Riem^2)"); sys.exit(1)
    sol = sol[0]
    resid = [sp.simplify(e.lhs.subs(sol) - e.rhs) for e in eqs[3:]]
    return sp.factor(sol[c1]), sp.factor(sol[c2]), sp.factor(sol[c3]), all(r == 0 for r in resid)




def _worker(a):
    return run_background(*a)


def main():
    global out, args
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=4)
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--seed', type=int, default=20261104)
    ap.add_argument('--bundles', default='scalar,vector', help="scalar,vector (step 1), graviton,ghost (step 2; 8 Nov) or gravscalar,ghostscalar (step 3a; 9 Nov)")
    ap.add_argument('--log', default=None)
    ap.add_argument('--backend', default='sympy', choices=['sympy', 'flint'], help='polynomial arithmetic: sympy (reference) or flint (python-flint fmpq_mpoly; same exact rationals, much faster)')
    ap.add_argument('--resume', action='store_true', help='continue a stopped run from <log>.ckpt.json, appending to the log')
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        import os
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            if args.resume and os.path.exists(args.log + '.ckpt.json'):
                out = open(args.log, 'a'); print(f"# --resume: appending to {args.log}, checkpoint {args.log}.ckpt.json")
            elif args.resume:
                print(f"--resume given but no checkpoint {args.log}.ckpt.json exists; nothing to resume (to start over:  del {args.log.replace('/', chr(92))})"); sys.exit(1)
            else:
                print(f"refusing to overwrite the non-empty log {args.log} (the 4 Nov truncation). To continue a stopped run: --resume. To start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        else:
            out = open(args.log, 'w')
    bundles = args.bundles.split(',')
    say(f"# rung3_dewitt_bundle: a_2 on {args.samples} random algebraic curvature tensors, seed {args.seed}, bundles {bundles}, jobs {args.jobs}" + (f", backend {args.backend}" if args.backend != 'sympy' else ''))
    say(version_line())
    rng = random.Random(args.seed)
    Rlist = [riemann_sample(rng) for _ in range(args.samples)]
    NOFIT = ('gravscalar', 'ghostscalar')            # backgrounds with a scalar field: pointwise comparison, no three-structure fit
    need_phi = any(b in NOFIT for b in bundles)
    if need_phi and any(b not in NOFIT for b in bundles):
        print("gravscalar/ghostscalar are run on their own (their backgrounds carry phi-tilde): do not mix with the other bundles"); sys.exit(1)
    rng2 = random.Random(args.seed + 1)            # a separate stream, so the curvature tensors are those of the earlier steps
    def phi_sample():
        return [(m, sp.Rational(rng2.randint(-3, 3), rng2.randint(1, 3))) for m in itertools.product(range(5), repeat=4) if 1 <= sum(m) <= 4]
    Plist = [phi_sample() if need_phi else None for _ in range(args.samples)]
    t1 = time.time()
    # checkpoint: each finished background is appended to <log>.ckpt.json (or rung3_dewitt_bundle.ckpt.json); a restart with the
    # same seed, samples and bundles skips them (8 Nov, the author's ask: stop and resume overnight)
    import json, os
    ckpt_path = (args.log or 'rung3_dewitt_bundle') + '.ckpt.json'
    done = {}
    if os.path.exists(ckpt_path) and args.resume:
        with open(ckpt_path) as f:
            ck = json.load(f)
        if ck.get('seed') == args.seed and ck.get('bundles') == bundles:
            done = {r['label']: r for r in ck['reps']}
            say(f"# resumed from {ckpt_path}: {sorted(done)} already done")
        else:
            say(f"# {ckpt_path} is for another seed/bundles; ignored")

    def save(rep):
        done[rep['label']] = rep
        ser = []
        for r in done.values():
            ser.append({k: (str(v) if not isinstance(v, (bool, str, int)) else v) for k, v in r.items()})
        tmp = ckpt_path + '.tmp'
        with open(tmp, 'w') as f:
            json.dump({'seed': args.seed, 'bundles': bundles, 'reps': ser}, f)
        os.replace(tmp, ckpt_path)

    def restore(r):
        out_ = {}
        for k, v in r.items():
            if isinstance(v, (bool, int)) or k == 'label':
                out_[k] = v
            else:
                val = sp.sympify(v)
                out_[k] = tuple(val) if k == 'inv' else val
        return out_
    todo = [(Rk, f'sample{i}', bundles, args.backend, Plist[i]) for i, Rk in enumerate(Rlist) if f'sample{i}' not in done]
    reps = [restore(r) for r in done.values() if isinstance(r.get('inv'), str)]
    if args.jobs > 1 and todo:
        import multiprocessing as mp
        with mp.Pool(args.jobs) as pool:
            for rep in pool.imap_unordered(_worker, todo):
                reps.append(rep); save(rep); say(f"# {rep['label']} done at {time.time() - t1:.0f} s", flush_ckpt=ckpt_path)
    else:
        for Rk, lab, bl, be, pj in todo:
            rep = run_background(Rk, lab, bl, be, pj); reps.append(rep); save(rep); say(f"# {lab} done at {time.time() - t1:.0f} s", flush_ckpt=ckpt_path)
    reps.sort(key=lambda r: r['label'])
    say(f"# backgrounds done in {time.time() - t1:.0f} s")
    ok = True
    DIM0 = {'scalar': 1, 'vector': 4, 'ghost': 4, 'graviton': 10, 'gravscalar': 11, 'ghostscalar': 4}
    for r in reps:
        checks = {k: v for k, v in r.items() if isinstance(v, bool)}
        say(f"{r['label']}: checks {checks}; Box R(0) = {r['boxR_at_0']}; " + '; '.join(f"Tr a1_{b}(0) = {r[f'a1_{b}']}" for b in bundles))
        Rs_, Ric2_, Riem2_ = r['inv']
        for b in bundles:
            GBv = Riem2_ - 4 * Ric2_ + Rs_**2
            pr = sp.expand(2 * r[f'hv340_{b}'] + sp.Rational(DIM0[b], 180) * GBv)        # a_2 = 2 (3.40) + (d/180) GB, pointwise
            if b in NOFIT:
                bx = r[f'boxTrE_{b}']; pr2 = sp.expand(pr + bx / 6); df = sp.expand(r[f'a2_{b}'] - pr2)
                q4 = sp.expand(r[f'a1_{b}'] - (sp.Rational(DIM0[b], 6) * Rs_ + r[f'trE_{b}'])) == 0
                say(f"   Tr a2_{b}(0) = {r[f'a2_{b}']}   Q2 predicted 2*(3.40) + (d/180)*GB + (1/6)*Box(Tr E) = {pr2}   diff {df}")
                say(f"      [Box(Tr E)(0) = {bx}; without the Box E term the diff would be {sp.expand(r[f'a2_{b}'] - pr)}; (3.40) here = {r[f'hv340_{b}']}, by the sibling's build_E = {r[f'hvM_{b}']}; Q4 a_1: {'MATCH' if q4 else 'MISMATCH'}]")
                ok &= (df == 0) and q4 and sp.expand(r[f'hv340_{b}'] - r[f'hvM_{b}']) == 0
                continue
            say(f"   Tr a2_{b}(0) = {r[f'a2_{b}']}   predicted 2*(3.40) + (d/180)*GB = {pr}   diff {sp.expand(r[f'a2_{b}'] - pr)}")
        ok &= all(checks.values())
        if 'gravscalar' in bundles and 'ghostscalar' in bundles:
            GBv = Riem2_ - 4 * Ric2_ + Rs_**2
            totP = sp.expand((r['a2_gravscalar'] - 2 * r['a2_ghostscalar']) / 2)
            totM = sp.expand(r['hvM_gravscalar'] - 2 * r['hvM_ghostscalar'])
            td = sp.expand(sp.Rational(11 - 8, 360) * GBv + (r['boxTrE_gravscalar'] - 2 * r['boxTrE_ghostscalar']) / 12)
            say(f"   total = (h, phi) - 2 ghost [v^2 = {r['scalar_bg'][0]}, tr W = {r['scalar_bg'][1]}]: method P (1/2) Tr a_2 = {totP}; method M (3.40) = {totM}; P - M = {sp.expand(totP - totM)} = (3/360) GB + (1/12) Box(Tr E_h - 2 Tr E_gh) = {td}: {'MATCH' if sp.expand(totP - totM - td) == 0 else 'MISMATCH'}")
            ok &= sp.expand(totP - totM - td) == 0
    if need_phi:
        say(f"\nQ1 (E(0), Omega^2 = the sibling script's, from the printed (4.31)-(4.33)): {'PASS' if all(r.get('Q1_E0_equals_M_gravscalar', True) and r.get('Q1_E0_equals_M_ghostscalar', True) for r in reps) else 'FAIL'};  Q2, Q4 per background above")
        say(f"\nstatus={'PASS' if ok else 'FAIL'}"); say(f"END rung3_dewitt_bundle {time.time() - t0:.0f} s"); sys.exit(0)
    if len(reps) < 3:
        say("fewer than 3 backgrounds: no fit"); say(f"END rung3_dewitt_bundle {time.time() - t0:.0f} s"); sys.exit(0)
    # ---- fits
    # method M's pre-GB coefficients, COMPUTED here from (3.40) on the same backgrounds (9 Nov: the 8 Nov table was typed by hand
    # from a log and the graviton row was wrong -- status=FAIL on a run whose physics passed; rule 11 record in the header)
    def hv340_coeffs(bname):
        rows = [(r[f'hv340_{bname}'],) + tuple(r['inv']) for r in reps]
        c1, c2, c3 = sp.symbols('m_R2 m_Ric2 m_Riem2')
        eqs = [sp.Eq(c1 * a**2 + c2 * b + c3 * c, v) for v, a, b, c in rows]
        sol = sp.solve(eqs[:3], [c1, c2, c3], dict=True)[0]
        return sp.factor(sol[c1]), sp.factor(sol[c2]), sp.factor(sol[c3])
    DIM = {'scalar': 1, 'vector': 4, 'ghost': 4, 'graviton': 10}
    M = {b: hv340_coeffs(b) + (DIM[b],) for b in bundles}
    predP = {b: tuple(sp.expand(2 * m + sp.Rational(DIM[b], 180) * k) for m, k in zip(M[b][:3], (1, -4, 1))) for b in bundles}
    res = {}
    for b in bundles:
        cR2, cRic2, cRiem2, fit_ok = fit(reps, f'a2_{b}')
        res[b] = (cR2, cRic2, cRiem2)
        gbR2, gbRic2 = sp.expand(cR2 - cRiem2), sp.expand(cRic2 + 4 * cRiem2)
        say(f"\n## {b}: Tr a_2(0) = {cR2} R^2 + {cRic2} Ric^2 + {cRiem2} Riem^2   (residuals zero: {fit_ok});  after GB: {gbR2} R^2 + {gbRic2} Ric^2")
        mR2, mRic2, mRiem2, d = M[b]
        say(f"   method M (3.40) on the same backgrounds, pre-GB: {mR2} R^2 + {mRic2} Ric^2 + {mRiem2} Riem^2")
        diff = [sp.expand(c / 2 - m) for c, m in zip((cR2, cRic2, cRiem2), (mR2, mRic2, mRiem2))]
        want = [sp.Rational(d, 360) * k for k in (1, -4, 1)]
        p3 = all(sp.expand(x - y) == 0 for x, y in zip(diff, want))
        say(f"   (1/2) Tr a_2 - (3.40) = {diff}   predicted (d/360)(1,-4,1) = {want}: {'MATCH' if p3 else 'MISMATCH'}")
        say(f"   after GB, (1/2) Tr a_2 = {sp.expand(gbR2 / 2)} R^2 + {sp.expand(gbRic2 / 2)} Ric^2  vs method M {sp.expand(mR2 - mRiem2)} R^2 + {sp.expand(mRic2 + 4 * mRiem2)} Ric^2")
        pP = all(sp.expand(c - p) == 0 for c, p in zip((cR2, cRic2, cRiem2), predP[b]))
        say(f"   P2 (a_2 = 2M + (d/180)(1,-4,1), M computed here): {'MATCH' if pP else 'MISMATCH'}")
        if b == 'scalar':
            p1s = all(sp.expand(c - p) == 0 for c, p in zip((cR2, cRic2, cRiem2), (xi**2 / 2 - xi / 6 + sp.Rational(1, 72), -sp.Rational(1, 180), sp.Rational(1, 180))))
            say(f"   P1 (rung 1's [a_2] with E = -xi R): {'MATCH' if p1s else 'MISMATCH'}"); ok &= p1s
        ok &= fit_ok and p3 and pP
    if 'graviton' in bundles and 'ghost' in bundles:
        pg = [sp.expand(v - 2 * s_) for v, s_ in zip(res['graviton'], res['ghost'])]
        say(f"\n## pure gravity = graviton - 2 ghost: Tr a_2 = {pg};  after GB /2: {sp.expand((pg[0] - pg[2]) / 2)} R^2 + {sp.expand((pg[1] + 4 * pg[2]) / 2)} Ric^2  (5.24: 1/120, 7/20)")
    if 'vector' in bundles and 'scalar' in bundles:
        mx = [sp.expand(v - 2 * s.subs(xi, 0)) for v, s in zip(res['vector'], res['scalar'])]
        say(f"\n## maxwell = vector - 2 scalar(0): Tr a_2 = {mx}   predicted (1/180)(-25, 88, -13) = {[sp.Rational(k, 180) for k in (-25, 88, -13)]}; after GB /2: {sp.expand((mx[0] - mx[2]) / 2)} R^2 + {sp.expand((mx[1] + 4 * mx[2]) / 2)} Ric^2 (method M: -1/30, 1/10)")
    if 'vector' in bundles:
        p4 = all(r.get('a0_x2_vector', False) and r.get('a0_x2_coeff') == sp.Rational(1, 6) for r in reps)
        say(f"P4 {'PASS' if p4 else 'FAIL'}: a_0's x^2 part is c R^a_{{mu b nu}} x^mu x^nu with c = {[r.get('a0_x2_coeff') for r in reps]} (re-scored: +1/6; the header said -1/6)")
    say(f"\nstatus={'PASS' if ok else 'FAIL'}")
    say(f"END rung3_dewitt_bundle {time.time() - t0:.0f} s")


if __name__ == '__main__':      # the 4 Nov lesson: on Windows a pool re-imports this file, and an unguarded body forks forever
    main()
