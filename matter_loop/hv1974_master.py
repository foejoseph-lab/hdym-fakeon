#!/usr/bin/env python3
"""
Rung 3, step 1 (method M): 't Hooft & Veltman 1974 eq. (3.40) as code -- the one-loop counter-Lagrangian of any second-order
operator on a bundle over a four-dimensional background, evaluated pointwise on random algebraic curvature tensors and fitted to
R^2, Ric^2, Riem^2, then reduced by Gauss-Bonnet (their appendix B).  Summary and plan: ../notes/2026-11-04_rung3_thooft_veltman_summary.md.

The formula (real fields, eq. 3.39-3.40; 1/eps = 1/(8 pi^2 (n-4))):
    L  = sqrt(g) ( -1/2 d phi g d phi + phi N^mu d_mu phi + 1/2 phi M phi )
    dL = (sqrt(g)/eps) Tr { (1/24) Y^{mu nu} Y_{mu nu} + 1/4 (M - N^mu N_mu - R/6)^2 + (1/120) (Ric^2 - R^2/3) 1 },
    Y_{mu nu} = D_mu N_nu - D_nu N_mu + N_mu N_nu - N_nu N_mu.
In covariant language (the bundle connection absorbed into D): Y = the bundle curvature Omega_{mu nu} = [D_mu, D_nu] and
M - N.N -> the endomorphism E of the operator -(D^2 + E).  CONVENTION TRAP (3 Nov; the justification corrected by the review of
4 Nov): their curvature scalar has the OPPOSITE sign to the usual (MTW/Gilkey) one.  The evidence is p. 79: the coupling that kills the
cross term is +(1/12) R_HV phi^2, i.e. M = R_HV/6, while in standard sign the conformal coupling is E = -R_std/6; hence R_HV = -R_std
and their (M - R/6)^2 is the standard (E + R/6)^2.  (The fragment of (A.12) quoted in an earlier draft decides nothing: the standard
Riemann contains -d_beta Gamma^mu_{nu alpha} too.)  In standard sign (3.40) is one half of Gilkey's a_2 = (1/360) Tr(180 E^2 + 60 R E
+ 5 R^2 - 2 Ric^2 + 2 Riem^2 + 30 Omega^2) MODULO the Gauss-Bonnet density: the difference is -(dim/360)(Riem^2 - 4 Ric^2 + R^2),
since the paper never fixes the Riem^2 coefficient (appendix B).  So the pre-GB rows of the log are NOT a_2/2; only the post-GB
numbers are.  This script works in STANDARD curvature sign throughout.

Bundles implemented (Feynman gauge where a gauge is needed; ghosts as separate scalar entries with the FP minus sign):
    scalar(xi):  dim 1, Omega = 0, E = -xi R                (L = -1/2 (d phi)^2 - 1/2 xi R phi^2 in standard sign)
    vector:      dim 4, Omega^a_{b mu nu} = R^a_{b mu nu}, E^a_b = -R^a_b   (-1/4 F^2 - 1/2 (D.A)^2 -> -(D^2 delta - Ric))
    maxwell:     vector - 2 * scalar(0)                       (the two FP ghosts of the Lorenz gauge)

PREDICTIONS (rule 25; written 4 Nov 19:30 before the first run; nothing fetched):
  P1: scalar(0) -> (1/240) R^2 + (1/120) Ric^2 after Gauss-Bonnet [v, their 5.23]; scalar(xi) -> R^2 coefficient
      1/4 (xi - 1/6)^2 - 1/360, Ric^2 coefficient 1/120 [rung 1's six-way number, results note of 31 Oct, line "p^4"; the paper
      gives only xi = 0 and the xi = 1/6 statement].
  P2 [obvious to 1974]: maxwell is conformally invariant in four dimensions, so after Gauss-Bonnet its counterterm is proportional to
      Ric^2 - R^2/3 (= C^2/2 mod GB): the R^2 coefficient is exactly -1/3 of the Ric^2 coefficient.
  P3 [memory, marked]: the Ric^2 coefficient of maxwell is 1/10 = 12 x the scalar's 1/120 (the 1:12 ratio of the Weyl-squared
      coefficients of a conformal scalar and a Maxwell field; Duff's 1977 trace-anomaly table is where I remember it from, which is
      three years ahead of the rung -- marked as memory, not derived).  Derived here by hand before the run, with standard signs and
      GB: (1/24) Tr Omega^2 = -(1/24) Riem^2 -> -(1/6) Ric^2 + (1/24) R^2; (1/4) Tr (E + R/6)^2 = (1/4) Ric^2 - R^2/18;
      4 x (1/120)(Ric^2 - R^2/3); minus 2 scalars: total (1/10) Ric^2 - (1/30) R^2.  The run decides whether the hand sum is right.
  P4 (reported, not a prediction -- the review of 4 Nov pointed out it is circular with gb()): the Riem^2 coefficient is nonzero for
      vector (from Tr Omega^2 = -Riem^2) and zero for scalar.
  On the "residuals all zero" line: hv340 returns Tr E^2, R Tr E, Tr Omega^2, R^2, Ric^2, which for these bundles are identically
      Ric^2, -R^2, -Riem^2, R^2, Ric^2 on any algebraic curvature tensor, so three backgrounds determine the coefficients exactly and
      the residuals on the rest vanish by construction: that line checks the index bookkeeping, not the formula (review, 4 Nov).
      The check of the formula is P1 (rung 1, six regulators) and P2 (conformal invariance, which owes nothing to memory).
  Single-method numbers are not checks (rule 12): everything here is method M until the DeWitt recursion on the bundle (method P,
  next session) reproduces it.  P1 is the exception, being rung 1's number six ways over.

Run:  python3 hv1974_master.py --samples 8 --log results/hv1974_master.log        [seconds]
"""
import sys, time, argparse, itertools, random
import sympy as sp

ap = argparse.ArgumentParser()
ap.add_argument('--samples', type=int, default=8)
ap.add_argument('--seed', type=int, default=20261104)
ap.add_argument('--log', default=None)
args = ap.parse_args()
t0 = time.time()
out = open(args.log, 'w') if args.log else None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


N = 4
idx = range(N)
xi = sp.Symbol('xi')


def riemann_sample(rng, nmat=8):
    """Random rational algebraic curvature tensor (sum of Kulkarni-Nomizu squares); rung1_dewitt.py's generator, copied.
    Standard sign: R_{abcd} = A_ac A_bd - A_ad A_bc has R_{abab} > 0 for the unit sphere."""
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
    Ric = {(b, d): sum(R[(a, b, a, d)] for a in idx) for b in idx for d in idx}     # R_{bd} = R^a_{bad}, Euclidean frame
    Rs = sp.expand(sum(Ric[(a, a)] for a in idx))
    Ric2 = sp.expand(sum(Ric[(a, b)] ** 2 for a in idx for b in idx))
    Riem2 = sp.expand(sum(v ** 2 for v in R.values()))
    return Ric, Rs, Ric2, Riem2


def tr(A):
    return sum(A[i][i] for i in range(len(A)))


def matmul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def hv340(dim, Omega, E, Rs, Ric2):
    """(3.40) in standard curvature sign at a point: Omega[mu][nu] and E are dim x dim matrices (lists), Rs = R, Ric2 = Ric^2."""
    I = [[sp.Integer(i == j) for j in range(dim)] for i in range(dim)]
    Y2 = sum(tr(matmul(Omega[m][n], Omega[m][n])) for m in idx for n in idx)             # Euclidean frame: up = down
    X = [[E[i][j] + Rs / 6 * I[i][j] for j in range(dim)] for i in range(dim)]
    X2 = tr(matmul(X, X))
    return sp.expand(sp.Rational(1, 24) * Y2 + sp.Rational(1, 4) * X2 + sp.Rational(1, 120) * dim * (Ric2 - Rs ** 2 / 3))


def bundle(name, R, Ric, Rs, Ric2):
    if name.startswith('scalar'):
        x = sp.sympify(name[len('scalar('):-1]) if '(' in name else sp.Integer(0)
        Om = [[[[sp.Integer(0)]] for _ in idx] for _ in idx]
        return hv340(1, Om, [[-x * Rs]], Rs, Ric2)
    if name == 'vector':
        Om = [[[[R[(a, b, m, n)] for b in idx] for a in idx] for n in idx] for m in idx]   # Omega^a_b mu nu = R^a_{b mu nu}
        E = [[-Ric[(a, b)] for b in idx] for a in idx]
        return hv340(4, Om, E, Rs, Ric2)
    if name == 'maxwell':
        return sp.expand(bundle('vector', R, Ric, Rs, Ric2) - 2 * bundle('scalar(0)', R, Ric, Rs, Ric2))
    raise ValueError(name)


def fit(rows):
    """rows: list of (value, R, Ric2, Riem2); exact least-squares-free solve: pick 3 independent rows, verify on the rest."""
    c1, c2, c3 = sp.symbols('c_R2 c_Ric2 c_Riem2')
    eqs = [sp.Eq(c1 * Rs ** 2 + c2 * Ric2 + c3 * Riem2, v) for v, Rs, Ric2, Riem2 in rows]
    sol = sp.solve(eqs[:3], [c1, c2, c3], dict=True)
    assert sol, "first three backgrounds degenerate"
    sol = sol[0]
    resid = [sp.simplify(e.lhs.subs(sol) - e.rhs) for e in eqs[3:]]
    return sol[c1], sol[c2], sol[c3], resid


def gb(cR2, cRic2, cRiem2):
    """Riem^2 -> 4 Ric^2 - R^2 (appendix B, four dimensions only)."""
    return sp.expand(cR2 - cRiem2), sp.expand(cRic2 + 4 * cRiem2)


say(f"# hv1974_master: (3.40) on {args.samples} random algebraic curvature tensors, seed {args.seed}; standard curvature sign")
rng = random.Random(args.seed)
backs = [riemann_sample(rng) for _ in range(args.samples)]
ok = True
results = {}
for name in ['scalar(0)', 'scalar(1/6)', 'scalar(xi)', 'vector', 'maxwell']:
    rows = []
    for R in backs:
        Ric, Rs, Ric2, Riem2 = invariants(R)
        rows.append((bundle(name, R, Ric, Rs, Ric2), Rs, Ric2, Riem2))
    cR2, cRic2, cRiem2, resid = fit(rows)
    fit_ok = all(r == 0 for r in resid)
    gR2, gRic2 = gb(cR2, cRic2, cRiem2)
    results[name] = (cR2, cRic2, cRiem2, gR2, gRic2)
    say(f"\n## {name}: fit R^2 {cR2}  Ric^2 {cRic2}  Riem^2 {cRiem2}   (residuals on the other {len(resid)} backgrounds all zero: {fit_ok})")
    say(f"   after Gauss-Bonnet:  R^2 {gR2}   Ric^2 {gRic2}")
    ok &= fit_ok

# ---- P1
s0 = results['scalar(0)']; sx = results['scalar(xi)']
p1 = (s0[3] == sp.Rational(1, 240) and s0[4] == sp.Rational(1, 120)
      and sp.expand(sx[3] - (sp.Rational(1, 4) * (xi - sp.Rational(1, 6)) ** 2 - sp.Rational(1, 360))) == 0 and sx[4] == sp.Rational(1, 120))
say(f"\nP1 {'PASS' if p1 else 'FAIL'}: scalar(0) = R^2/240 + Ric^2/120 (their 5.23); scalar(xi) R^2 = 1/4 (xi-1/6)^2 - 1/360, Ric^2 = 1/120 (rung 1)")
say(f"   scalar(1/6) = {results['scalar(1/6)'][3]} R^2 + {results['scalar(1/6)'][4]} Ric^2  (conformal: R^2/Ric^2 = {results['scalar(1/6)'][3] / results['scalar(1/6)'][4]})")
# ---- P2, P3
mx = results['maxwell']
p2 = sp.expand(mx[3] + mx[4] / 3) == 0
p3 = mx[4] == sp.Rational(1, 10) and mx[3] == -sp.Rational(1, 30)
say(f"P2 {'PASS' if p2 else 'FAIL'}: maxwell conformal, R^2 = -Ric^2/3   (R^2 {mx[3]}, Ric^2 {mx[4]})")
say(f"P3 {'PASS' if p3 else 'FAIL'}: maxwell Ric^2 = 1/10 = 12 x 1/120, R^2 = -1/30 [memory + hand sum; single method until P]")
# ---- P4
v = results['vector']
p4 = (v[2] != 0 and s0[2] == 0 and sp.expand(v[3] - v[0] + v[2]) == 0 and sp.expand(v[4] - v[1] - 4 * v[2]) == 0)
say(f"P4 {'PASS' if p4 else 'FAIL'}: Riem^2 coefficient vector {v[2]}, scalar {s0[2]}; GB shift as stated")
ok &= p1 and p2 and p3 and p4
say(f"\nstatus={'PASS' if ok else 'FAIL'}  (P1 {p1}, P2 {p2}, P3 {p3}, P4 {p4})")
say(f"END hv1974_master {time.time() - t0:.0f} s")
