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
(8 Nov: the body is under a main guard -- on Windows multiprocessing re-imports the file and the unguarded 4 Nov version forked
workers forever, each truncating the log; the log is now refused if non-empty, and a line is printed per finished background.)
"""
import time, argparse, itertools, random, sys
import sympy as sp

out = None


def say(*a):
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


# ----------------------------------------------------------------- rung1_dewitt.py, verbatim
def trunc(expr, deg=DEG):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(expr, *X)
    return sp.Add(*[c * sp.Mul(*[x**e for x, e in zip(X, mon)]) for mon, c in P.terms() if sum(mon) <= deg])


def homog(expr, d):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(expr, *X)
    return sp.Add(*[c * sp.Mul(*[x**e for x, e in zip(X, mon)]) for mon, c in P.terms() if sum(mon) == d])


def at0(expr):
    return sp.expand(expr).subs({x: 0 for x in X})


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
    dg = [[[trunc(sp.diff(g[a][b], X[c])) for c in idx] for b in idx] for a in idx]
    G = [[[None] * N for _ in idx] for _ in idx]
    for mu in idx:
        for a in idx:
            for b in idx:
                G[mu][a][b] = trunc(sp.Rational(1, 2) * sum(ginv[mu][l] * (dg[l][b][a] + dg[l][a][b] - dg[a][b][l]) for l in idx))
    return G


def riemann_field(G, g):
    Rup = {}
    for r, s, m, n in itertools.product(idx, repeat=4):
        Rup[(r, s, m, n)] = trunc(sp.diff(G[r][n][s], X[m]) - sp.diff(G[r][m][s], X[n])
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
    return sp.expand(sum(A[i][i] for i in range(len(A))))


def box_bundle(f, Nc, ginv, sqrtg, isqrtg, deg):
    """Box_N f = (1/sqrt g) d_mu (sqrt g g^{mu nu} V_nu) + g^{mu nu} N_mu V_nu,  V_nu = d_nu f + N_nu f; f, N_mu d x d matrices."""
    V = [madd(mmap(lambda e: sp.diff(e, X[nu]), f), mmul(Nc[nu], f, deg + 1)) for nu in idx]
    d = len(f)
    tot = mzero(d)
    for mu in idx:
        inner = mzero(d)
        for nu in idx:
            inner = madd(inner, mmap(lambda e, mu=mu, nu=nu: trunc(sqrtg * ginv[mu][nu] * e, deg + 1), V[nu]))
        tot = madd(tot, mmap(lambda e, mu=mu: sp.diff(e, X[mu]), inner))
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


def run_background(R, label, bundles):
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
    rep['boxR_at_0'] = sp.expand(sum(sp.diff(Rsf, X[m], X[m]) for m in idx).subs({x: 0 for x in X}))
    rep['Dhalf_x2'] = sp.expand(homog(Dhalf, 2) - sp.Rational(1, 12) * sum(Ric[(a, b)] * X[a] * X[b] for a in idx for b in idx)) == 0
    # Ricci field with one index up: R^a_b(x) = g^{ac} R_{cb}, R_{cb} = R^r_{c r b}
    Ricf = {(c, b): trunc(sum(Rup[(r, c, r, b)] for r in idx), 2) for c in idx for b in idx}
    Ricup = [[trunc(sum(ginv[a][c] * Ricf[(c, b)] for c in idx), 2) for b in idx] for a in idx]
    for name in bundles:
        if name == 'scalar':
            d = 1; Nc = [mzero(1) for _ in idx]; Ef = [[-xi * Rsf]]
        elif name == 'vector':
            d = 4; Nc = [[[G[a][mu][b] for b in idx] for a in idx] for mu in idx]          # (N_mu)^a_b = Gamma^a_{mu b}
            Ef = [[-Ricup[a][b] for b in idx] for a in idx]
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
                        c = sp.cancel(homog(a0[a][b], 2) / basis[a][b])
            rep['a0_x2_coeff'] = c
            rep['a0_x2_vector'] = c is not None and c.is_Rational and all(sp.expand(homog(a0[a][b], 2) - c * basis[a][b]) == 0 for a in idx for b in idx)
        f0 = mmap(lambda e: trunc(Dhalf * e, 4), a0)
        rhs1 = madd(box_bundle(f0, Nc, ginv, sqrtg, isqrtg, 2), mmul(Ef, f0, 2))
        rhs1 = mmap(lambda e: trunc(iDhalf * e, 2), rhs1)
        a1 = transport(1, rhs1, Nc, d, 2)
        f1 = mmap(lambda e: trunc(Dhalf * e, 2), a1)
        rhs2 = madd(box_bundle(f1, Nc, ginv, sqrtg, isqrtg, 0), mmul(Ef, f1, 0))
        a2_0 = mmap(at0, mscale(rhs2, sp.Rational(1, 2)))
        rep[f'a1_{name}'] = mtr(mmap(at0, a1))
        rep[f'a2_{name}'] = mtr(a2_0)
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
    ap.add_argument('--bundles', default='scalar,vector')
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        import os
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log} (the 4 Nov truncation); move it first"); sys.exit(1)
        out = open(args.log, 'w')
    bundles = args.bundles.split(',')
    say(f"# rung3_dewitt_bundle: a_2 on {args.samples} random algebraic curvature tensors, seed {args.seed}, bundles {bundles}, jobs {args.jobs}")
    rng = random.Random(args.seed)
    Rlist = [riemann_sample(rng) for _ in range(args.samples)]
    t1 = time.time()
    reps = []
    if args.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(args.jobs) as pool:
            for rep in pool.imap_unordered(_worker, [(Rk, f'sample{i}', bundles) for i, Rk in enumerate(Rlist)]):
                reps.append(rep); say(f"# {rep['label']} done at {time.time() - t1:.0f} s")
        reps.sort(key=lambda r: r['label'])
    else:
        for i, Rk in enumerate(Rlist):
            reps.append(run_background(Rk, f'sample{i}', bundles)); say(f"# sample{i} done at {time.time() - t1:.0f} s")
    say(f"# backgrounds done in {time.time() - t1:.0f} s")
    ok = True
    for r in reps:
        checks = {k: v for k, v in r.items() if isinstance(v, bool)}
        say(f"{r['label']}: checks {checks}; Box R(0) = {r['boxR_at_0']}; " + '; '.join(f"Tr a1_{b}(0) = {r[f'a1_{b}']}" for b in bundles))
        Rs_, Ric2_, Riem2_ = r['inv']
        for b in bundles:
            pr = {'scalar': (xi**2 / 2 - xi / 6 + sp.Rational(1, 72)) * Rs_**2 - Ric2_ / 180 + Riem2_ / 180,
                  'vector': -Rs_**2 / 9 + sp.Rational(43, 90) * Ric2_ - sp.Rational(11, 180) * Riem2_}[b]
            say(f"   Tr a2_{b}(0) = {r[f'a2_{b}']}   predicted {sp.expand(pr)}   diff {sp.expand(r[f'a2_{b}'] - pr)}")
        ok &= all(checks.values())
    if len(reps) < 3:
        say("fewer than 3 backgrounds: no fit"); say(f"END rung3_dewitt_bundle {time.time() - t0:.0f} s"); sys.exit(0)
    # ---- fits
    M = {'scalar': (xi**2 / 4 - xi / 12 + sp.Rational(1, 240), sp.Rational(1, 120), sp.Integer(0), 1),
         'vector': (-sp.Rational(1, 15), sp.Rational(17, 60), -sp.Rational(1, 24), 4)}            # hv1974_master.py pre-GB, and d
    predP = {'scalar': (xi**2 / 2 - xi / 6 + sp.Rational(1, 72), -sp.Rational(1, 180), sp.Rational(1, 180)),
             'vector': (-sp.Rational(1, 9), sp.Rational(43, 90), -sp.Rational(11, 180))}
    res = {}
    for b in bundles:
        cR2, cRic2, cRiem2, fit_ok = fit(reps, f'a2_{b}')
        res[b] = (cR2, cRic2, cRiem2)
        gbR2, gbRic2 = sp.expand(cR2 - cRiem2), sp.expand(cRic2 + 4 * cRiem2)
        say(f"\n## {b}: Tr a_2(0) = {cR2} R^2 + {cRic2} Ric^2 + {cRiem2} Riem^2   (residuals zero: {fit_ok});  after GB: {gbR2} R^2 + {gbRic2} Ric^2")
        mR2, mRic2, mRiem2, d = M[b]
        diff = [sp.expand(c / 2 - m) for c, m in zip((cR2, cRic2, cRiem2), (mR2, mRic2, mRiem2))]
        want = [sp.Rational(d, 360) * k for k in (1, -4, 1)]
        p3 = all(sp.expand(x - y) == 0 for x, y in zip(diff, want))
        say(f"   (1/2) Tr a_2 - (3.40) = {diff}   predicted (d/360)(1,-4,1) = {want}: {'MATCH' if p3 else 'MISMATCH'}")
        say(f"   after GB, (1/2) Tr a_2 = {sp.expand(gbR2 / 2)} R^2 + {sp.expand(gbRic2 / 2)} Ric^2  vs method M {sp.expand(mR2 - mRiem2)} R^2 + {sp.expand(mRic2 + 4 * mRiem2)} Ric^2")
        pP = all(sp.expand(c - p) == 0 for c, p in zip((cR2, cRic2, cRiem2), predP[b]))
        say(f"   {'P1' if b == 'scalar' else 'P2'} ({'rung 1' if b == 'scalar' else 'memory, Gilkey'}): {'MATCH' if pP else 'MISMATCH'}")
        ok &= fit_ok and p3 and (pP if b == 'scalar' else True)
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
