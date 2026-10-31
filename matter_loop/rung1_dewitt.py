#!/usr/bin/env python3
"""
Rung 1, technique (b2): DeWitt's 1965 covariant recursion for the heat-kernel coefficients, run in Riemann
normal coordinates on a generic background (no expansion in h).  Setup and predictions:
../notes/2026-10-31_rung1_b2_dewitt_recursion_setup.md.

Operator D = -(g^{mu nu} nabla_mu nabla_nu + E) on a real scalar, E = -xi R - m^2, kept as one constant symbol
E.  Under the ansatz below R(x) - R(0) has no x^1 part and a traceless x^2 part (the code reports Box R(0) and
the number of nonzero x^2 coefficients of R(x)); only the trace enters [a_2], through Box E = -xi Box R, and it
is zero on these backgrounds, so treating E as constant is exact here.
Ansatz K = (4 pi t)^{-n/2} Delta^{1/2} exp(-sigma/2t) sum_k t^k a_k, recursion
    sigma^{;mu} nabla_mu a_k + k a_k = Delta^{-1/2} (Box + E)(Delta^{1/2} a_{k-1}),  a_0 = 1.
In normal coordinates about the origin: sigma = x^2/2, sigma^{;mu} = x^mu, Delta = 1/sqrt(g), so the left
side is (x.d + k) a_k and a_k = sum_d P_d/(d+k) over the homogeneous parts P_d of the right side.

Input metric (R = +2 on the unit S^2; nabla R = 0 at the origin; the derivative terms of the general
normal-coordinate expansion are dropped, which only changes the total-derivative content of [a_2]):
    g_{mu nu} = delta - (1/3) R_{mu a nu b} x^a x^b + (2/45) R_{a mu b l} R_{c nu d l} x^a x^b x^c x^d + O(x^5).
Checks made on this input before the recursion: Gamma^mu_{ab} x^a x^b = O(x^4); Riemann(0) = input;
d Riemann(0) = 0 (trivially, the metric is even in x; checked anyway); R(0) of the ansatz = contraction of the
input; Delta^{1/2} = 1 + (1/12) Ric_{ab} x^a x^b + O(x^4).  Box R(0) of the ansatz is computed
and fitted as a fourth structure ([a_2] contains it), not assumed to vanish.  Within this family of metrics
Box R(0) is itself quadratic in R, i.e. a combination of R^2, Ric^2, Riem^2, so it cannot be fitted
independently: the curvature-squared coefficients are unambiguous only on a metric with Box R(0) = 0, which is
what the standard x^4 term delivers (checked on every background).  --quartic 0 is kept as the control that
shows what goes wrong without it.

PREDICTIONS (written before the first run; see the setup note section 4):
    [a_1] = E + R/6
    [a_2] = (1/2) E^2 + (1/6) R E + c_R2 R^2 + c_Ric2 Ric^2 + c_Riem2 Riem^2 + c_box Box R   with
    c_Ric2 + 4 c_Riem2 = 1/60  and  c_R2 - c_Riem2 = 1/120   (forced by the six O(h^2) logs);
    expected c_R2 = 1/72, c_Ric2 = -1/180, c_Riem2 = 1/180, c_box = 1/30;  unit S^4: 29/15;  S^2 x S^2: 11/45  (E = 0).
    E is a constant symbol here (m^2 and xi R(0)); the Box E total derivative is not represented and not fitted.
Modes:
    --samples N --jobs J N random rational algebraic curvature tensors (Kulkarni-Nomizu squares), fitted
    --symbolic           the 20 independent components of R_{abcd} are symbols; the fit is an exact identity.
                         NOT RUN (no log); kept as an option, expected to be much slower than the samples.
Run:
    python3 rung1_dewitt.py --samples 6 --jobs 6 --log results/rung1_dewitt_samples.log
"""
import time, argparse, itertools, random, sys
from fractions import Fraction
import sympy as sp

ap = argparse.ArgumentParser()
ap.add_argument('--symbolic', action='store_true')
ap.add_argument('--samples', type=int, default=0)
ap.add_argument('--jobs', type=int, default=1)
ap.add_argument('--seed', type=int, default=20261031)
ap.add_argument('--quartic', default='2/45', help='coefficient of the R.R x^4 term of the normal-coordinate metric')
ap.add_argument('--log', default=None)
args = ap.parse_args()
if not args.symbolic and args.samples == 0:
    args.samples = 4
t0 = time.time()
out = open(args.log, 'w') if args.log else None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


N = 4
X = sp.symbols('x1:5', real=True)
E = sp.Symbol('E')
DEG = 4
idx = range(N)
PAIRS = [(a, b) for a in idx for b in idx if a < b]


# ----------------------------------------------------------------- polynomial helpers (truncate in x)
def trunc(expr, deg=DEG):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(expr, *X)
    return sp.Add(*[c * sp.Mul(*[x**e for x, e in zip(X, mon)])
                    for mon, c in P.terms() if sum(mon) <= deg])


def homog(expr, d):
    expr = sp.expand(expr)
    if expr == 0:
        return sp.Integer(0)
    P = sp.Poly(expr, *X)
    return sp.Add(*[c * sp.Mul(*[x**e for x, e in zip(X, mon)]) for mon, c in P.terms() if sum(mon) == d])


def at0(expr):
    return sp.expand(expr).subs({x: 0 for x in X})


# ----------------------------------------------------------------- the Riemann tensor at the origin
def riemann_symbolic():
    """20 symbols: symmetric 6x6 in pair indices, first Bianchi used to eliminate S[(0,3),(1,2)]."""
    S = {}
    for i, P in enumerate(PAIRS):
        for j, Q in enumerate(PAIRS):
            if i <= j:
                S[(P, Q)] = sp.Symbol(f'R{P[0]+1}{P[1]+1}{Q[0]+1}{Q[1]+1}')
    # R_{0123} + R_{0231} + R_{0312} = 0  ->  S[(01),(23)] - S[(02),(13)] + S[(03),(12)] = 0
    S[((0, 3), (1, 2))] = S[((0, 2), (1, 3))] - S[((0, 1), (2, 3))]

    def comp(a, b, c, d):
        if a == b or c == d:
            return sp.Integer(0)
        s = 1
        if a > b:
            a, b, s = b, a, -s
        if c > d:
            c, d, s = d, c, -s
        P, Q = (a, b), (c, d)
        if PAIRS.index(P) > PAIRS.index(Q):
            P, Q = Q, P
        return s * S[(P, Q)]
    R = {(a, b, c, d): comp(a, b, c, d) for a in idx for b in idx for c in idx for d in idx}
    syms = sorted({v for k, v in S.items() if k != ((0, 3), (1, 2))}, key=lambda s: s.name)
    return R, syms


def riemann_sample(rng, nmat=8):
    """Random rational algebraic curvature tensor: sum of Kulkarni-Nomizu squares of symmetric matrices."""
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


def riemann_named(name):
    R = {k: sp.Integer(0) for k in itertools.product(idx, repeat=4)}
    if name == 'S4':
        for a, b, c, d in itertools.product(idx, repeat=4):
            R[(a, b, c, d)] = sp.Integer((a == c) * (b == d) - (a == d) * (b == c))
    elif name == 'S2xS2':
        blk = lambda a: a // 2
        for a, b, c, d in itertools.product(idx, repeat=4):
            if blk(a) == blk(b) == blk(c) == blk(d):
                R[(a, b, c, d)] = sp.Integer((a == c) * (b == d) - (a == d) * (b == c))
    return R


def invariants(R):
    Ric = {(b, d): sum(R[(a, b, a, d)] for a in idx) for b in idx for d in idx}
    Rs = sum(Ric[(a, a)] for a in idx)
    Ric2 = sum(Ric[(a, b)]**2 for a in idx for b in idx)
    Riem2 = sum(v**2 for v in R.values())
    return Ric, sp.expand(Rs), sp.expand(Ric2), sp.expand(Riem2)


# ----------------------------------------------------------------- geometry from the ansatz metric
def geometry(R):
    g = [[None] * N for _ in idx]
    for mu in idx:
        for nu in idx:
            q2 = -sp.Rational(1, 3) * sum(R[(mu, a, nu, b)] * X[a] * X[b] for a in idx for b in idx)
            q4 = sp.Rational(args.quartic) * sum(R[(a, mu, b, l)] * R[(l, c, nu, d)] * X[a] * X[b] * X[c] * X[d]
                                          for a in idx for b in idx for c in idx for d in idx for l in idx)
            g[mu][nu] = trunc(sp.Integer(mu == nu) + q2 + q4)
    A = [[trunc(g[m][n] - sp.Integer(m == n)) for n in idx] for m in idx]
    A2 = [[homog(A[m][n], 2) for n in idx] for m in idx]
    A2sq = [[trunc(sum(A2[m][k] * A2[k][n] for k in idx)) for n in idx] for m in idx]
    ginv = [[trunc(sp.Integer(m == n) - A[m][n] + A2sq[m][n]) for n in idx] for m in idx]
    L = trunc(sum(A[m][m] for m in idx) - sp.Rational(1, 2) * sum(A2sq[m][m] for m in idx))   # ln det g
    sqrtg = trunc(1 + L / 2 + L**2 / 8)
    isqrtg = trunc(1 - L / 2 + L**2 / 8)
    Dhalf = trunc(1 - L / 4 + L**2 / 32)        # Delta^{1/2} = g^{-1/4}
    iDhalf = trunc(1 + L / 4 + L**2 / 32)
    return g, ginv, sqrtg, isqrtg, Dhalf, iDhalf


def christoffel(g, ginv):
    dg = [[[trunc(sp.diff(g[a][b], X[c])) for c in idx] for b in idx] for a in idx]   # dg[a][b][c] = d_c g_ab
    G = [[[None] * N for _ in idx] for _ in idx]
    for mu in idx:
        for a in idx:
            for b in idx:
                G[mu][a][b] = trunc(sp.Rational(1, 2) * sum(ginv[mu][l] * (dg[l][b][a] + dg[l][a][b] - dg[a][b][l])
                                                            for l in idx))
    return G


def riemann_field(G, g):
    """R^rho_{sigma mu nu} = d_mu G^rho_{nu sigma} - d_nu G^rho_{mu sigma} + G G - G G, then lowered; as
    polynomials correct through degree 2."""
    Rup = {}
    for r, s, m, n in itertools.product(idx, repeat=4):
        Rup[(r, s, m, n)] = trunc(sp.diff(G[r][n][s], X[m]) - sp.diff(G[r][m][s], X[n])
                                  + sum(G[r][m][l] * G[l][n][s] - G[r][n][l] * G[l][m][s] for l in idx), 2)
    Rlow = {}
    for a, s, m, n in itertools.product(idx, repeat=4):
        Rlow[(a, s, m, n)] = trunc(sum(g[a][r] * Rup[(r, s, m, n)] for r in idx), 2)
    return Rlow


def laplacian(f, ginv, sqrtg, isqrtg):
    tot = 0
    for m in idx:
        inner = trunc(sqrtg * sum(ginv[m][n] * sp.diff(f, X[n]) for n in idx))
        tot += sp.diff(inner, X[m])
    return trunc(isqrtg * tot)


# ----------------------------------------------------------------- one background: checks + recursion
def run_background(R, label, do_checks=True):
    rep = {}
    Ric, Rs, Ric2, Riem2 = invariants(R)
    g, ginv, sqrtg, isqrtg, Dhalf, iDhalf = geometry(R)
    if do_checks:
        G = christoffel(g, ginv)
        # 1. normal-coordinate condition
        nc = [trunc(sum(G[m][a][b] * X[a] * X[b] for a in idx for b in idx), 3) for m in idx]
        rep['normal_coords'] = all(v == 0 for v in nc)
        # 2. curvature of the ansatz at the origin and its first derivative
        Rf = riemann_field(G, g)
        rep['riemann_at_0'] = all(sp.expand(at0(Rf[k]) - R[k]) == 0 for k in Rf)
        # 2b. Box R at the origin from the ansatz (R(x) through x^2); [a_2] contains (1/30) Box R, so it is
        #     fitted, not assumed away: the ansatz with nabla R = 0 does not make Box R vanish at order R.R.
        Rsf = trunc(sum(ginv[a][c] * ginv[b][d] * Rf[(a, b, c, d)] for a in idx for b in idx for c in idx for d in idx), 2)
        rep['R_at_0_matches'] = sp.expand(at0(Rsf) - Rs) == 0
        rep['dRiemann_at_0'] = all(at0(sp.diff(Rf[k], X[c])) == 0 for k in Rf for c in idx)
        rep['boxR_at_0'] = sp.expand(sum(sp.diff(Rsf, X[m], X[m]) for m in idx).subs({x: 0 for x in X}))
        R2part = homog(Rsf, 2)
        rep['R_x2_nonzero_coeffs'] = 0 if R2part == 0 else len(sp.Poly(R2part, *X).terms())
        # 3. Delta^{1/2}
        rep['Dhalf_x2'] = sp.expand(homog(Dhalf, 2) - sp.Rational(1, 12) * sum(Ric[(a, b)] * X[a] * X[b]
                                                                                for a in idx for b in idx)) == 0
    # recursion
    rhs1 = trunc(iDhalf * laplacian(Dhalf, ginv, sqrtg, isqrtg) + E, 2)
    a1 = sum(homog(rhs1, d) / (d + 1) for d in range(3))
    a1_0 = at0(a1)
    rep['a1_at_0'] = sp.expand(a1_0)
    rep['a1_check'] = sp.expand(a1_0 - (E + Rs / 6)) == 0
    f = trunc(Dhalf * a1, 2)
    a2_0 = sp.expand(sp.Rational(1, 2) * (at0(laplacian(f, ginv, sqrtg, isqrtg)) + E * a1_0))
    rep['a2_at_0'] = a2_0
    rep['inv'] = (Rs, Ric2, Riem2)
    rep['label'] = label
    return rep


def fit(reps):
    """[a_2] = cEE E^2 + cRE R E + c1 R^2 + c2 Ric^2 + c3 Riem^2; exact solve over the backgrounds given."""
    cEE, cRE, c1, c2, c3, c4 = sp.symbols('cEE cRE c1 c2 c3 c4')
    have_box = any(r['boxR_at_0'] != 0 for r in reps)
    eqs = []
    for r in reps:
        Rs, Ric2, Riem2 = r['inv']
        resid = sp.expand(r['a2_at_0'] - (cEE * E**2 + cRE * Rs * E + c1 * Rs**2 + c2 * Ric2 + c3 * Riem2
                                          + c4 * r['boxR_at_0']))
        P = sp.Poly(resid, E, *[s for s in resid.free_symbols if s not in (cEE, cRE, c1, c2, c3, c4, E)])
        eqs += [c for c in P.coeffs()]
    if not have_box:
        eqs = [e.subs(c4, 0) for e in eqs]
    sol = sp.solve(eqs, [cEE, cRE, c1, c2, c3] + ([c4] if have_box else []), dict=True)
    return sol, len(eqs), have_box


# ----------------------------------------------------------------- main
name = 'rung1 dewitt recursion ' + ('symbolic' if args.symbolic else f'samples={args.samples}') + (f' quartic={args.quartic}' if args.quartic != '2/45' else '')
say(f'# {name}')
say('# predictions: [a1] = E + R/6; [a2] = E^2/2 + R E/6 + R^2/72 - Ric^2/180 + Riem^2/180 + Box R/30; S4: 29/15; S2xS2: 11/45')
say(f'# x^4 term of the ansatz metric: coefficient {args.quartic}')

reps = []
if args.symbolic:
    R, syms = riemann_symbolic()
    say(f'# generic Riemann: {len(syms)} independent symbols')
    t1 = time.time()
    rep = run_background(R, 'symbolic')
    say(f'# symbolic background done in {time.time()-t1:.1f}s')
    reps.append(rep)
else:
    rng = random.Random(args.seed)
    Rlist = [riemann_sample(rng) for _ in range(args.samples)]
    if args.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(args.jobs) as pool:
            reps = pool.starmap(run_background, [(Rk, f'sample{i}') for i, Rk in enumerate(Rlist)])
    else:
        reps = [run_background(Rk, f'sample{i}') for i, Rk in enumerate(Rlist)]

status = 'OK'
for r in reps:
    checks = {k: v for k, v in r.items() if isinstance(v, bool)}
    say(f"{r['label']}: checks {checks}")
    say(f"   [a1] = {r['a1_at_0']}   Box R(0) = {r.get('boxR_at_0')}   nonzero x^2 coefficients of R(x): {r.get('R_x2_nonzero_coeffs')}")
    if not r['label'].startswith('symbolic'):
        say(f"   R = {r['inv'][0]}, Ric^2 = {r['inv'][1]}, Riem^2 = {r['inv'][2]},  [a2] = {r['a2_at_0']}")
    if not all(checks.values()):
        status = 'INPUT-CHECK-FAILED'

sol, neq, have_box = fit(reps)
say(f'# fit over {len(reps)} background(s), {neq} equations: {sol}')
if not have_box:
    say('# Box R(0) = 0 on every background with this ansatz: the Box R coefficient is not determined by this run (see --quartic 0)')
if len(sol) != 1 or any(k not in sol[0] for k in sp.symbols('cEE cRE c1 c2 c3') + ((sp.Symbol('c4'),) if have_box else ())):
    status = 'FIT-FAILED'
    say('# fit underdetermined or failed')
else:
    s = sol[0]
    c1, c2, c3 = s[sp.Symbol('c1')], s[sp.Symbol('c2')], s[sp.Symbol('c3')]
    con1 = sp.nsimplify(c2 + 4 * c3); con2 = sp.nsimplify(c1 - c3)
    say(f'# constraints from the six O(h^2) logs: c_Ric2 + 4 c_Riem2 = {con1} (want 1/60);  c_R2 - c_Riem2 = {con2} (want 1/120)')
    pred = {'cEE': sp.Rational(1, 2), 'cRE': sp.Rational(1, 6), 'c1': sp.Rational(1, 72),
            'c2': -sp.Rational(1, 180), 'c3': sp.Rational(1, 180), 'c4': sp.Rational(1, 30)}
    hit = all(sp.simplify(s[sp.Symbol(k)] - v) == 0 for k, v in pred.items() if sp.Symbol(k) in s)
    say(f'# equals the predicted (Vassilevich 4.28) values: {hit}')
    if con1 != sp.Rational(1, 60) or con2 != sp.Rational(1, 120):
        status = 'CONSTRAINT-VIOLATED'
    elif not hit:
        status = 'CONSTRAINTS-OK_NEW-NUMBERS'
    # Riem^2 = W^2 + 2 Ric^2 - R^2/3 ; Euler = Riem^2 - 4 Ric^2 + R^2 ; solve for the (W^2, R^2, Euler) split
    w, r2, eu = sp.symbols('w r2 eu')
    split = sp.solve([sp.Eq(w + eu, c3), sp.Eq(-2 * w - 4 * eu, c2), sp.Eq(w / 3 + r2 + eu, c1)], [w, r2, eu])
    say(f'# (W^2, R^2, Euler) split at E = 0: {split}   (W^2 = Riem^2 - 2Ric^2 + R^2/3, Euler = Riem^2 - 4Ric^2 + R^2)')

# named backgrounds, recursion only (their input checks are the same code path)
for nm, want in (('S4', sp.Rational(29, 15)), ('S2xS2', sp.Rational(11, 45))):
    r = run_background(riemann_named(nm), nm, do_checks=True)
    v = r['a2_at_0'].subs(E, 0)
    checks = {k: val for k, val in r.items() if isinstance(val, bool)}
    say(f"{nm}: checks {checks}  Box R(0) = {r['boxR_at_0']}  [a2]_(E=0) = {v}  (predicted {want}) {'ok' if v == want else 'MISMATCH'}")
    if v != want and status == 'OK':
        status = 'NAMED-BACKGROUND-MISMATCH'

say(f'# END {name} elapsed={time.time()-t0:.1f}s status={status}')
