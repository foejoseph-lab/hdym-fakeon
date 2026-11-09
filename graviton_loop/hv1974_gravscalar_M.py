#!/usr/bin/env python3
"""
Rung 3, step 3a (method M): gravity + a massless scalar, 't Hooft-Veltman 1974 sections 4-6 with phi-tilde != 0 -- the 11-dimensional
bundle (h_{ab}, phi) of (4.30) with X (4.31), Y (4.32), Z (4.33) as printed (re-read at 4x on 9 Nov), the P^{-1} replacement (5.3),
the ghost (5.19) with its phi-tilde term, through (3.40), pointwise on random algebraic data: a curvature tensor, a gradient
v_mu = D_mu phi-tilde and a symmetric W_{mu nu} = D_mu D_nu phi-tilde (independent off shell).  Targets: the off-shell total (5.22)
and the on-shell (6.5), (203/80) R^2.  Plan: ../notes/2026-11-04_rung3_thooft_veltman_summary.md section 2, step 3a.

Conventions (hv1974_graviton_M.py, established 8 Nov): standard curvature sign; the paper's Ric and R are MINUS the standard ones,
its R^{b m}_{a n} the standard Riemann; the phi-tilde terms carry no curvature and are taken as printed.  E_std = P^{-1} Xbar etc. with
that map (E_std = M_HV, since -R_HV/6 = +R_std/6 in (3.40)).  On shell (6.3) in standard sign: Ric = (1/2) v v, R = (1/2) v^2,
tr W = Box phi-tilde = 0.

The operator on Phi = (h, phi) from (4.30) = sqrt g( -1/2 Dh P Dh - 1/2 (D phi)^2 + 1/2 h X h + phi Y.h + 1/2 Z phi^2 ):
    E (h, phi) = ( P^{-1}(Xbar h + phi Y),  <Y, h> + Z phi ),   <Y, h> = sum_ab Y_ab h_ab,
    X^{b m}_{a n} = 2( -1/2 d^b_n v^m v_a + 1/4 d^b_a v^m v_n - 1/16 d^b_a d^m_n v^2 + 1/8 d^b_n d^m_a v^2
                       - [ -1/8 d^b_a d^m_n R + 1/4 d^b_n d^m_a R - 1/2 d^b_n Ric^m_a + 1/2 d^b_a Ric^m_n ] + 1/2 R^{b m}_{a n} )
    Y_ab = 1/2 d_ab tr W - W_ab,    Z = -v^2,
    ghost: E = +Ric - v v  (their (5.19): Box eta - R^{HV} eta - (d phi d phi) eta), complex, FP sign -2.
    Omega: the h-block's induced curvature, 0 on phi;  (1/120) x 11 x (Ric^2 - R^2/3).
Fit basis (dimension 4, even in phi-tilde): R^2, Ric^2, Riem^2, R v^2, Ric_mn v^m v^n, (v^2)^2, (tr W)^2, W.W.  Reductions
(total derivatives, as the paper's "apart from total derivatives"): Riem^2 -> 4 Ric^2 - R^2 (Gauss-Bonnet);
W.W -> (tr W)^2 - Ric_mn v^m v^n (D_m(W^mn v_n) - D_n(tr W v^n) = W.W - (tr W)^2 + Ric v v; the Ricci identity).

PREDICTIONS (rule 25; 9 Nov 17:30, before the first run; [v] = the paper under the convention map):
  P1 [v, (5.22) in standard sign] after the reductions: R^2 9/720 = 1/80, Ric^2 43/120, (v^2)^2 1/2, R v^2 +1/12 (their -1/12 with
      R_HV = -R_std), (tr W)^2 2, Ric v v 0.
      RE-SCORED after the first run (9 Nov 17:45): the run gave (tr W)^2 = 1/2 and everything else as predicted.  The paper's own
      (5.17), 2 Tr(Ybar P^-1 Ybar) = 2 (Box phi)^2 + 4 R^ab d phi d phi, carried through (3.35)'s factor 1/2 and the halving of the
      doubling, gives 1/2 (Box phi)^2 in the total, and its R^ab d phi d phi cancels against the ghost's -- which is why (5.22) has no
      such term, as ours has none.  The printed 2 (D_mu D^mu phi)^2 in (5.22) is inconsistent with (5.15)-(5.18), (5.21) and (3.35),
      which give 1/2 -- the slip is in the final addition, not in the pieces; 2 would need both the (3.35) factor 1/2 and the halving
      skipped for this one term.  It vanishes on shell and for pure gravity, so (5.24) and (6.5) never depended on it.  P1 = the five
      other coefficients + P1b, the (tr W)^2 coefficient scored against 1/2.
      FRESH REVIEW (9 Nov, separate agent, from the images): assembled (5.22) by hand from (5.14)-(5.18), (5.21) -> 1/2, every other
      coefficient as here; checked by sympy on ten fresh backgrounds that this script's X and Y reproduce the printed intermediates
      (5.15), (5.16), (5.17) exactly; mixed-block normalization, Ricci-identity sign, phi-tilde terms of (4.31), on-shell signs all
      confirmed (Ric = -vv/2 would give 529/960, not 203/320, so (6.5) pins the sign).  Caveat kept: the 1/2 rests on Y (4.32) as
      printed -- only (5.17) tests Y's normalization, (5.24) and (6.5) are blind to it -- so deriving (4.32) from the action (step 3b)
      is what takes the slip claim beyond "consistent with the paper's own pieces".
  P2 [v, (6.5)] on shell: substituting Ric = v v/2, R = v^2/2, tr W = 0 into P1 gives (1/320 + 43/480 + 1/2 + 1/24)(v^2)^2 =
      203/320 (v^2)^2 = (203/80) R^2 -- checked by hand before the run: the two printed results are consistent under the map.
  P3 [v, (5.21)] the ghost alone after the reductions, in standard sign: -(17/60) R^2 - (7/30) Ric^2 + (1/6) R v^2 + Ric v v
      - 1/2 (v^2)^2.  (The first version of this line had Ric v v with a minus sign -- my sign conversion, wrong: their -{... +
      R^ab d phi d phi ...} with Ric_HV = -Ric_std is +Ric_std v v.  Re-read at 4x and re-scored 9 Nov 17:45; the run had +1.)
  P4 [v] setting v = W = 0 reproduces hv1974_graviton_M.py's pure gravity plus the scalar loop (5.23): total (1/120 + 1/240) R^2 +
      (7/20 + 1/120) Ric^2 = (1/80) R^2 + (43/120) Ric^2 -- the same as P1's curvature part, as it must be.
  P5 [reported]: the pre-reduction coefficients, which the paper does not print.
Status PASS if P1-P4 hold.  Single method until the recursion runs this bundle (method P, next).

Run:  python3 hv1974_gravscalar_M.py --samples 12 --log results/hv1974_gravscalar_M.log        [seconds]
"""
import sys, time, argparse, itertools, random, platform
import sympy as sp

N = 4
idx = range(N)
PAIRS = [(a, b) for a in idx for b in idx if a <= b]
out = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def version_line():
    return f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}; {sys.executable}"


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


def scalar_sample(rng):
    v = sp.Matrix([sp.Rational(rng.randint(-3, 3), rng.randint(1, 3)) for _ in idx])
    W = sp.zeros(4, 4)
    for a in idx:
        for b in idx:
            if a <= b:
                W[a, b] = W[b, a] = sp.Rational(rng.randint(-3, 3), rng.randint(1, 3))
    return v, W


def invariants(R, v, W):
    Ric = sp.Matrix(4, 4, lambda b, d: sum(R[(a, b, a, d)] for a in idx))
    Rs = sp.expand(Ric.trace())
    inv = {'R2': Rs ** 2, 'Ric2': sp.expand((Ric * Ric).trace()), 'Riem2': sp.expand(sum(x ** 2 for x in R.values())),
           'Rv2': sp.expand(Rs * (v.T * v)[0]), 'Ricvv': sp.expand((v.T * Ric * v)[0]), 'v4': sp.expand((v.T * v)[0] ** 2),
           'trW2': sp.expand(W.trace() ** 2), 'WW': sp.expand((W * W).trace())}
    return Ric, Rs, inv


def basis_elem(p):
    a, b = p
    h = sp.zeros(4, 4); h[a, b] += 1
    if a != b:
        h[b, a] += 1
    return h


def coords(h):
    return [h[a, b] for (a, b) in PAIRS]


def op_matrix11(f):
    """11 x 11 matrix of a linear map on (h, phi): columns = images of the 10 basis matrices and of phi = 1."""
    cols = []
    for p in PAIRS:
        hh, ph = f(basis_elem(p), sp.Integer(0))
        cols.append(coords(hh) + [ph])
    hh, ph = f(sp.zeros(4, 4), sp.Integer(1))
    cols.append(coords(hh) + [ph])
    return sp.Matrix(11, 11, lambda i, j: cols[j][i])


def op_matrix10(f):
    cols = [coords(f(basis_elem(p))) for p in PAIRS]
    return sp.Matrix(10, 10, lambda i, j: cols[j][i])


def build_E(R, Ric, Rs, v, W, with_scalar=True):
    d = lambda i, j: sp.Integer(i == j)
    v2 = (v.T * v)[0]
    trW = W.trace()
    Gm = sp.diag(*[1 if a == b else 2 for (a, b) in PAIRS])

    def X_of(h):
        Xh = sp.zeros(4, 4)
        for a in idx:
            for b in idx:
                s = 0
                for m in idx:
                    for n in idx:
                        X = 2 * (-sp.Rational(1, 2) * d(b, n) * v[m] * v[a] + sp.Rational(1, 4) * d(b, a) * v[m] * v[n]
                                 - sp.Rational(1, 16) * d(b, a) * d(m, n) * v2 + sp.Rational(1, 8) * d(b, n) * d(m, a) * v2
                                 - (-sp.Rational(1, 8) * d(b, a) * d(m, n) * Rs + sp.Rational(1, 4) * d(b, n) * d(m, a) * Rs
                                    - sp.Rational(1, 2) * d(b, n) * Ric[m, a] + sp.Rational(1, 2) * d(b, a) * Ric[m, n])
                                 + sp.Rational(1, 2) * R[(b, m, a, n)])
                        s += X * h[n, m]
                Xh[a, b] = s
        return (Xh + Xh.T) / 2
    Xbar = op_matrix10(X_of)
    Xbar = (Xbar + Gm.inv() * Xbar.T * Gm) / 2            # the self-adjoint part, as a 10x10 in the basis
    Pinv = op_matrix10(lambda h: 2 * h - sp.eye(4) * h.trace())

    def mat_from_coords(c):
        h = sp.zeros(4, 4)
        for (a, b), val in zip(PAIRS, c):
            h[a, b] += val
            if a != b:
                h[b, a] += val
        return h
    Y = sp.Rational(1, 2) * sp.eye(4) * trW - W
    Z = -v2

    def E_of(h, ph):
        xh = mat_from_coords(list(Xbar * sp.Matrix(coords(h))))
        if not with_scalar:
            ph = sp.Integer(0)
        top = mat_from_coords(list(Pinv * sp.Matrix(coords(xh + ph * Y))))
        bot = sum(Y[a, b] * h[a, b] for a in idx for b in idx) + Z * ph
        if not with_scalar:
            bot = sp.Integer(0)
        return top, bot
    E = op_matrix11(E_of)
    # Omega^2 on the h-block
    def Om(m, n):
        def f(h):
            return sp.Matrix(4, 4, lambda a, b: -sum(R[(l, a, m, n)] * h[l, b] + R[(l, b, m, n)] * h[a, l] for l in idx))
        return f
    Om2 = sp.expand(sum((op_matrix10(Om(m, n)) * op_matrix10(Om(m, n))).trace() for m in idx for n in idx))
    return E, Om2


def hv340(dim, Om2, E, Rs, Ric2):
    X = E + Rs / 6 * sp.eye(dim)
    return sp.expand(sp.Rational(1, 24) * Om2 + sp.Rational(1, 4) * (X * X).trace() + sp.Rational(1, 120) * dim * (Ric2 - Rs ** 2 / 3))


KEYS = ['R2', 'Ric2', 'Riem2', 'Rv2', 'Ricvv', 'v4', 'trW2', 'WW']


def fit(rows):
    cs = sp.symbols('c0:8')
    eqs = [sp.Eq(sum(c * inv[k] for c, k in zip(cs, KEYS)), val) for val, inv in rows]
    sol = sp.solve(eqs[:8], cs, dict=True)
    if not sol:
        say("FIT-FAILED: degenerate backgrounds"); sys.exit(1)
    sol = sol[0]
    resid = all(sp.simplify(e.lhs.subs(sol) - e.rhs) == 0 for e in eqs[8:])
    return {k: sol[c] for c, k in zip(cs, KEYS)}, resid


def reduce(c):
    """Gauss-Bonnet and the Ricci identity: Riem^2 -> 4 Ric^2 - R^2; W.W -> (tr W)^2 - Ric v v."""
    r = dict(c)
    r['R2'] = sp.expand(c['R2'] - c['Riem2']); r['Ric2'] = sp.expand(c['Ric2'] + 4 * c['Riem2']); r['Riem2'] = 0
    r['trW2'] = sp.expand(c['trW2'] + c['WW']); r['Ricvv'] = sp.expand(c['Ricvv'] - c['WW']); r['WW'] = 0
    return r


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=12)
    ap.add_argument('--seed', type=int, default=20261109)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        import os
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; delete it to start over"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# hv1974_gravscalar_M: the (h, phi) bundle (11) and the ghost through (3.40) on {args.samples} random (curvature, v, W) backgrounds, seed {args.seed}")
    say(version_line())
    rng = random.Random(args.seed)
    rows = {'total': [], 'ghost': [], 'hphi': [], 'total_v0': []}
    for i in range(args.samples):
        R = riemann_sample(rng); v, W = scalar_sample(rng)
        Ric, Rs, inv = invariants(R, v, W)
        E, Om2 = build_E(R, Ric, Rs, v, W)
        hphi = hv340(11, Om2, E, Rs, inv['Ric2'])
        Eg = Ric - v * v.T
        ghost = sp.expand(-2 * hv340(4, -inv['Riem2'], Eg, Rs, inv['Ric2']))
        rows['hphi'].append((hphi, inv)); rows['ghost'].append((ghost, inv)); rows['total'].append((sp.expand(hphi + ghost), inv))
        # P4: the same backgrounds with v = W = 0
        v0, W0 = sp.zeros(4, 1), sp.zeros(4, 4)
        Ric0, Rs0, inv0 = invariants(R, v0, W0)
        E0, Om20 = build_E(R, Ric0, Rs0, v0, W0)
        tot0 = sp.expand(hv340(11, Om20, E0, Rs0, inv0['Ric2']) - 2 * hv340(4, -inv0['Riem2'], Ric0, Rs0, inv0['Ric2']))
        rows['total_v0'].append((tot0, inv0))
        say(f"sample{i}: R = {Rs}, v^2 = {(v.T * v)[0]}, tr W = {W.trace()}  [{time.time() - t0:.0f} s]")
    res = {}
    for k in ('hphi', 'ghost', 'total'):
        c, ok_fit = fit(rows[k]); res[k] = (c, reduce(c), ok_fit)
        say(f"\n## {k}: pre-reduction {c}  (residuals zero: {ok_fit})")
        say(f"   reduced (GB, Ricci identity): { {kk: vv for kk, vv in res[k][1].items() if vv != 0} }")
    # P4 with the v = 0 rows: fit only the curvature structures
    c4 = sp.symbols('d0:3')
    eqs = [sp.Eq(c4[0] * inv['R2'] + c4[1] * inv['Ric2'] + c4[2] * inv['Riem2'], val) for val, inv in rows['total_v0']]
    s4 = sp.solve(eqs[:3], c4, dict=True)[0]
    p4 = (sp.expand(s4[c4[0]] - s4[c4[2]]) == sp.Rational(1, 80)) and (sp.expand(s4[c4[1]] + 4 * s4[c4[2]]) == sp.Rational(43, 120)) \
        and all(sp.simplify(e.lhs.subs(s4) - e.rhs) == 0 for e in eqs[3:])
    tgt1 = {'R2': sp.Rational(1, 80), 'Ric2': sp.Rational(43, 120), 'v4': sp.Rational(1, 2), 'Rv2': sp.Rational(1, 12), 'Ricvv': 0}
    red = res['total'][1]
    p1 = all(sp.expand(red[k] - tgt1[k]) == 0 for k in tgt1)
    say(f"\nP1 (5.22) in standard sign after the reductions, five coefficients {tgt1}: {'PASS' if p1 else 'FAIL'}")
    p1b = sp.expand(red['trW2'] - sp.Rational(1, 2)) == 0
    say(f"P1b (Box phi)^2 coefficient: {red['trW2']}; from the paper's own (5.17) through (3.35) and the halving: 1/2; printed in (5.22): 2  -> {'PASS' if p1b else 'FAIL'} against 1/2")
    p1 = p1 and p1b
    # on shell: Ric -> vv/2, R -> v^2/2, tr W -> 0 in the reduced total: coefficient of (v^2)^2
    on = sp.expand(red['R2'] / 4 + red['Ric2'] / 4 + red['Rv2'] / 2 + red['Ricvv'] / 2 + red['v4'])
    p2 = on == sp.Rational(203, 320)
    say(f"P2 (6.5) on shell: total = {on} (v^2)^2 = {sp.nsimplify(on * 4)} R^2   (203/80 R^2 predicted): {'PASS' if p2 else 'FAIL'}")
    tgt3 = {'R2': -sp.Rational(17, 60), 'Ric2': -sp.Rational(7, 30), 'Rv2': sp.Rational(1, 6), 'Ricvv': 1, 'v4': -sp.Rational(1, 2), 'trW2': 0}
    redg = res['ghost'][1]
    p3 = all(sp.expand(redg[k] - tgt3[k]) == 0 for k in tgt3)
    say(f"P3 (5.21) ghost in standard sign {tgt3}: {'PASS' if p3 else 'FAIL'}")
    say(f"P4 v = W = 0 reproduces pure gravity + the scalar loop, (1/80, 43/120): {'PASS' if p4 else 'FAIL'}")
    ok = p1 and p2 and p3 and p4 and all(res[k][2] for k in res)
    say(f"\nstatus={'PASS' if ok else 'FAIL'}")
    say(f"END hv1974_gravscalar_M {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
