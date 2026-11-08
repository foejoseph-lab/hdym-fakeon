#!/usr/bin/env python3
"""
Rung 3, step 2a (method M): the graviton and its ghost as two more bundles for 't Hooft-Veltman's (3.40) -- their operator as
written ((4.30)-(4.31) at phi-tilde = 0, the P^{-1} replacement (5.3)-(5.4), the ghost (5.19)) evaluated pointwise on random
algebraic curvature tensors, their section 5 reproduced piece by piece before the sum.  Step 2b (deriving (4.31) from the action)
is the independent part and comes after.  Plan: ../notes/2026-11-04_rung3_thooft_veltman_summary.md section 2; the two
step-1 codes: ../matter_loop/hv1974_master.py (this formula), ../matter_loop/rung3_dewitt_bundle.py (method P, next for these bundles).

Conventions.  Standard curvature sign throughout (hv1974_master.py header: the paper's R is minus the standard one, so their
(M - R/6)^2 is the standard (E + R/6)^2, and since X is linear in curvature at phi-tilde = 0, E_std = -P^{-1} Xbar evaluated on
standard-sign tensors).  Euclidean frame at the point, indices up = down.  The h-bundle is the 10-dimensional space of symmetric
tensors h_{ab} in the basis E_aa, E_ab + E_ba (a < b); operators are 10 x 10 matrices in that basis (traces are basis-independent).
    X^{beta mu}_{alpha nu} = 2( -1/8 d^b_n d^m_a R + 1/4 d^b_a d^m_n R - 1/2 d^b_n R^m_a + 1/2 d^m_a R^b_n + 1/2 R^{b m}_{a n} )   (4.31)
    (Xbar h)_{ab} = sym_{ab} sum_{mn} X^{b m}_{a n} h_{n m};   P^{-1} h = 2 h - delta tr h   (5.4);   E_h = -P^{-1} Xbar
    (Omega_{mn} h)_{ab} = -R_{l a m n} h_{l b} - R_{l b m n} h_{a l}  (the induced connection's curvature; sign irrelevant in Tr Omega^2)
    ghost: the vector bundle with E = +Ric (their (5.19): Box eta^mu - R^{HV}_{a mu} eta^a = Box + Ric_std), complex (x2), FP sign (-).
    (3.40) real fields: dL = (sqrt g/eps) Tr{ (1/24) Omega^2 + 1/4 (E + R/6)^2 + (1/120)(Ric^2 - R^2/3) 1 }.

PREDICTIONS (rule 25; written 8 Nov 15:00 before the first run; all [v] from the paper's section 5 unless marked):
  P1  Tr_10 Omega^2 = -6 Riem^2  [(5.14): 6R^2 - 24 Ric^2 after GB, i.e. -6 Riem^2]   (also = (d + 2) x the vector's -Riem^2)
  P2  Tr (P^{-1} Xbar)^2 = 2 R^2 + 6 Ric^2  [(5.15) at phi-tilde = 0]
  P3  Tr (P^{-1} Xbar) = 6 R  [from (5.16): Tr(-1/3 P^{-1} Xbar R + R^2/36) = -31/18 R^2 over the 10 x 10 identity]
  P4  ghost (vector, E = +Ric, x(-2)) after GB = -(17/60) R^2 - (7/30) Ric^2  [(5.21) at phi-tilde = 0]
  P5  h-loop alone after GB = (7/24) R^2 + (7/12) Ric^2  [derived 8 Nov from (5.24) minus (5.21): the paper does not print it]
  P6  h-loop + ghost after GB = (1/120) R^2 + (7/20) Ric^2  [(5.24), pure gravity]
  P7  [the parametrisation check, obvious]: without the P^{-1} replacement (E = -Xbar) the numbers are wrong -- reported, not scored.
  Pre-GB coefficients (ours; the paper never fixes Riem^2): reported.
Status PASS if P1-P6 hold.  Single method until rung3_dewitt_bundle.py runs these bundles (method P); P2-P3 are the paper's
own intermediate traces, which is as close to a second method as a transcription can get.

RECORD (rule 11).  First run, 8 Nov 15:10, with E = -P^-1 Xbar from my FIRST TRANSCRIPTION of (4.31): P1 PASS, P4 PASS (Omega and
the ghost), P2, P3, P5, P6 FAIL -- Tr(P^-1 Xbar) = -10 R, no sign/index-order reading of the five terms AS TRANSCRIBED reproduces
(5.15)-(5.16) even modulo Gauss-Bonnet (results/hv1974_graviton_M_v1_transcribed431.log; --transcribed --misread reruns it).  The
operator was then DERIVED from the action (graviton_operator_derive.py; validated on Maxwell first):
    E h = -R h + (1/2) R delta tr h + (Ric h + h Ric) - delta tr(Ric h) - Ric tr h + 2 R_{a m b n} h_{m n},
Tr E = -6 R (their +6 R_HV), Tr E^2 = 2 R^2 + 6 Ric^2 + 3 GB ((5.15) modulo Gauss-Bonnet), (5.21), and (5.24).  The fresh review of
8 Nov derived E by hand independently (same six coefficients) and reconstructed which (4.31) would match it; re-reading the print
at 4x showed MY TRANSCRIPTION HAD TWO MISPLACED SUBSCRIPTS (the delta placements of the two R terms, and d^b_a R^m_n in the second
Ricci term).  Correctly transcribed, (4.31) equals the derived operator exactly on two seeds with E = P^-1 Xbar, the paper's Ric and
R read as minus the standard ones and their R^{b m}_{a n} as the standard Riemann.  The paper is consistent from operator to result;
the claim made against it between 15:10 and 18:30 on 8 Nov is withdrawn here and in the notes.  P2 re-scored "modulo GB", P3 as
Tr E = -6 R.  The default run uses the derived E; --transcribed uses the printed operator and must agree ("E == derived E: True").

Run:  python3 hv1974_graviton_M.py --samples 6 --log results/hv1974_graviton_M.log        [seconds]
"""
import sys, time, argparse, itertools, random
import sympy as sp

N = 4
idx = range(N)
PAIRS = [(a, b) for a in idx for b in idx if a <= b]          # 10 basis elements of symmetric tensors
out = None
MISREAD = False


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def riemann_sample(rng, nmat=8):
    """rung1_dewitt.py's generator: random rational algebraic curvature tensor, standard sign."""
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
    Rs = sp.expand(sum(Ric[(a, a)] for a in idx))
    Ric2 = sp.expand(sum(Ric[(a, b)] ** 2 for a in idx for b in idx))
    Riem2 = sp.expand(sum(v ** 2 for v in R.values()))
    return Ric, Rs, Ric2, Riem2


def basis_elem(p):
    a, b = p
    h = sp.zeros(4, 4)
    h[a, b] += 1
    if a != b:
        h[b, a] += 1
    return h


def coords(h):
    """coordinates of a symmetric matrix in the basis E_aa, E_ab + E_ba."""
    return [h[a, b] for (a, b) in PAIRS]


def op_matrix(f):
    """10 x 10 matrix of a linear map f on symmetric 4 x 4 matrices, columns = images of the basis elements."""
    cols = [coords(f(basis_elem(p))) for p in PAIRS]
    return sp.Matrix(10, 10, lambda i, j: cols[j][i])


def graviton_ops(R, Ric, Rs):
    d = lambda i, j: sp.Integer(i == j)

    def X_of(h, misread=False):
        """(4.31) at phi-tilde = 0 as PRINTED (re-read at 4x on 8 Nov 18:30):
             X^{b m}_{a n} = 2( -1/8 d^b_a d^m_n R + 1/4 d^b_n d^m_a R - 1/2 d^b_n R^m_a + 1/2 d^b_a R^m_n + 1/2 R^{b m}_{a n} ),
           with the paper's Ric and R equal to MINUS the standard ones and R^{b m}_{a n} the standard Riemann (the sector signs that
           make it equal the operator derived from the action, graviton_operator_derive.py; both seeds).  misread=True is the first
           transcription of 8 Nov 15:00 (two subscripts misplaced), kept as the record."""
        Xh = sp.zeros(4, 4)
        for a in idx:
            for b in idx:
                s = 0
                for m in idx:
                    for n in idx:
                        if misread:
                            X = 2 * (-sp.Rational(1, 8) * d(b, n) * d(m, a) * Rs + sp.Rational(1, 4) * d(b, a) * d(m, n) * Rs
                                     - sp.Rational(1, 2) * d(b, n) * Ric[(m, a)] + sp.Rational(1, 2) * d(m, a) * Ric[(b, n)]
                                     + sp.Rational(1, 2) * R[(b, m, a, n)])
                        else:
                            X = 2 * (-(-sp.Rational(1, 8) * d(b, a) * d(m, n) * Rs + sp.Rational(1, 4) * d(b, n) * d(m, a) * Rs
                                       - sp.Rational(1, 2) * d(b, n) * Ric[(m, a)] + sp.Rational(1, 2) * d(b, a) * Ric[(m, n)])
                                     + sp.Rational(1, 2) * R[(b, m, a, n)])
                        s += X * h[n, m]
                Xh[a, b] = s
        return (Xh + Xh.T) / 2

    def Pinv(h):
        return 2 * h - sp.eye(4) * h.trace()

    def Omega(m, n):
        def f(h):
            o = sp.zeros(4, 4)
            for a in idx:
                for b in idx:
                    o[a, b] = -sum(R[(l, a, m, n)] * h[l, b] + R[(l, b, m, n)] * h[a, l] for l in idx)
            return o
        return f

    Xbar = op_matrix(lambda h: X_of(h, misread=MISREAD))
    Gm = sp.diag(*[1 if a == b else 2 for (a, b) in PAIRS])        # Gram matrix of the basis for the inner product sum h_ab k_ab
    Xbar = (Xbar + Gm.inv() * Xbar.T * Gm) / 2                        # the quadratic form sees only the self-adjoint part (review, 8 Nov)
    PinvM = op_matrix(Pinv)
    PX = PinvM * Xbar
    Om = {(m, n): op_matrix(Omega(m, n)) for m in idx for n in idx}
    RicM = sp.Matrix(4, 4, lambda a, b: Ric[(a, b)])
    # the operator DERIVED from the action (graviton_operator_derive.py, 8 Nov), in closed form:
    E_derived = op_matrix(lambda h: -Rs * h + Rs / 2 * sp.eye(4) * h.trace() + (RicM * h + h * RicM) - sp.eye(4) * (RicM * h).trace()
                          - RicM * h.trace() + 2 * sp.Matrix(4, 4, lambda a, b: sum(R[(a, m, b, n)] * h[m, n] for m in idx for n in idx)))
    return Xbar, PX, Om, E_derived


def hv340(dim, Om2, E, Rs, Ric2):
    X = E + Rs / 6 * sp.eye(dim)
    return sp.expand(sp.Rational(1, 24) * Om2 + sp.Rational(1, 4) * (X * X).trace() + sp.Rational(1, 120) * dim * (Ric2 - Rs ** 2 / 3))


def vector_ghost(R, Ric, Rs, Ric2, Riem2):
    """(3.40) for the vector bundle with E = +Ric, times -2 (complex FP ghost)."""
    E = sp.Matrix(4, 4, lambda a, b: Ric[(a, b)])
    Om2 = -Riem2                                               # Tr Omega^2 on vectors (hv1974_master.py, verified 4 Nov)
    return sp.expand(-2 * hv340(4, Om2, E, Rs, Ric2))


def fit(rows):
    c1, c2, c3 = sp.symbols('c_R2 c_Ric2 c_Riem2')
    eqs = [sp.Eq(c1 * Rs ** 2 + c2 * Ric2 + c3 * Riem2, v) for v, Rs, Ric2, Riem2 in rows]
    sol = sp.solve(eqs[:3], [c1, c2, c3], dict=True)
    if not sol:
        say("FIT-FAILED: degenerate backgrounds"); sys.exit(1)
    sol = sol[0]
    resid = all(sp.simplify(e.lhs.subs(sol) - e.rhs) == 0 for e in eqs[3:])
    return sol[c1], sol[c2], sol[c3], resid


def gb(c):
    return sp.expand(c[0] - c[2]), sp.expand(c[1] + 4 * c[2])


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=6)
    ap.add_argument('--seed', type=int, default=20261108)
    ap.add_argument('--log', default=None)
    ap.add_argument('--transcribed', action='store_true', help='use P^-1 Xbar from the printed (4.31) instead of the derived E')
    ap.add_argument('--misread', action='store_true', help='with --transcribed: the first, wrong transcription of 8 Nov 15:00 (the record)')
    args = ap.parse_args()
    global MISREAD
    MISREAD = args.misread
    t0 = time.time()
    if args.log:
        import os
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; move it first"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# hv1974_graviton_M: the h-bundle (10) and the ghost (4) through (3.40) on {args.samples} random algebraic curvature tensors, seed {args.seed}; "
        f"E = {('P^-1 Xbar from the printed (4.31)' + (' [first, wrong transcription]' if args.misread else '')) if args.transcribed else 'the operator derived from the action (graviton_operator_derive.py)'}")
    rng = random.Random(args.seed)
    rows = {k: [] for k in ('h', 'ghost', 'total', 'h_noPinv')}
    p123 = True
    for i in range(args.samples):
        R = riemann_sample(rng)
        Ric, Rs, Ric2, Riem2 = invariants(R)
        Xbar, PX, Om, E_derived = graviton_ops(R, Ric, Rs)
        Om2 = sp.expand(sum((Om[(m, n)] * Om[(m, n)]).trace() for m in idx for n in idx))
        trPX2 = sp.expand((PX * PX).trace()); trPX = sp.expand(PX.trace())
        Eh_ = E_derived if not args.transcribed else (PX if not MISREAD else -PX)
        c1 = sp.expand(Om2 + 6 * Riem2) == 0
        GBv = sp.expand(Riem2 - 4 * Ric2 + Rs ** 2)
        trE2, trE = sp.expand((Eh_ * Eh_).trace()), sp.expand(Eh_.trace())
        c2 = sp.cancel((trE2 - (2 * Rs ** 2 + 6 * Ric2)) / GBv) == 3
        c3 = (trE == -6 * Rs)
        eq_derived = (Eh_ - E_derived) == sp.zeros(10, 10)
        p123 &= c1 and c2 and c3
        say(f"sample{i}: R = {Rs}, Ric^2 = {Ric2}, Riem^2 = {Riem2};  Tr Omega^2 + 6 Riem^2 = {sp.expand(Om2 + 6 * Riem2)} [P1 {c1}];  "
            f"[Tr E^2 - (2R^2 + 6Ric^2)]/GB = {sp.cancel((trE2 - 2 * Rs**2 - 6 * Ric2) / GBv)} [P2 {c2}];  Tr E / R = {sp.cancel(trE / Rs)} [P3 {c3}];  E == derived E: {eq_derived}")
        Eh = E_derived if not args.transcribed else (PX if not MISREAD else -PX)
        h_val = hv340(10, Om2, Eh, Rs, Ric2)
        g_val = vector_ghost(R, Ric, Rs, Ric2, Riem2)
        rows['h'].append((h_val, Rs, Ric2, Riem2)); rows['ghost'].append((g_val, Rs, Ric2, Riem2))
        rows['total'].append((sp.expand(h_val + g_val), Rs, Ric2, Riem2))
        rows['h_noPinv'].append((hv340(10, Om2, -Xbar, Rs, Ric2), Rs, Ric2, Riem2))
    res = {}
    for k in rows:
        c = fit(rows[k]); res[k] = c
        say(f"\n## {k}: pre-GB R^2 {c[0]}  Ric^2 {c[1]}  Riem^2 {c[2]}  (residuals zero: {c[3]});  after GB: R^2 {gb(c)[0]}  Ric^2 {gb(c)[1]}")
    p4 = gb(res['ghost']) == (-sp.Rational(17, 60), -sp.Rational(7, 30))
    p5 = gb(res['h']) == (sp.Rational(7, 24), sp.Rational(7, 12))
    p6 = gb(res['total']) == (sp.Rational(1, 120), sp.Rational(7, 20))
    say(f"\nP1-P3 (5.14), (5.15) mod GB, (5.16) on every background: {'PASS' if p123 else 'FAIL'}")
    say(f"P4 ghost = -(17/60) R^2 - (7/30) Ric^2 (5.21): {'PASS' if p4 else 'FAIL'}")
    say(f"P5 h-loop alone = (7/24) R^2 + (7/12) Ric^2 (derived): {'PASS' if p5 else 'FAIL'}")
    say(f"P6 pure gravity = (1/120) R^2 + (7/20) Ric^2 (5.24): {'PASS' if p6 else 'FAIL'}")
    say(f"P7 (reported) without the P^-1 replacement the h-loop after GB would be R^2 {gb(res['h_noPinv'])[0]}, Ric^2 {gb(res['h_noPinv'])[1]}")
    ok = p123 and p4 and p5 and p6
    say(f"\nstatus={'PASS' if ok else 'FAIL'}")
    say(f"END hv1974_graviton_M {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
