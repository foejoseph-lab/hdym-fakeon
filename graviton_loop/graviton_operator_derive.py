#!/usr/bin/env python3
"""
Rung 3, step 2b: the gauge-fixed graviton fluctuation operator DERIVED from the action, not transcribed -- the second variation
of 't Hooft-Veltman's L = sqrt(g) R (4.1, up to the overall sign, which drops out) plus the gauge term -1/2 C^2 with C of (4.27) at
phi-tilde = 0, on a generic background in Riemann normal coordinates, as the Euler-Lagrange operator at the origin:
    O h = s ( P D^2 h + Xbar h ),   P = 1/2 g g - 1/4 g g (5.2),   E = P^{-1} Xbar  (Gilkey's endomorphism of -(D^2 + E)).
E is independent of the overall sign and normalisation s of the action; the RELATIVE sign of the gauge term is fixed by minimality
(no mixed-index second-derivative terms), which is measured, not assumed.
Why: step 2a (hv1974_graviton_M.py) found no reading of the printed (4.31) that reproduces both (5.15) and (5.16).

Method (one point, exact rationals).  Background g(x) = delta - (1/3) R x x + (2/45) R.R x x x x on a random algebraic curvature
tensor (standard sign).  For a polynomial base field b(x) the Euler-Lagrange expression at the origin is obtained with a probe
b + dl (H0 + H1 x + 1/2 H2 x x) (150 jet symbols, first order in dl): with G0, G1, G2 the coefficient functions of the symbols in
L2(x), the Lagrangian's jet-derivatives along the base field are Phi0 = G0, Phi1_m = -x^m G0 + G1_m, Phi2_mn = 1/2 x^m x^n G0
- x^(n G1_m) + G2_mn (the polynomial coefficients c and the jet j at x are related by c0 = j0 - j1 x + 1/2 j2 xx, c1 = j1 - j2 x,
c2 = j2), and EL(0) = [Phi0 - d_m Phi1_m + d_m d_n Phi2_mn](0) with plain x-derivatives -- checked on L = (h')^2/2, b = x^2/2 -> -1.
Base fields: e_k (10 constant symmetric matrices) give O e_k = s(P D^2 e_k + Xbar e_k) with (D^2 e_k)(0) the connection term,
computed from the ansatz's Gamma; x_m^2/2 e_q and x_m x_n e_q give the second-derivative coefficients (s and minimality).

PREDICTIONS (rule 25; 8 Nov 16:30, before the first run) -- kept as written, the record follows:
  P1 [obvious]: with the right relative sign of the gauge term, EL[x_m^2/2 e_q](0) = s P e_q for every m and EL[x_m x_n e_q](0) = 0
      (minimal operator); the connection term is -(1/3)(Ric h + h Ric) [hand, normal coordinates].
  P2 [(5.15)-(5.16), the paper's own traces, via the derived E]: Tr E = 6 R and Tr E^2 = 2 R^2 + 6 Ric^2, up to the sign of E.
  P3 [(5.24)]: (3.40) with the derived E plus the ghost gives (1/120) R^2 + (7/20) Ric^2 after Gauss-Bonnet; h-loop alone
      (7/24) R^2 + (7/12) Ric^2.
  P4 [reported]: the derived Xbar against the transcribed (4.31), term by term.

RECORD OF THE FIRST RUNS (rule 11; sandbox, 8 Nov 17:00-17:40, seed 20261108):
  - The machinery was validated on Maxwell in Feynman gauge before being trusted for the graviton (/tmp/maxwell_test.py, kept as
    maxwell_operator_check.py): minimal, s = -1, E = -Ric exactly, connection term +(1/3) Ric on constant vectors.
  - P1's connection-term sign was WRONG by hand: it is +(1/3)(Ric h + h Ric) (d_mu Gamma^c_{mu a}(0) = -(1/3) Ric^c_a); the code
    computes it from Gamma, so nothing downstream used the hand sign.  Re-scored.
  - A bug of mine, caught by the Maxwell check: for a symmetric field the Euler-Lagrange derivative with respect to the independent
    off-diagonal variable h_{ab} (a < b) is 2 (O h)_{ab}, and the first version stored it unhalved -- E then fitted no covariant
    structure and s read 2 from an off-diagonal probe.  Fixed (val/2); s = 1.
  - With the fix the derived E is EXACTLY a combination of the six curvature structures that map symmetric tensors to symmetric
    tensors:  E h = -R h + (1/2) R delta tr h + (Ric h + h Ric) - delta tr(Ric h) - Ric tr h + 2 R_{a m b n} h_{m n}.
    Tr E = -6 R (their +6 R with R_HV = -R_std: (5.16) MATCH); Tr E^2 = 2 R^2 + 6 Ric^2 + 3 GB, GB = Riem^2 - 4 Ric^2 + R^2:
    (5.15) holds MODULO Gauss-Bonnet (P2 re-scored to "mod GB"; the paper's "apart from total derivatives" after (5.14) applies).
    (3.40) with this E plus the ghost = (1/120) R^2 + (7/20) Ric^2 + (7/12) GB: (5.24) EXACTLY after Gauss-Bonnet; the h-loop
    alone = (7/24) R^2 + (7/12) Ric^2 + (1/2) GB.  P3 holds (the multi-background fit is hv1974_graviton_M.py with this E).
  - P4, as first written 17:40: "the printed (4.31) ... is not the operator of section 5".  WITHDRAWN 18:30: the fresh review derived
    E by hand (same six coefficients) and reconstructed which (4.31) would match; re-reading the print at 4x showed my transcription
    had two misplaced subscripts.  Correctly transcribed, (4.31) IS the derived operator (hv1974_graviton_M.py header).  What
    remains ours: the derivation as an independent confirmation, and the convention statement (their Ric, R minus standard; their
    Riemann standard) which the matter-loop sign argument could not fix for the Riemann sector.
  - Review of 8 Nov, folded in: P3's comparison with the `expected` closed form is circular (the values came from this code's
    earlier run) -- the independent confirmation is the reviewer's hand derivation and the agreement with the printed (4.31); the
    minimality probe was one traceless element, so a trace-sector probe (x_0 x_1 e_(0,0)) is added; the Maxwell check is a script
    for the author to run.
Status PASS if P1 (minimal) holds and the structure fit of E is exact.

Run:  python3 graviton_operator_derive.py --log results/graviton_operator_derive.log      [sandbox timing in the header record]
"""
import sys, time, argparse, itertools, random
import sympy as sp

N = 4
idx = range(N)
X = sp.symbols('x1:5', real=True)
eps, dl = sp.symbols('eps dl')
PAIRS = [(a, b) for a in idx for b in idx if a <= b]
out = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


# ---------------------------------------------------------------- sparse polynomial ring (8 Nov: sympy expand() was too slow)
from sympy import ring, QQ
JET0 = [f'H0_{a}{b}' for (a, b) in PAIRS]
JET1 = [f'H1_{a}{b}_{m}' for (a, b) in PAIRS for m in idx]
JET2 = [f'H2_{a}{b}_{m}{n}' for (a, b) in PAIRS for m in idx for n in idx if m <= n]
RING, *GENS = ring('x1,x2,x3,x4,eps,dl', QQ)
XR = GENS[:4]; EPS, DL = GENS[4], GENS[5]
ONE = RING(1)


def q(r):
    r = sp.Rational(r)
    return QQ(int(r.p), int(r.q))


def trunc(p, xdeg=2, edeg=2, ddeg=1):
    """keep monomials with x-degree <= xdeg, eps-degree <= edeg, dl-degree <= ddeg (ring elements)."""
    if not isinstance(p, type(ONE)):
        p = RING(p)
    return RING({m: c for m, c in p.terms() if sum(m[:4]) <= xdeg and m[4] <= edeg and m[5] <= ddeg})


def at0(p):
    """value at x = 0 (as a ring element in the remaining generators)."""
    return RING({m: c for m, c in p.terms() if sum(m[:4]) == 0})


def coeff_dl(p):
    """coefficient of dl^1 in p, as a ring element without dl."""
    return RING({m[:5] + (0,): c for m, c in p.terms() if m[5] == 1})


def to_sympy(p):
    return p.as_expr()


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
    Rs = sp.expand(sum(Ric[(a, a)] for a in idx))
    Ric2 = sp.expand(sum(Ric[(a, b)] ** 2 for a in idx for b in idx))
    Riem2 = sp.expand(sum(v ** 2 for v in R.values()))
    return Ric, Rs, Ric2, Riem2


def background(R):
    """4x4 list of ring elements: g = delta - (1/3) R x x + (2/45) R.R x x x x."""
    g = [[None] * 4 for _ in idx]
    for mu in idx:
        for nu in idx:
            q2 = sum((RING(q(-sp.Rational(1, 3) * R[(mu, a, nu, b)])) * XR[a] * XR[b] for a in idx for b in idx), RING(0))
            q4 = sum((RING(q(sp.Rational(2, 45) * R[(a, mu, b, l)] * R[(l, c, nu, d)])) * XR[a] * XR[b] * XR[c] * XR[d]
                      for a in idx for b in idx for c in idx for d in idx for l in idx), RING(0))
            g[mu][nu] = RING(1 if mu == nu else 0) + q2 + q4
    return g


def mat_mul(A, B, xdeg):
    return [[trunc(sum((A[i][k] * B[k][j] for k in idx), RING(0)), xdeg) for j in idx] for i in idx]


def mat_add(A, B, s=1):
    return [[A[i][j] + s * B[i][j] for j in idx] for i in idx]


def eye():
    return [[RING(1 if i == j else 0) for j in idx] for i in idx]


def series_inverse(g, xdeg):
    A = mat_add(g, eye(), -1)
    inv = eye(); term = eye()
    for _ in range(4):
        term = mat_mul([[-t for t in row] for row in term], A, xdeg)
        inv = mat_add(inv, term)
    return [[trunc(e, xdeg) for e in row] for row in inv]


def mtrace(A):
    return sum((A[i][i] for i in idx), RING(0))


def basis_elem(p):
    a, b = p
    h = sp.zeros(4, 4); h[a, b] += 1
    if a != b:
        h[b, a] += 1
    return h


def coords(h):
    return [h[a, b] for (a, b) in PAIRS]


def op_matrix(f):
    cols = [coords(f(basis_elem(p))) for p in PAIRS]
    return sp.Matrix(10, 10, lambda i, j: cols[j][i])


def sym_to_ring_matrix(M):
    return [[RING(q(M[i, j])) if not isinstance(M[i, j], type(ONE)) else M[i, j] for j in idx] for i in idx]


class Geometry:
    def __init__(self, R):
        self.R = R
        self.gb = background(R)
        self.gbinv = series_inverse(self.gb, 2)
        dgb = [[[trunc(self.gb[a][b].diff(XR[c]), 3) for c in idx] for b in idx] for a in idx]
        self.Gb = [[[trunc(sum((self.gbinv[m][l] * (dgb[l][b][a] + dgb[l][a][b] - dgb[a][b][l]) for l in idx), RING(0)) / 2, 3)
                     for b in idx] for a in idx] for m in idx]
        A = mat_add(self.gb, eye(), -1); A2 = mat_mul(A, A, 2)
        logdet = trunc(mtrace(A) - mtrace(A2) / 2, 2)
        self.sqrtgb = trunc(1 + logdet / 2 + logdet ** 2 / 8, 2)

    def lagrangian_L2(self, h, gauge_sign):
        """h: 4x4 ring elements (contains eps^0 only; multiplied by eps here).  Returns the eps^2 coefficient (ring, no eps)."""
        T = trunc
        G = [[self.gb[a][b] + EPS * h[a][b] for b in idx] for a in idx]
        Ginv = series_inverse(G, 2)
        A = mat_add(G, eye(), -1); A2 = mat_mul(A, A, 2); A3 = mat_mul(A2, A, 2)
        logdet = T(mtrace(A) - mtrace(A2) / 2 + mtrace(A3) / 3, 2)
        sqrtg = T(1 + logdet / 2 + logdet ** 2 / 8 + logdet ** 3 / 48, 2)
        dG = [[[T(G[a][b].diff(XR[c]), 3) for c in idx] for b in idx] for a in idx]
        Gam3 = [[[T(sum((Ginv[m][l] * (dG[l][b][a] + dG[l][a][b] - dG[a][b][l]) for l in idx), RING(0)) / 2, 3)
                  for b in idx] for a in idx] for m in idx]
        Gam1 = [[[T(Gam3[m][a][b], 1) for b in idx] for a in idx] for m in idx]
        Ric = [[None] * 4 for _ in idx]
        for b in idx:
            for c in idx:
                Ric[b][c] = T(sum((Gam3[a][b][c].diff(XR[a]) - Gam3[a][b][a].diff(XR[c]) for a in idx), RING(0))
                              + sum((Gam1[a][a][d] * Gam1[d][b][c] - Gam1[a][c][d] * Gam1[d][b][a] for a in idx for d in idx), RING(0)), 2)
        Rsc = T(sum((Ginv[b][c] * Ric[b][c] for b in idx for c in idx), RING(0)), 2)
        L = T(sqrtg * Rsc, 2)
        Gb, gbinv = self.Gb, self.gbinv
        Dh = [[[T(h[a][b].diff(XR[c]) - sum((Gb[l][c][a] * h[l][b] + Gb[l][c][b] * h[a][l] for l in idx), RING(0)), 3)
                for c in idx] for b in idx] for a in idx]
        Dtrh = [T(sum((gbinv[a][b] * Dh[a][b][c] for a in idx for b in idx), RING(0)), 3) for c in idx]
        Cl = [T(sum((gbinv[n][r] * Dh[m][n][r] for n in idx for r in idx), RING(0)) - Dtrh[m] / 2, 3) for m in idx]
        C2 = T(self.sqrtgb * sum((gbinv[m][n] * Cl[m] * Cl[n] for m in idx for n in idx), RING(0)), 2)
        L = T(L + gauge_sign * QQ(-1, 2) * EPS ** 2 * C2, 2)
        return RING({m[:4] + (0,) + m[5:]: c for m, c in L.terms() if m[4] == 2})

    def connection_term(self):
        dGam0 = {(m, c, a): to_sympy(at0(self.Gb[c][m][a].diff(XR[m]))) for m in idx for c in idx for a in idx}

        def f(h):
            o = sp.zeros(4, 4)
            for a in idx:
                for b in idx:
                    o[a, b] = -sum(dGam0[(m, c, a)] * h[c, b] + dGam0[(m, c, b)] * h[a, c] for m in idx for c in idx)
            return o
        return f


def probe_fields():
    """the 150 jet directions as numeric monomial fields: (kind, key, 4x4 ring matrix)."""
    P = []
    for (a, b) in PAIRS:
        E = [[RING(1 if {i, j} == {a, b} else 0) for j in idx] for i in idx] if a != b else [[RING(1 if (i, j) == (a, a) else 0) for j in idx] for i in idx]
        P.append(('0', (a, b), E))
        for m in idx:
            P.append(('1', (a, b, m), [[E[i][j] * XR[m] for j in idx] for i in idx]))
        for m in idx:
            for n in idx:
                if m <= n:
                    w = QQ(1, 2) if m == n else QQ(1)
                    P.append(('2', (a, b, m, n), [[E[i][j] * XR[m] * XR[n] * w for j in idx] for i in idx]))
    return P


def el_at_origin(geo, base, gauge_sign, probes):
    """base: sympy 4x4 of polynomials in X; EL_{ab}(0) by 150 numeric probes."""
    base_r = [[ring_from_sympy(base[a, b]) for b in idx] for a in idx]
    L2b = geo.lagrangian_L2(base_r, gauge_sign)           # not needed for the dl-linear part, kept for a sanity check
    G = {}
    for kind, key, p in probes:
        h = [[base_r[a][b] + DL * p[a][b] for b in idx] for a in idx]
        G[(kind, key)] = coeff_dl(geo.lagrangian_L2(h, gauge_sign))
    EL = sp.zeros(4, 4)
    for (a, b) in PAIRS:
        G0 = G[('0', (a, b))]
        G1 = [G[('1', (a, b, m))] for m in idx]
        G2 = {(m, n): G[('2', (a, b, m, n))] for m in idx for n in idx if m <= n}
        expr = G0
        for m in idx:
            expr -= (-XR[m] * G0 + G1[m]).diff(XR[m])
        for m in idx:
            for n in idx:
                mm, nn = min(m, n), max(m, n)
                w = QQ(1) if m == n else QQ(1, 2)
                F2 = XR[m] * XR[n] * G0 / 2 - (XR[n] * G1[m] + XR[m] * G1[n]) / 2 + w * G2[(mm, nn)]
                expr += F2.diff(XR[m]).diff(XR[n])
        val = to_sympy(at0(expr))
        if a != b:
            val = val / 2          # the independent variable h_{ab} (a < b) enters the action as h_{ab} + h_{ba}: dS/dh_{ab} = 2 (O h)_{ab}
        EL[a, b] = val
        if a != b:
            EL[b, a] = val
    return EL


def ring_from_sympy(e):
    e = sp.expand(e)
    if e == 0:
        return RING(0)
    P = sp.Poly(e, *X)
    out_ = RING(0)
    for mon, c in P.terms():
        t = RING(q(c))
        for xi, m in zip(XR, mon):
            t = t * xi ** m
        out_ += t
    return out_


def selfcheck_1d():
    """L = (h')^2/2, base h = x^2/2: EL = -h'' = -1, through the same jet formula."""
    x = sp.Symbol('x'); c0, c1, c2 = sp.symbols('c0 c1 c2')
    h = x ** 2 / 2 + dl * (c0 + c1 * x + c2 * x ** 2 / 2)
    L = sp.expand(sp.diff(h, x) ** 2 / 2).coeff(dl, 1)
    G0, G1, G2 = L.coeff(c0, 1), L.coeff(c1, 1), L.coeff(c2, 1)
    EL = G0 - sp.diff(-x * G0 + G1, x) + sp.diff(x ** 2 / 2 * G0 - x * G1 + G2, x, 2)
    return sp.simplify(EL.subs(x, 0))


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261108)
    ap.add_argument('--gauge_sign', type=int, default=1)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        import os
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; move it first"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# graviton_operator_derive: seed {args.seed}, gauge_sign {args.gauge_sign}; 1-d self-check of the EL formula: {selfcheck_1d()} (want -1)")
    rng = random.Random(args.seed)
    R = riemann_sample(rng)
    Ric, Rs, Ric2, Riem2 = invariants(R)
    say(f"# background: R = {Rs}, Ric^2 = {Ric2}, Riem^2 = {Riem2}")
    geo = Geometry(R)
    jets = probe_fields()
    RicM = sp.Matrix(4, 4, lambda a, b: Ric[(a, b)])
    Pop = op_matrix(lambda h: h / 2 - sp.eye(4) * h.trace() / 4)
    Pinv = op_matrix(lambda h: 2 * h - sp.eye(4) * h.trace())
    # ---- P1: second-derivative coefficients and normalisation
    t1 = time.time()
    q = PAIRS[1]                                   # e_(0,1), an off-diagonal element, as the probe field
    EL_diag = [el_at_origin(geo, X[m] ** 2 / 2 * basis_elem(q), args.gauge_sign, jets) for m in (0, 2)]
    EL_mixed = el_at_origin(geo, X[0] * X[2] * basis_elem(q), args.gauge_sign, jets)
    EL_trace = el_at_origin(geo, X[0] * X[1] * basis_elem(PAIRS[0]), args.gauge_sign, jets)      # trace-sector probe (review, 8 Nov)
    Pq = (basis_elem(q) / 2 - sp.eye(4) * basis_elem(q).trace() / 4)
    s_list = []
    for ELm in EL_diag:
        ratio = {sp.cancel(ELm[a, b] / Pq[a, b]) for a in idx for b in idx if Pq[a, b] != 0}
        zero_ok = all(ELm[a, b] == 0 for a in idx for b in idx if Pq[a, b] == 0)
        s_list.append((ratio, zero_ok))
    say(f"# second-derivative probes ({time.time() - t1:.0f} s): EL[x_m^2/2 e_q] / (P e_q) = {s_list};  EL[x_0 x_2 e_q] = {coords(EL_mixed)};  EL[x_0 x_1 e_00] = {coords(EL_trace)}")
    p1 = (all(len(r) == 1 and z for r, z in s_list) and len({tuple(r)[0] for r, z in s_list}) == 1 and all(v == 0 for v in coords(EL_mixed))
          and all(v == 0 for v in coords(EL_trace)))
    s = tuple(s_list[0][0])[0] if p1 else None
    say(f"P1 minimal with s = {s}: {'PASS' if p1 else 'FAIL'}")
    if not p1:
        say("  (not minimal: rerun with --gauge_sign -1)"); say(f"END graviton_operator_derive {time.time() - t0:.0f} s"); return
    # ---- O on constant fields
    cols = []
    for k, p in enumerate(PAIRS):
        t1 = time.time()
        EL = el_at_origin(geo, basis_elem(p), args.gauge_sign, jets)
        cols.append(coords(EL))
        say(f"# column {k} e_{p} done in {time.time() - t1:.0f} s")
    O = sp.Matrix(10, 10, lambda i, j: cols[j][i])
    D2c = op_matrix(geo.connection_term())
    conn_ok = D2c == op_matrix(lambda h: (RicM * h + h * RicM) / 3)
    say(f"# connection term = +(1/3)(Ric h + h Ric) [the header's hand sign was wrong; computed from Gamma]: {conn_ok}")
    Xbar = O / s - Pop * D2c
    E = Pinv * Xbar
    say(f"# E = P^-1 Xbar (10x10, basis E_aa, E_ab+E_ba):\n{E.tolist()}")
    trE, trE2 = sp.expand(E.trace()), sp.expand((E * E).trace())
    GB = sp.expand(Riem2 - 4 * Ric2 + Rs ** 2)
    say(f"# Tr E = {sp.cancel(trE / Rs)} R;   Tr E^2 - (2R^2 + 6Ric^2) = {sp.cancel((trE2 - 2 * Rs**2 - 6 * Ric2) / GB)} * GB")
    p2 = (sp.cancel(trE / Rs) == -6) and (sp.cancel((trE2 - 2 * Rs ** 2 - 6 * Ric2) / GB) == 3)
    say(f"P2 (5.16) exactly and (5.15) modulo Gauss-Bonnet [re-scored]: {'PASS' if p2 else 'FAIL'}")
    # structure fit
    structs = {'R*1': op_matrix(lambda h: Rs * h), 'R*dd': op_matrix(lambda h: Rs * sp.eye(4) * h.trace()),
               'Ric.h+h.Ric': op_matrix(lambda h: RicM * h + h * RicM), 'd(x)Ric': op_matrix(lambda h: sp.eye(4) * (RicM * h).trace()),
               'Ric(x)d': op_matrix(lambda h: RicM * h.trace()),
               'Riem': op_matrix(lambda h: sp.Matrix(4, 4, lambda a, b: sum(R[(a, m, b, n)] * h[m, n] for m in idx for n in idx)))}
    names = list(structs); cs = sp.symbols('c0:6')
    M = sum((c * structs[n] for c, n in zip(cs, names)), sp.zeros(10, 10))
    sol = sp.solve([sp.Eq(M[i, j], E[i, j]) for i in range(10) for j in range(10)], cs, dict=True)
    fitE = {n: sol[0][c] for c, n in zip(cs, names)} if sol else None
    expected = {'R*1': -1, 'R*dd': sp.Rational(1, 2), 'Ric.h+h.Ric': 1, 'd(x)Ric': -1, 'Ric(x)d': -1, 'Riem': 2}
    say(f"# E as a combination of the six covariant structures: {fitE}   (the 8 Nov derived values: {expected})")
    p3 = fitE is not None and all(sp.simplify(fitE[k] - expected[k]) == 0 for k in expected)
    say(f"P3 structure fit exact and equal to the closed form of 8 Nov (circular against this code's own first run; the independent checks are the reviewer's hand derivation and the printed (4.31)): {'PASS' if p3 else 'FAIL'}")
    # (3.40) on this background with the derived E
    def Om_sym(m, n):
        def f(h):
            return sp.Matrix(4, 4, lambda a, b: -sum(R[(l, a, m, n)] * h[l, b] + R[(l, b, m, n)] * h[a, l] for l in idx))
        return f
    Om2 = sp.expand(sum((op_matrix(Om_sym(m, n)) * op_matrix(Om_sym(m, n))).trace() for m in idx for n in idx))
    Xm = E + Rs / 6 * sp.eye(10)
    h_val = sp.expand(Om2 / 24 + (Xm * Xm).trace() / 4 + sp.Rational(10, 120) * (Ric2 - Rs ** 2 / 3))
    Xg = RicM + Rs / 6 * sp.eye(4)
    ghost = sp.expand(-2 * (-Riem2 / 24 + (Xg * Xg).trace() / 4 + sp.Rational(4, 120) * (Ric2 - Rs ** 2 / 3)))
    tot = sp.expand(h_val + ghost)
    say(f"# (3.40): h-loop + ghost - [(1/120) R^2 + (7/20) Ric^2] = {sp.cancel((tot - Rs**2 / 120 - sp.Rational(7, 20) * Ric2) / GB)} * GB;  "
        f"h-loop - [(7/24) R^2 + (7/12) Ric^2] = {sp.cancel((h_val - sp.Rational(7, 24) * Rs**2 - sp.Rational(7, 12) * Ric2) / GB)} * GB   (one background: reported; the fit is hv1974_graviton_M.py)")
    # P4: the transcribed (4.31) in structure form
    def X431(h):
        d = lambda i, j: sp.Integer(i == j)
        Xh = sp.zeros(4, 4)
        for a in idx:
            for b in idx:
                Xh[a, b] = sum(2 * (-sp.Rational(1, 8) * d(b, n) * d(m, a) * Rs + sp.Rational(1, 4) * d(b, a) * d(m, n) * Rs
                                    - sp.Rational(1, 2) * d(b, n) * Ric[(m, a)] + sp.Rational(1, 2) * d(m, a) * Ric[(b, n)]
                                    + sp.Rational(1, 2) * R[(b, m, a, n)]) * h[n, m] for m in idx for n in idx)
        return (Xh + Xh.T) / 2
    X4 = op_matrix(X431)
    Gm = sp.diag(*[1 if a == b else 2 for (a, b) in PAIRS])
    X4s = (X4 + Gm.inv() * X4.T * Gm) / 2
    for sign, lab in ((1, '+'), (-1, '-')):
        E4 = sign * Pinv * X4s
        sol4 = sp.solve([sp.Eq(M[i, j], E4[i, j]) for i in range(10) for j in range(10)], cs, dict=True)
        say(f"# P4 (the record): {lab}P^-1 Xbar from the FIRST, WRONG transcription of (4.31) in structure form: {({n: sol4[0][c] for c, n in zip(cs, names)} if sol4 else None)}")
    def X431_print(h):
        d = lambda i, j: sp.Integer(i == j)
        Xh = sp.zeros(4, 4)
        for a in idx:
            for b in idx:
                Xh[a, b] = sum(2 * (-(-sp.Rational(1, 8) * d(b, a) * d(m, n) * Rs + sp.Rational(1, 4) * d(b, n) * d(m, a) * Rs
                                      - sp.Rational(1, 2) * d(b, n) * Ric[(m, a)] + sp.Rational(1, 2) * d(b, a) * Ric[(m, n)])
                                    + sp.Rational(1, 2) * R[(b, m, a, n)]) * h[n, m] for m in idx for n in idx)
        return (Xh + Xh.T) / 2
    X4p = op_matrix(X431_print); X4ps = (X4p + Gm.inv() * X4p.T * Gm) / 2
    say(f"# P4: (4.31) as PRINTED (re-read 8 Nov 18:30; their Ric, R = minus standard, their Riemann standard): P^-1 Xbar == derived E: {(Pinv * X4ps - E) == sp.zeros(10, 10)}")
    say(f"status={'PASS' if (p1 and p3) else 'FAIL'}  (P1 {p1}, P2 {p2}, P3 {p3})")
    say(f"END graviton_operator_derive {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
