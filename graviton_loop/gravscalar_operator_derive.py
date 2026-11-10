#!/usr/bin/env python3
"""
Rung 3, step 3b (part one): the gauge-fixed fluctuation operator of gravity + a massless scalar DERIVED from the action, not
transcribed -- the second variation of 't Hooft-Veltman's L = sqrt(g)(-R - 1/2 d phi-bar g d phi-bar) (4.1) [their R is minus the
standard one: in standard sign  sqrt(g)(R - 1/2 (d phi-bar)^2)], g-bar = g + h, phi-bar = phi-tilde + phi, plus the gauge term
-1/2 C^2 with C_mu = D^nu h_{mu nu} - 1/2 D_mu h - c phi d_mu phi-tilde (4.27; c = 1 in the paper), on a background carrying a
curvature tensor AND a scalar field phi-tilde(x), as the Euler-Lagrange operator on Phi = (h_ab, phi) at the origin of Riemann normal
coordinates:
    O Phi = s Pbig ( D^2 Phi + E Phi ),   Pbig = diag(P, 1),  P = 1/2 g g - 1/4 g g (5.2),   E = Gilkey's endomorphism of -(D^2 + E).
The method is graviton_operator_derive.py's (step 2b: numeric jet probes, exact rationals; validated there on Maxwell and on the
graviton), rewritten on python-flint (fmpq_mpoly) and extended to 11 components; the step-2 script is left untouched.  With
phi-tilde = 0 this script must reproduce step 2's derived E (run first, as the validation on a known case).
What it is for: hv1974_gravscalar_M.py (method M) and rung3_dewitt_bundle.py (method P) both take E from the PRINTED (4.31)-(4.33);
the fresh review of 9 Nov made that the open point of step 3a.  Here E comes from the action.

Backgrounds.  g(x) = delta - (1/3) R x x + (2/45) R.R x x x x on an algebraic curvature tensor; phi-tilde(x) = v.x + 1/2 W x x + cubic
(random).  The paper expands with the classical fields ON SHELL ((6.3); in standard sign Ric = 1/2 v v, R = 1/2 v^2, tr W = 0), so
its printed X, Y, Z may have been simplified with the field equations.  Two kinds of background are therefore run:
  on shell at the origin: random v, random traceless W, random Weyl part, Ricci part set to 1/2 v v;
  off shell: curvature, v, W all independent.
The comparison target is the sibling's build_E (hv1974_gravscalar_M.py: the printed (4.31)-(4.33) under the convention map of 8 Nov),
imported and evaluated at the same (Riemann, v, W) -- computed in place, nothing typed.

PREDICTIONS (rule 25; written 9 Nov 19:45, before the first run; nothing fetched):
  D0 [validation on a known case]: with phi-tilde = 0 the h-block of E is step 2's closed form
      -R h + 1/2 R d tr h + (Ric h + h Ric) - d tr(Ric h) - Ric tr h + 2 Riem h, the scalar decouples with E_phiphi = 0, and s = 1.
  D1 [derived by hand: the cross term h^{mn} v_m d_n phi - 1/2 h v.d phi of the matter Lagrangian is cancelled by the gauge term's
      c phi v^m (D^n h_mn - 1/2 D_m h) after one integration by parts iff c = 1]: with c = 1 the operator has NO first-derivative
      terms at the origin (O[x_m e_k](0) = 0 for all 44 linear base fields) and its second-derivative part is Pbig delta^{mn};
      with c = 0 the first-derivative terms are there (run with --c 0 as the control).
  D2 [the content]: ON SHELL the derived E equals the printed one as an 11 x 11 matrix.
  D3 [can fail; it decides the (5.22) claim of 9 Nov]: the mixed block is  (E phi)_ab = -2 W_ab,  (E h)_phi = 1/2 tr W tr h - W.h
      -- i.e. Y (4.32) with the print's normalisation -- and E_phiphi = -v^2.  If Y's normalisation comes out different, the
      (Box phi)^2 coefficient 1/2 claimed against the printed 2 of (5.22) falls with it.
  D4 [a question, not a prediction]: OFF SHELL, does the derived E equal the print?  If not, is the difference built from the field
      equations (EOM_ab = Ric_ab - 1/2 v_a v_b, and tr W)?  I do not know whether the print used them.
  D5 [by dimension]: the cubic jets of phi-tilde do not enter E(0).  (Shown by E equalling the print, which has no cubic jets, on
      backgrounds that carry them; --nocubic is there for a direct comparison of two logs.)
Status PASS if D0, D1, D2, D3 hold on at least one on-shell and one off-shell background; D4 and D5 are reported.

SMOKE RECORD (rule 11; sandbox 9 Nov 20:00-20:20, 2 cores, ~4-7 min per background; three partial runs, one cut by an environment
restart, none of them a complete status=PASS of this version): D0 holds; one on-shell and one off-shell background: s = 1, no
first-derivative terms, E == printed (11 x 11) on BOTH -- D4's answer is yes, the print is the off-shell operator, no field equations
were used; D3 holds.  Control --c 0: 22 nonzero first-derivative entries, (E phi)_ab != -2W, E_phiphi != -v^2, h-block and the row
unchanged.
FRESH REVIEW (9 Nov, separate agent; 5 runs, mutations on scratch copies; repository untouched):
  * No error found in the derivation; the comparison with build_E is a real second source for the (h, phi) block of E (derive() uses
    nothing from the sibling; the reviewer derived the phi-tilde terms of X, Y = 1/2 g tr W - W and Z = -c^2 v^2 by hand as a third,
    and predicted the c = 0 and c = -1 failure patterns before running them).
  * Mutations: --c -1 caught (D1, D3); scalar kinetic weight 1 instead of 1/2 caught through s (it CRASHED instead of printing FAIL:
    fixed below); phi-tilde's quadratic jet W instead of W/2 caught (D3, D4); metric inverse to x^2 instead of x^3: no change (step
    2's order was adequate).  Pool under a forced 'spawn' start method (Windows' way): passes.
  * SCOPE -- the ghost is NOT derived here.  Its E = Ric - v v (5.19) is still from the print in code (the reviewer's hand variation
    of the gauge condition gives D^2 eta + Ric eta - c (eta.v) v, agreeing with the print and with no W term).  So this script makes
    E two-source for the (h, phi) loop and settles the (Box phi)^2 coefficient (it depends only on the mixed block, both halves now
    derived, and on the absence of mixed first-derivative terms, D1); the ghost row (5.21), and with it the other coefficients of
    (5.22) and 203/80, wait for the ghost operator derived in code.
  * Folded in: s = None no longer crashes; PASS needs at least one on-shell and one off-shell background, and D2 requires the
    field equations to hold at the origin on the backgrounds it counts; exit code 1 on FAIL; v is drawn with no zero component.
  * After the fold-in (sandbox 9 Nov 21:05, --onshell 1 --offshell 1 --jobs 2, 779 s): one complete run, status=PASS, D0-D3 hold, D4 True.
  * Kept as limits: "no first-derivative terms" is measured at the origin (everywhere by covariance); the randomized probes miss a
    given nonzero coefficient pattern with probability of order 1/20 per combination, ~1/500 for the two (not "~1"); only 3 of the 44
    diagonal second-derivative pairs are probed by default -- one background with --full closes both; E itself comes from the 11
    constant probes only and is untouched by any of this.

Run:  python3 gravscalar_operator_derive.py --log results/gravscalar_operator_derive.log
"""
import sys, os, time, argparse, itertools, random, platform
import sympy as sp
import flint

N = 4
idx = range(N)
PAIRS = [(a, b) for a in idx for b in idx if a <= b]
out = None
CTX = flint.fmpq_mpoly_ctx.get(('x1', 'x2', 'x3', 'x4', 'eps', 'dl'), 'degrevlex')
GENS = CTX.gens()
XR = GENS[:4]; EPS, DL = GENS[4], GENS[5]
ZERO = CTX.from_dict({})
ONE = CTX.constant(flint.fmpq(1))


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def version_line():
    return f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}, python-flint {flint.__version__}; {sys.executable}"


def q(r):
    r = sp.Rational(r)
    return flint.fmpq(int(r.p), int(r.q))


def C(r):
    return CTX.constant(q(r))


def T(p, xdeg, edeg=2, ddeg=1):
    """keep monomials with x-degree <= xdeg, eps-degree <= edeg, dl-degree <= ddeg"""
    return CTX.from_dict({m: c for m, c in p.to_dict().items() if m[0] + m[1] + m[2] + m[3] <= xdeg and m[4] <= edeg and m[5] <= ddeg})


def d(p, m):
    return p.derivative(m)


def at0_rational(p):
    """the constant term of a polynomial in x only, as a sympy Rational"""
    dd = p.to_dict()
    c = dd.get((0, 0, 0, 0, 0, 0), None)
    return sp.Rational(int(c.p), int(c.q)) if c is not None else sp.Integer(0)


def psum(it):
    tot = ZERO
    for e in it:
        tot = tot + e
    return tot


def from_sympy(e):
    """a sympy polynomial in x1..x4 -> ring element"""
    e = sp.expand(e)
    if e == 0:
        return ZERO
    P = sp.Poly(e, *SX)
    return CTX.from_dict({mon + (0, 0): q(c) for mon, c in P.terms()})


SX = sp.symbols('x1:5', real=True)


def riemann_sample(rng, nmat=8):
    R = {k: sp.Integer(0) for k in itertools.product(idx, repeat=4)}
    for _ in range(nmat):
        A = [[None] * N for _ in idx]
        for a in idx:
            for b in idx:
                if a <= b:
                    A[a][b] = A[b][a] = sp.Rational(rng.randint(-3, 3), rng.randint(1, 3))
        w = sp.Rational(rng.randint(-2, 2) or 1, rng.randint(1, 2))
        for a, b, c, dd in itertools.product(idx, repeat=4):
            R[(a, b, c, dd)] += w * (A[a][c] * A[b][dd] - A[a][dd] * A[b][c])
    return R


def ricci_of(R):
    Ric = sp.Matrix(4, 4, lambda b, dd: sum(R[(a, b, a, dd)] for a in idx))
    return Ric, sp.expand(Ric.trace())


def with_ricci(R, Ric_target):
    """the same Weyl part, the Ricci part replaced: R - g^S(R) + g^S(target), S = 1/2 (Ric - R g/6), (g^S)_abcd = g_ac S_bd - g_ad S_bc + g_bd S_ac - g_bc S_ad"""
    Ric0, Rs0 = ricci_of(R)
    S0 = (Ric0 - Rs0 * sp.eye(4) / 6) / 2
    S1 = (Ric_target - Ric_target.trace() * sp.eye(4) / 6) / 2
    dS = S1 - S0
    dl_ = lambda i, j: sp.Integer(i == j)
    return {(a, b, c, e): sp.expand(R[(a, b, c, e)] + dl_(a, c) * dS[b, e] - dl_(a, e) * dS[b, c] + dl_(b, e) * dS[a, c] - dl_(b, c) * dS[a, e])
            for a, b, c, e in itertools.product(idx, repeat=4)}


def scalar_sample(rng, cubic=True):
    v = sp.Matrix([sp.Rational(rng.choice((-3, -2, -1, 1, 2, 3)), rng.randint(1, 3)) for _ in idx])          # no zero component (review, 9 Nov)
    W = sp.zeros(4, 4)
    for a in idx:
        for b in idx:
            if a <= b:
                W[a, b] = W[b, a] = sp.Rational(rng.randint(-3, 3), rng.randint(1, 3))
    cub = {m: (sp.Rational(rng.randint(-3, 3), rng.randint(1, 3)) if cubic else sp.Integer(0)) for m in itertools.product(range(4), repeat=4) if sum(m) == 3}
    return v, W, cub


def basis_elem(p):
    a, b = p
    h = sp.zeros(4, 4); h[a, b] += 1
    if a != b:
        h[b, a] += 1
    return h


def coords(h):
    return [h[a, b] for (a, b) in PAIRS]


def op_matrix10(f):
    cols = [coords(f(basis_elem(p))) for p in PAIRS]
    return sp.Matrix(10, 10, lambda i, j: cols[j][i])


def mat_mul(A, B, xdeg):
    return [[T(psum(A[i][k] * B[k][j] for k in idx), xdeg) for j in idx] for i in idx]


def eye():
    return [[ONE if i == j else ZERO for j in idx] for i in idx]


def series_inverse(g, xdeg):
    A = [[g[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]
    inv = eye(); term = eye()
    for _ in range(6):
        term = mat_mul([[-t for t in row] for row in term], A, xdeg)
        inv = [[inv[i][j] + term[i][j] for j in idx] for i in idx]
    return inv


class Background:
    def __init__(self, R, v, W, cub):
        self.R, self.v, self.W = R, v, W
        g = [[None] * 4 for _ in idx]
        for mu in idx:
            for nu in idx:
                q2 = psum(C(-sp.Rational(1, 3) * R[(mu, a, nu, b)]) * XR[a] * XR[b] for a in idx for b in idx)
                q4 = psum(C(sp.Rational(2, 45) * R[(a, mu, b, l)] * R[(l, c, nu, e)]) * XR[a] * XR[b] * XR[c] * XR[e]
                          for a in idx for b in idx for c in idx for e in idx for l in idx)
                g[mu][nu] = (ONE if mu == nu else ZERO) + q2 + q4
        self.gb = g
        self.gbinv = series_inverse(g, 3)
        dgb = [[[d(g[a][b], c) for c in idx] for b in idx] for a in idx]
        self.Gb = [[[T(psum(self.gbinv[m][l] * (dgb[l][b][a] + dgb[l][a][b] - dgb[a][b][l]) for l in idx) * C(sp.Rational(1, 2)), 3)
                     for b in idx] for a in idx] for m in idx]
        A = [[g[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]; A2 = mat_mul(A, A, 2)
        logdet = T(psum(A[i][i] for i in idx) - psum(A2[i][i] for i in idx) * C(sp.Rational(1, 2)), 2)
        self.sqrtgb = T(ONE + logdet * C(sp.Rational(1, 2)) + logdet * logdet * C(sp.Rational(1, 8)), 2)
        # phi-tilde = v.x + 1/2 W x x + cubic: in normal coordinates d_a phi(0) = v_a, D_a D_b phi(0) = W_ab (Gamma(0) = 0)
        phi = psum(C(v[a]) * XR[a] for a in idx) + psum(C(W[a, b] / 2) * XR[a] * XR[b] for a in idx for b in idx)
        phi = phi + psum(C(c) * XR[0]**m[0] * XR[1]**m[1] * XR[2]**m[2] * XR[3]**m[3] for m, c in cub.items())
        self.phit = phi
        self.dphit = [d(phi, a) for a in idx]

    def L2(self, h, ph, c_gauge, gauge_sign=1):
        """h: 4x4 ring elements, ph: ring element (the fluctuations, eps^0; multiplied by eps here).  The eps^2 coefficient of
        sqrt(G)(R[G] - 1/2 G^{mn} d Phi d Phi) + gauge_sign (-1/2) sqrt(g) g^{mn} C_m C_n, to x^2."""
        G = [[self.gb[a][b] + EPS * h[a][b] for b in idx] for a in idx]
        Ginv = series_inverse(G, 3)
        A = [[G[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]; A2 = mat_mul(A, A, 2); A3 = mat_mul(A2, A, 2)
        tr = lambda M: psum(M[i][i] for i in idx)
        logdet = T(tr(A) - tr(A2) * C(sp.Rational(1, 2)) + tr(A3) * C(sp.Rational(1, 3)), 2)
        sqrtg = T(ONE + logdet * C(sp.Rational(1, 2)) + logdet * logdet * C(sp.Rational(1, 8)) + logdet * logdet * logdet * C(sp.Rational(1, 48)), 2)
        dG = [[[T(d(G[a][b], c), 3) for c in idx] for b in idx] for a in idx]
        Gam = [[[T(psum(Ginv[m][l] * (dG[l][b][a] + dG[l][a][b] - dG[a][b][l]) for l in idx) * C(sp.Rational(1, 2)), 3)
                 for b in idx] for a in idx] for m in idx]
        Ric = [[None] * 4 for _ in idx]
        for b in idx:
            for c in idx:
                Ric[b][c] = T(psum(d(Gam[a][b][c], a) - d(Gam[a][b][a], c) for a in idx)
                              + psum(Gam[a][a][e] * Gam[e][b][c] - Gam[a][c][e] * Gam[e][b][a] for a in idx for e in idx), 2)
        Rsc = T(psum(Ginv[b][c] * Ric[b][c] for b in idx for c in idx), 2)
        dPhi = [self.dphit[a] + EPS * d(ph, a) for a in idx]
        kin = T(psum(Ginv[m][n] * dPhi[m] * dPhi[n] for m in idx for n in idx), 2)
        L = T(sqrtg * (Rsc - kin * C(sp.Rational(1, 2))), 2)
        Gb, gbinv = self.Gb, self.gbinv
        Dh = [[[T(d(h[a][b], c) - psum(Gb[l][c][a] * h[l][b] + Gb[l][c][b] * h[a][l] for l in idx), 3)
                for c in idx] for b in idx] for a in idx]
        Dtrh = [T(psum(gbinv[a][b] * Dh[a][b][c] for a in idx for b in idx), 3) for c in idx]
        Cl = [T(psum(gbinv[n][r] * Dh[m][n][r] for n in idx for r in idx) - Dtrh[m] * C(sp.Rational(1, 2)) - C(c_gauge) * ph * self.dphit[m], 3) for m in idx]
        C2 = T(self.sqrtgb * psum(gbinv[m][n] * Cl[m] * Cl[n] for m in idx for n in idx), 2)
        L = T(L + C(-sp.Rational(gauge_sign, 2)) * EPS * EPS * C2, 2)
        return CTX.from_dict({m[:4] + (0,) + m[5:]: c for m, c in L.to_dict().items() if m[4] == 2})

    def connection_term(self):
        """(D^2 h)(0) on a constant symmetric h: -(d_m Gamma^c_{m a})(0) h_cb - (a <-> b)"""
        dGam0 = {(m, c, a): at0_rational(T(d(self.Gb[c][m][a], m), 0)) for m in idx for c in idx for a in idx}

        def f(h):
            o = sp.zeros(4, 4)
            for a in idx:
                for b in idx:
                    o[a, b] = -sum(dGam0[(m, c, a)] * h[c, b] + dGam0[(m, c, b)] * h[a, c] for m in idx for c in idx)
            return o
        return f


# ---- the 11 field components: 10 symmetric pairs + the scalar ('s')
COMPS = PAIRS + ['s']


def comp_field(comp, mono):
    """the field (h 4x4, phi) equal to `mono` (a ring element) in component comp, zero elsewhere"""
    h = [[ZERO for _ in idx] for _ in idx]; ph = ZERO
    if comp == 's':
        ph = mono
    else:
        a, b = comp
        h[a][b] = mono
        if a != b:
            h[b][a] = mono
    return h, ph


def probe_fields():
    P = []
    for comp in COMPS:
        P.append(('0', comp, None, comp_field(comp, ONE)))
        for m in idx:
            P.append(('1', comp, m, comp_field(comp, XR[m])))
        for m in idx:
            for n in idx:
                if m <= n:
                    w = C(sp.Rational(1, 2)) if m == n else ONE
                    P.append(('2', comp, (m, n), comp_field(comp, XR[m] * XR[n] * w)))
    return P


def coeff_dl(p):
    return CTX.from_dict({m[:5] + (0,): c for m, c in p.to_dict().items() if m[5] == 1})


def el_at_origin(bg, base, c_gauge, probes, gauge_sign=1):
    """base = (h 4x4, phi) ring elements in x; returns the 11 Euler-Lagrange components at the origin (sympy Rationals, order COMPS)"""
    bh, bp = base
    G = {}
    for kind, comp, key, (ph_h, ph_p) in probes:
        h = [[bh[a][b] + DL * ph_h[a][b] for b in idx] for a in idx]
        G[(kind, comp, key)] = coeff_dl(bg.L2(h, bp + DL * ph_p, c_gauge, gauge_sign))
    half = C(sp.Rational(1, 2))
    outv = []
    for comp in COMPS:
        G0 = G[('0', comp, None)]
        G1 = [G[('1', comp, m)] for m in idx]
        expr = G0
        for m in idx:
            expr = expr - d(G1[m] - XR[m] * G0, m)
        for m in idx:
            for n in idx:
                mm, nn = min(m, n), max(m, n)
                w = ONE if m == n else half
                F2 = XR[m] * XR[n] * G0 * half - (XR[n] * G1[m] + XR[m] * G1[n]) * half + w * G[('2', comp, (mm, nn))]
                expr = expr + d(d(F2, m), n)
        val = at0_rational(T(expr, 0))
        if comp != 's' and comp[0] != comp[1]:
            val = val / 2          # the independent variable h_ab (a < b) enters as h_ab + h_ba: dS/dh_ab = 2 (O h)_ab
        outv.append(val)
    return outv


# ---- base fields as plain data (so a pool can carry them): a list of (component index, x-exponents (4), rational coefficient)
def base_from_spec(spec):
    h = [[ZERO for _ in idx] for _ in idx]; ph = ZERO
    for k, mon, c in spec:
        term = C(c) * XR[0]**mon[0] * XR[1]**mon[1] * XR[2]**mon[2] * XR[3]**mon[3]
        if COMPS[k] == 's':
            ph = ph + term
        else:
            a, b = COMPS[k]
            h[a][b] = h[a][b] + term
            if a != b:
                h[b][a] = h[b][a] + term
    return h, ph


_W = {}


def _init(Rd, v, W, cub, c_gauge):
    _W['bg'] = Background(Rd, v, W, cub); _W['probes'] = probe_fields(); _W['c'] = c_gauge


def _el(spec):
    return el_at_origin(_W['bg'], base_from_spec(spec), _W['c'], _W['probes'])


def mono(m, n=None):
    e = [0, 0, 0, 0]; e[m] += 1
    if n is not None:
        e[n] += 1
    return tuple(e)


def derive(Rd, v, W, cub, c_gauge, rng, jobs=1, full=False):
    """returns s, flags, and the 11x11 E = Pbig^-1 O/s - connection on constants.
    Second- and first-derivative parts: by default a randomized identity test (O is linear in the base field, so a nonzero
    coefficient pattern survives a random rational combination except with probability of order 1/20): three diagonal quadratic probes for s, two random
    combinations of all 66 mixed quadratic base fields (must give 0), two random combinations of all 44 linear base fields (must
    give 0).  --full probes every base field separately."""
    Pop = op_matrix10(lambda h: h / 2 - sp.eye(4) * h.trace() / 4)
    Pbig = sp.diag(Pop, sp.Matrix([[1]]))
    rr = lambda: sp.Rational(rng.randint(1, 5) * rng.choice((-1, 1)), rng.randint(1, 3))
    tasks = []          # (tag, spec)
    for k in range(11):
        tasks.append((('const', k), [(k, (0, 0, 0, 0), 1)]))
    if full:
        for k in range(11):
            for m in idx:
                tasks.append((('lin', k, m), [(k, mono(m), 1)]))
                tasks.append((('diag', k, m), [(k, mono(m, m), sp.Rational(1, 2))]))
            for m in idx:
                for n in idx:
                    if m < n:
                        tasks.append((('mixed', k, m, n), [(k, mono(m, n), 1)]))
    else:
        for (k, m) in ((1, 0), (0, 2), (10, 3)):                                   # e_(0,1), e_(0,0), the scalar
            tasks.append((('diag', k, m), [(k, mono(m, m), sp.Rational(1, 2))]))
        for j in range(2):
            tasks.append((('lin', 'random', j), [(k, mono(m), rr()) for k in range(11) for m in idx]))
            tasks.append((('mixed', 'random', j), [(k, mono(m, n), rr()) for k in range(11) for m in idx for n in idx if m < n]))
    specs = [t[1] for t in tasks]
    if jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs, initializer=_init, initargs=(Rd, v, W, cub, c_gauge)) as pool:
            vals = pool.map(_el, specs, chunksize=1)
    else:
        _init(Rd, v, W, cub, c_gauge)
        vals = [_el(sp_) for sp_ in specs]
    res = {t[0]: sp.Matrix(val) for t, val in zip(tasks, vals)}
    s_set = set(); second_ok = True; first_zero = True; n_first_nonzero = 0
    for tag, col in res.items():
        if tag[0] == 'diag':
            want = Pbig[:, tag[1]]
            ratios = {sp.cancel(col[i] / want[i]) for i in range(11) if want[i] != 0}
            second_ok &= len(ratios) == 1 and all(col[i] == 0 for i in range(11) if want[i] == 0)
            s_set |= ratios
        elif tag[0] == 'mixed':
            second_ok &= all(x == 0 for x in col)
        elif tag[0] == 'lin':
            nz = sum(1 for x in col if x != 0); n_first_nonzero += nz; first_zero &= nz == 0
    second_ok &= len(s_set) == 1
    s = list(s_set)[0] if len(s_set) == 1 else None
    O = sp.Matrix(11, 11, lambda i, j: res[('const', j)][i])
    conn = sp.diag(op_matrix10(_W['bg'].connection_term() if jobs == 1 else Background(Rd, v, W, cub).connection_term()), sp.Matrix([[0]]))
    E = (Pbig.inv() * O / s - conn) if s is not None else None
    return s, second_ok, first_zero, n_first_nonzero, E, sorted(s_set, key=str)


def step2_E(R):
    Ric, Rs = ricci_of(R)
    return op_matrix10(lambda h: -Rs * h + Rs / 2 * sp.eye(4) * h.trace() + Ric * h + h * Ric - sp.eye(4) * (Ric * h).trace() - Ric * h.trace()
                       + 2 * sp.Matrix(4, 4, lambda a, b: sum(R[(a, m, b, n)] * h[m, n] for m in idx for n in idx)))


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261109)
    ap.add_argument('--c', type=int, default=1, help='coefficient of phi d phi-tilde in the gauge condition (the paper: 1; 0 is the control for D1)')
    ap.add_argument('--onshell', type=int, default=2)
    ap.add_argument('--offshell', type=int, default=2)
    ap.add_argument('--nocubic', action='store_true', help='phi-tilde without cubic jets (D5: E must not change)')
    ap.add_argument('--full', action='store_true', help='probe every linear and quadratic base field separately instead of the randomized test (about 9 times the work)')
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--skip_d0', action='store_true', help='skip the phi-tilde = 0 validation (smokes only; a run without D0 never prints status=PASS)')
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# gravscalar_operator_derive: seed {args.seed}, c = {args.c}, {args.onshell} on-shell + {args.offshell} off-shell backgrounds, cubic jets {'off' if args.nocubic else 'on'}, derivative-part probes {'all' if args.full else 'randomized'}, jobs {args.jobs}")
    say(version_line())
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import hv1974_gravscalar_M as GM
    rng = random.Random(args.seed)
    rng_probe = random.Random(args.seed + 2)
    ok = True
    # ---- D0: phi-tilde = 0
    R = riemann_sample(rng)
    if not args.skip_d0:
        s, second_ok, first_zero, nfz, E, sset = derive(R, sp.zeros(4, 1), sp.zeros(4, 4), {}, args.c, rng_probe, args.jobs, args.full)
        if E is None:
            say(f"D0 FAIL: the diagonal second-derivative probes do not give one normalisation s: {sset}"); say("\nstatus=FAIL"); say(f"END gravscalar_operator_derive {time.time() - t0:.0f} s"); sys.exit(1)
        d0 = second_ok and first_zero and s == 1 and E[:10, :10] == step2_E(R) and E[10, :] == sp.zeros(1, 11) and E[:, 10] == sp.zeros(11, 1)
        say(f"D0 phi-tilde = 0: s = {s}, second-derivative part Pbig delta: {second_ok}, no first-derivative terms: {first_zero}, h-block == step 2's closed form: {E[:10, :10] == step2_E(R)}, scalar decoupled with E = 0: {E[10, :] == sp.zeros(1, 11) and E[:, 10] == sp.zeros(11, 1)}  -> {'PASS' if d0 else 'FAIL'}  [{time.time() - t0:.0f} s]")
        ok &= d0
    # ---- backgrounds with the scalar
    results = []
    for kind, n in (('on shell', args.onshell), ('off shell', args.offshell)):
        for i in range(n):
            R = riemann_sample(rng); v, W, cub = scalar_sample(rng, cubic=not args.nocubic)
            if kind == 'on shell':
                W = W - W.trace() * sp.eye(4) / 4
                R = with_ricci(R, v * v.T / 2)
            Ric, Rs = ricci_of(R)
            s, second_ok, first_zero, nfz, E, sset = derive(R, v, W, cub, args.c, rng_probe, args.jobs, args.full)
            if E is None:
                say(f"{kind} {i}: FAIL: the diagonal second-derivative probes (e_(0,1), e_(0,0), the scalar) do not give one normalisation s: {sset} -- the relative weight of the h and phi rows is off; no E"); say("\nstatus=FAIL"); say(f"END gravscalar_operator_derive {time.time() - t0:.0f} s"); sys.exit(1)
            E_print, _ = GM.build_E(R, Ric, Rs, v, W)
            Ed = E.applyfunc(sp.expand); Ep = E_print.applyfunc(sp.expand)
            eq_full = Ed == Ep
            Wc = sp.Matrix(coords(W))
            mixed_col = Ed[:10, 10] == -2 * Wc                                                          # (E phi)_ab = -2 W_ab
            brow_want = sp.Matrix([[sum((W.trace() / 2 * sp.eye(4) - W)[a, b] * basis_elem(p)[a, b] for a in idx for b in idx) for p in PAIRS]])
            mixed_row = Ed[10, :10] == brow_want
            zz = sp.expand(Ed[10, 10] + (v.T * v)[0]) == 0
            hblock = Ed[:10, :10] == Ep[:10, :10]
            eom = (Ric - v * v.T / 2); onsh = eom == sp.zeros(4, 4) and W.trace() == 0
            say(f"{kind} {i}: R = {Rs}, v^2 = {(v.T * v)[0]}, tr W = {W.trace()}, field equations hold at 0: {onsh};  s = {s}, second-derivative part: {second_ok}, no first-derivative terms: {first_zero}"
                + ("" if first_zero else f" (nonzero entries: {nfz})")
                + f";  E == printed: {eq_full}  [h-block {hblock}; (E phi)_ab = -2W: {mixed_col}; (E h)_phi = <Y,h>: {mixed_row}; E_phiphi = -v^2: {zz}]  [{time.time() - t0:.0f} s]")
            results.append((kind, R, Ric, Rs, v, W, Ed, Ep, s, second_ok, first_zero, eq_full, mixed_col and mixed_row and zz, onsh))
    on = [r for r in results if r[0] == 'on shell']; off = [r for r in results if r[0] == 'off shell']
    d1 = all(r[8] == 1 and r[9] and r[10] for r in results)
    d2 = all(r[11] and r[13] for r in on)
    d3 = all(r[12] for r in results)
    say(f"\nD1 (c = {args.c}) minimal: s = 1, second-derivative part Pbig delta, no first-derivative terms, every background: {'PASS' if d1 else 'FAIL'}")
    say(f"D2 on shell, derived E == printed E (11 x 11): {('PASS' if d2 else 'FAIL') if on else 'no on-shell background run'}")
    say(f"D3 mixed block and E_phiphi as the print's Y (4.32), Z (4.33), every background: {'PASS' if d3 else 'FAIL'}")
    d4 = all(r[11] for r in off)
    say(f"D4 off shell, derived E == printed E: {d4 if off else 'no off-shell background run'}")
    if not d4:
        # is the difference built from the field equations?  fit (derived - printed) on the h-block to the six structures with
        # Ric -> EOM = Ric - 1/2 v v, R -> tr EOM, plus  v v-type terms are not allowed
        for r in off:
            kind, R, Ric, Rs, v, W, Ed, Ep = r[:8]
            D = (Ed - Ep)[:10, :10]
            Em = Ric - v * v.T / 2; et = Em.trace()
            structs = {'trEOM*1': op_matrix10(lambda h: et * h), 'trEOM*dd': op_matrix10(lambda h: et * sp.eye(4) * h.trace()),
                       'EOM.h+h.EOM': op_matrix10(lambda h: Em * h + h * Em), 'd(x)EOM': op_matrix10(lambda h: sp.eye(4) * (Em * h).trace()),
                       'EOM(x)d': op_matrix10(lambda h: Em * h.trace())}
            names = list(structs); cs = sp.symbols('c0:5')
            Mx = sum((c * structs[nm] for c, nm in zip(cs, names)), sp.zeros(10, 10))
            sol = sp.solve([sp.Eq(Mx[i, j], D[i, j]) for i in range(10) for j in range(10)], cs, dict=True)
            say(f"   off shell: (derived - printed) h-block as field-equation structures: {({nm: sol[0].get(c, c) for c, nm in zip(cs, names)} if sol else 'NOT of this form')};  other blocks equal: {(Ed - Ep)[10, :] == sp.zeros(1, 11) and (Ed - Ep)[:, 10] == sp.zeros(11, 1)}")
    ok &= d1 and d2 and d3
    complete = bool(on) and bool(off) and not args.skip_d0
    if not ok:
        say("\nstatus=FAIL")
    elif not complete:
        say("\nstatus=INCOMPLETE (PASS needs D0 and at least one on-shell and one off-shell background)")
    else:
        say("\nstatus=PASS")
    say(f"END gravscalar_operator_derive {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
