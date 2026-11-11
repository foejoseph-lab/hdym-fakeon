#!/usr/bin/env python3
"""
Rung 3, step 3b: THE PROVING RUN for a new method, before it is used on 't Hooft-Veltman's section 7 (the improved tensor).

The author's rule (10 Nov): a method that differs from the historical one gets its own column in the claims ledger, and for the
historical papers it must be PROVEN to give the same results as the method it stands beside before it is accepted; the easiest way
is an exact comparison on a simplified setting.

The new method ("conjugation").  If the fields of a quadratic Lagrangian are redefined linearly with a position-dependent matrix,
chi = J(x) phi, the fluctuation operator of the new Lagrangian is the old one conjugated,  Delta' = J^-1 Delta J,  so that
    connection  w'_m = J^-1 w_m J + J^-1 d_m J,     endomorphism  E' = J^-1 E J,     curvature  Omega' = J^-1 Omega J,
and every trace that enters the one-loop divergence is unchanged.  The method takes the divergence of the redefined theory from
the original one's, without deriving the redefined operator.  ('t Hooft-Veltman's "little theorem", their (3.6)-(3.7), is the
special case they prove; their own route to section 7 was to redefine fields and then compute directly, by Schoonschip.)

The method it must be proven against ("direct").  Derive the operator of the REDEFINED Lagrangian from that Lagrangian alone, with
no knowledge of J: second variation by jet probes (the method of graviton_operator_derive.py), here extended to give the operator's
coefficient FIELDS near the origin and not only their values at it.  For
    O phi = Kmn(x) d_m d_n phi + B_m(x) d_m phi + C(x) phi
the connection and endomorphism at the origin need K and B to first order in x and C at the origin.  The extension: translate the
background by a symbolic vector t kept to first order (every background field f(y + t), probes and base fields monomials in y),
so that the Euler-Lagrange operator "at the origin" is the operator at the point t.  From the measured fields:
    S = (1/4) g_mn Kmn  (the field-space metric times sqrt(g); the run checks Kmn = S g^mn, i.e. that the operator is of Laplace type),
    b = S^-1 B,  c = S^-1 C,     w^m = (1/2)(b^m + g^ab Gamma^m_ab),     E(0) = c - d_m w^m - w_m w_m,
    Omega_mn(0) = d_m w_n - d_n w_m + [w_m, w_n].
This extension is also what a direct computation of section 7 in the paper's own gauge would need, so this run is its first test.

The simplified settings (exact rationals; normal-coordinate metric on a random algebraic curvature tensor, as everywhere in rung 3):
  A  one scalar, L0 = sqrt(g)(-1/2 (d chi)^2 + 1/2 M(x) chi^2), chi = F(x) phi.        (abelian: Omega = 0; tests E and the first-derivative term)
  B  two scalars, M(x) a symmetric 2x2 matrix field, chi = J(x) phi, J a 2x2 polynomial matrix, not orthogonal.   (flat connection conjugated: Omega' = 0 must come out)
  C  a covector field, L0 = sqrt(g)(-1/2 D_m A_n D^m A^n + 1/2 A.M(x).A), A = J(x) phi, J 4x4.   (Omega = Riemann, non-abelian, conjugated)
J(x) = J0 (1 + X(x)) with J0 a random invertible rational matrix and X a random polynomial of degree 1-3 with X(0) = 0; M(x) a random
symmetric polynomial of degree <= 2.

PREDICTIONS (rule 25; written 10 Nov 20:35 ET, clock read, before any code below the docstring existed; nothing fetched):
  C0 [validation of the extension on the known case J = 1]: direct gives S = sqrt(g) 1 to first order, w = 0 (A, B) or the
      Levi-Civita connection on covectors (C; so Tr Omega^2 = -Riem^2, as the ghost bundle of step 2), and E(0) = M(0).
  C1 [the proof; derived by hand 10 Nov from O' = J^T O J]: for random J, direct gives, as matrices at the origin,
        w'_m = J^-1 w_m J + J^-1 d_m J,    E' = J^-1 E J,    Omega'_mn = J^-1 Omega_mn J,
      with the J = 1 quantities of C0 on the right.  In particular every second derivative of J cancels out of E' (c alone
      contains J^-1 Box J), and in setting B the curvature of the pure-gauge connection J^-1 dJ comes out zero.
  C2 [follows from C1]: Tr E', Tr E'^2, Tr Omega'^2 equal the J = 1 values, so the (3.40) number is unchanged.
  C3 [the recursion with a connection that is not the metric's; not exercised anywhere before]: DeWitt's recursion
      (rung3_dewitt_bundle.py's functions, imported, flint backend) run with the conjugated connection and endomorphism AS FIELDS
      gives Tr a_2(0) equal to the J = 1 value exactly, pointwise; and both equal 2 x (3.40) + (d/180) GB + (1/6) Box(Tr E)(0).
  Controls (must FAIL): --mutate half  halves the J^-1 dJ term in the conjugation formula;  --mutate esym  uses J^T E J for E'.
What a PASS will and will not mean.  It confirms, on random instances of these three settings, that conjugation and direct
derivation give the same operator data, and hence (C2) the same (3.40) number.  It does not prove it for a field redefinition that mixes the metric with the fields (the
Weyl map of section 7 also rescales g): there the same statement needs the gauge-fixing term carried along, and that is the next
run, not this one.  And it says nothing about finite parts.
Status PASS if C0-C3 hold in all three settings on every background.

WHAT KIND OF PROOF THIS IS (asked by the author, 10 Nov 20:34: "if it's rigor it's rigor").  Two different things, kept apart:
  (i) The mathematical statement has a rigorous proof on paper, three lines, no physics.  For chi = J phi the second variation is
      O' = J^T O J and the field-space metric is S' = J^T S J (the chain rule for a linear substitution; exact).  Hence
      Delta' = -S'^-1 O' = J^-1 Delta J.  Similar operators have conjugate heat kernels, K'(x, y) = J(x)^-1 K(x, y) J(y), so the
      diagonal traces agree at every x and every order: Tr a_k'(x) = Tr a_k(x).  That is the theorem; it does not need this run.
  (ii) This run is NOT a proof of the theorem.  It is an exact test that the CODE on both sides implements it: the direct
      derivation (which knows nothing of J) and the conjugation formulas give the same matrices in exact rational arithmetic on
      random instances.  The quantities are rational functions of finitely many jets, so agreement on random rational instances
      is polynomial identity testing: overwhelming evidence, not certainty -- an identity that failed would have to vanish on
      every instance drawn.  It is not Lean and not a symbolic proof.
  So: the method is rigorous by (i); what (ii) establishes is that our implementation of it, and of the direct method it stands
  beside, agree.  If the author wants (ii) raised to a proof in the simplified setting, setting A can be run with symbolic jets.

Run:  python3 conjugation_proving.py --jobs 6 --log results/conjugation_proving.log
      controls (must FAIL):  --settings B --backgrounds 1 --mutate half   and   --mutate esym

FRESH REVIEW (10 Nov night, separate agent; ~40 mutations on scratch copies; repository untouched).  The mathematics above is right
and every instance run passes (A, B on two seeds, C once; pool under a forced spawn start method identical to serial).  Folded in:
  * MUST-FIX: the first version passed its matrices through sympy's nsimplify, which is NOT the identity on exact rationals (it
    snaps values within ~1e-7 and now and then returns an irrational product), so "exact" was untrue and a correct run could have
    failed.  Removed everywhere; all entries are exact Rationals and are compared as such.  No number of the smoke changed.
  * What each check carries (the review's mutation table):  C1 and C2 test COVARIANCE under the redefinition -- they catch a wrong
    1/2 in w, a dropped S1 or w.w term, a wrong sign of d.w in E, c without S^-1, a wrong Omega (in B, C), an unshifted J;  they do
    NOT catch an error that is itself covariant (the Gamma term of w, a dropped sqrt(g)): those are caught by C0, which carries
    the absolute normalisation at J = 1 (now strengthened: S = 1 to first order, Omega zero as matrices in A, B, and equal to
    Riemann elementwise with its sign in C -- an overall sign of Omega passed every check before).  At J = 1 the data is
    degenerate (w = 0, S1 = 0), so the S1 term, the w.w term and the commutator are exercised only through the redefined side.
  * C3 as first written tested less than it said: Tr a_2(0) depends only on Omega(0), E(0) and Box Tr E, so dropping J^-1 dJ
    altogether, or leaving E unconjugated, still passed; and in C the covector connection's sign and the form of E were invisible.
    What C3 does test is the recursion with N(0) != 0 (the generalised transport: a dropped degree-1 term, a truncated a_0, a wrong
    J^-1 dJ in B, C are caught).  Added: the fields fed to the recursion must have, at the origin, exactly the w, E and Omega the
    direct method derived, for both the original and the conjugated operator -- this ties C3 to the direct side and sees the sign.
  * Unstated assumption, now stated: the local Jacobian det J of the field redefinition is dropped, as in dimensional regularization.
  * The instances are small (numerators -3..3, denominators 1..3): "overwhelming" is fair in practice, weak as a formal bound.
  * One stream of random numbers per setting, so --settings C alone draws the same instances as the full run; cases are printed as
    they finish.  The reference to 't Hooft-Veltman (3.6)-(3.7) is from the rung 3 note's summary, not re-read.
SMOKE RECORD (rule 11; sandbox 10 Nov 20:50-21:25 ET, one background per setting, single core): A 19 s, B 113 s, C 1015 s; C0-C3
hold in all three; both controls on B fail at C1 as intended (half: the connection; esym: the endomorphism).  Passed on the first
execution; the per-case pool was added afterwards and changes no number.  After the review's fold-in the instances differ (one
random stream per setting); re-run with the strengthened checks, sandbox 10 Nov, logs' file times 21:47 and 22:08 ET: A 24 s, B 140 s (two workers), C 1250 s
(alone), every check True including Omega == Riemann elementwise in C and the recursion's input fields against the direct data.
"""
import sys, os, time, argparse, itertools, random, platform
import sympy as sp
import flint

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'graviton_loop'))
idx = range(4)
out = None
CTX = flint.fmpq_mpoly_ctx.get(('y1', 'y2', 'y3', 'y4', 't1', 't2', 't3', 't4', 'dl'), 'degrevlex')
GEN = CTX.gens(); Y = GEN[:4]; TT = GEN[4:8]; DL = GEN[8]
ZERO = CTX.from_dict({}); ONE = CTX.constant(flint.fmpq(1))
MONS = {dg: [m for m in itertools.product(range(dg + 1), repeat=4) if sum(m) == dg] for dg in range(5)}


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def q(r):
    r = sp.Rational(r); return flint.fmpq(int(r.p), int(r.q))


def C(r):
    return CTX.constant(q(r))


def T(p, ydeg):
    """keep y-degree <= ydeg, t-degree <= 1, dl-degree <= 1"""
    return CTX.from_dict({m: c for m, c in p.to_dict().items() if sum(m[:4]) <= ydeg and sum(m[4:8]) <= 1 and m[8] <= 1})


def psum(it):
    tot = ZERO
    for e in it:
        tot = tot + e
    return tot


def ypoly(coeffs):
    """{monomial (4): Rational} -> polynomial in y"""
    return psum(C(c) * Y[0]**m[0] * Y[1]**m[1] * Y[2]**m[2] * Y[3]**m[3] for m, c in coeffs.items() if c != 0)


def shift(p, ydeg):
    """p(y) -> p(y + t) to first order in t"""
    return T(p + psum(TT[m] * p.derivative(m) for m in idx), ydeg)


def rat(c):
    return sp.Rational(int(c.p), int(c.q))


def t_expand(p):
    """a polynomial in t only (y = 0 part of p), as (value at t = 0, [d/dt_c]) of sympy Rationals"""
    v = sp.Integer(0); g = [sp.Integer(0)] * 4
    for m, c in p.to_dict().items():
        if sum(m[:4]) or m[8]:
            continue
        if sum(m[4:8]) == 0:
            v = rat(c)
        else:
            g[list(m[4:8]).index(1)] = rat(c)
    return v, g


# ----------------------------------------------------------------- random data (plain rationals; built into either ring later)
def rr(rng, lo=-3, hi=3, den=3):
    return sp.Rational(rng.randint(lo, hi), rng.randint(1, den))


def sample(rng, d):
    while True:
        J0 = sp.Matrix(d, d, lambda i, j: sp.Rational(rng.randint(-2, 2), rng.randint(1, 2)) + (2 if i == j else 0))
        if J0.det() != 0:
            break
    X = {m: sp.Matrix(d, d, lambda i, j: rr(rng) / 2) for dg in (1, 2, 3) for m in MONS[dg]}
    M = {}
    for dg in (0, 1, 2):
        for m in MONS[dg]:
            A = sp.Matrix(d, d, lambda i, j: rr(rng)); M[m] = (A + A.T) / 2
    return J0, X, M


# ----------------------------------------------------------------- the direct method: the operator of the redefined Lagrangian
class Direct:
    def __init__(self, R, setting, J0, X, M):
        self.setting = setting; d = self.d = J0.shape[0]
        gy = [[None] * 4 for _ in idx]
        for mu in idx:
            for nu in idx:
                q2 = psum(C(-sp.Rational(1, 3) * R[(mu, a, nu, b)]) * Y[a] * Y[b] for a in idx for b in idx)
                q4 = psum(C(sp.Rational(2, 45) * R[(a, mu, b, l)] * R[(l, c, nu, e)]) * Y[a] * Y[b] * Y[c] * Y[e]
                          for a in idx for b in idx for c in idx for e in idx for l in idx)
                gy[mu][nu] = (ONE if mu == nu else ZERO) + q2 + q4
        g = [[shift(gy[a][b], 3) for b in idx] for a in idx]                      # g(y + t), to y^3, first order in t
        A = [[g[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]
        mm = lambda P, Q: [[T(psum(P[i][k] * Q[k][j] for k in idx), 3) for j in idx] for i in idx]
        A2 = mm(A, A)
        self.ginv = [[(ONE if i == j else ZERO) - A[i][j] + A2[i][j] for j in idx] for i in idx]      # A = O(2) in (y, t): A^3 = O(6), beyond y^3 t^1
        tr = lambda P: psum(P[i][i] for i in idx)
        ld = T(tr(A) - tr(A2) * C(sp.Rational(1, 2)), 3)
        self.sg = T(ONE + ld * C(sp.Rational(1, 2)) + ld * ld * C(sp.Rational(1, 8)), 3)
        dg = [[[g[a][b].derivative(c) for c in idx] for b in idx] for a in idx]
        self.Gam = [[[T(psum(self.ginv[l][k] * (dg[k][n][m] + dg[k][m][n] - dg[m][n][k]) for k in idx) * C(sp.Rational(1, 2)), 2)
                      for n in idx] for m in idx] for l in idx]                  # Gam[l][m][n] = Gamma^l_{mn}
        Xy = [[ypoly({m: X[m][i, j] for m in X}) for j in range(d)] for i in range(d)]
        Jy = [[psum(C(J0[i, k]) * ((ONE if k == j else ZERO) + Xy[k][j]) for k in range(d)) for j in range(d)] for i in range(d)]
        self.J = [[shift(Jy[i][j], 3) for j in range(d)] for i in range(d)]
        self.M = [[shift(ypoly({m: M[m][i, j] for m in M}), 2) for j in range(d)] for i in range(d)]

    def lagrangian(self, phi):
        """the quadratic Lagrangian of the redefined fields, to y^2: L0(chi = J phi)"""
        d = self.d
        chi = [T(psum(self.J[i][j] * phi[j] for j in range(d)), 3) for i in range(d)]
        dchi = [[chi[i].derivative(m) for m in idx] for i in range(d)]
        if self.setting == 'C':            # chi is a covector: D_m chi_n = d_m chi_n - Gamma^l_{mn} chi_l
            Dc = [[T(dchi[n][m] - psum(self.Gam[l][m][n] * chi[l] for l in idx), 2) for n in idx] for m in idx]
            up = [[T(psum(self.ginv[m][p] * Dc[p][n] for p in idx), 2) for n in idx] for m in idx]
            up2 = [[T(psum(self.ginv[n][q_] * up[m][q_] for q_ in idx), 2) for n in idx] for m in idx]
            kin = T(psum(Dc[m][n] * up2[m][n] for m in idx for n in idx), 2)
        else:
            kin = T(psum(self.ginv[m][n] * dchi[i][m] * dchi[i][n] for i in range(d) for m in idx for n in idx), 2)
        pot = T(psum(chi[i] * self.M[i][j] * chi[j] for i in range(d) for j in range(d)), 2)
        return T(self.sg * (pot - kin) * C(sp.Rational(1, 2)), 2)

    def el(self, base):
        """Euler-Lagrange operator on the base field at y = 0, i.e. at the point t: d polynomials in t"""
        d = self.d
        G = {}
        for k in range(d):
            for kind, key, mono in PROBES:
                phi = [base[i] + (DL * mono if i == k else ZERO) for i in range(d)]
                L = self.lagrangian(phi)
                G[(k, kind, key)] = CTX.from_dict({m[:8] + (0,): c for m, c in L.to_dict().items() if m[8] == 1})
        half = C(sp.Rational(1, 2)); res = []
        for k in range(d):
            G0 = G[(k, 0, None)]; G1 = [G[(k, 1, m)] for m in idx]
            expr = G0
            for m in idx:
                expr = expr - (G1[m] - Y[m] * G0).derivative(m)
            for m in idx:
                for n in idx:
                    w = ONE if m == n else half
                    F2 = Y[m] * Y[n] * G0 * half - (Y[n] * G1[m] + Y[m] * G1[n]) * half + w * G[(k, 2, (min(m, n), max(m, n)))]
                    expr = expr + F2.derivative(m).derivative(n)
            res.append(T(expr, 0))
        return res

    def data(self):
        """S, w, E, Omega at the origin from the coefficient fields K, B, C measured to first order in t"""
        d = self.d
        def col_fields(mono):
            cols = []
            for k in range(d):
                r = self.el([mono if i == k else ZERO for i in range(d)])
                cols.append([t_expand(x) for x in r])
            v0 = sp.Matrix(d, d, lambda i, j: cols[j][i][0])
            v1 = [sp.Matrix(d, d, lambda i, j, c=c: cols[j][i][1][c]) for c in idx]
            return v0, v1
        C0, C1 = col_fields(ONE)
        Bf = [col_fields(Y[m]) for m in idx]
        Kf = {(m, n): col_fields(Y[m] * Y[n] * (C(sp.Rational(1, 2)) if m == n else ONE)) for m in idx for n in idx if m <= n}
        S0 = sum((Kf[(m, m)][0] for m in idx), sp.zeros(d, d)) / 4
        S1 = [sum((Kf[(m, m)][1][c] for m in idx), sp.zeros(d, d)) / 4 for c in idx]
        laplace = all(Kf[(m, m)][0] == S0 and all(Kf[(m, m)][1][c] == S1[c] for c in idx) for m in idx) and \
            all(Kf[(m, n)][0] == sp.zeros(d, d) and all(Kf[(m, n)][1][c] == sp.zeros(d, d) for c in idx) for (m, n) in Kf if m != n)
        Si = S0.inv()
        # gamma^m = g^{ab} Gamma^m_{ab}: zero at the origin; its first derivative d_c gamma^m(0) from the t-linear part
        gam1 = [[t_expand(T(psum(self.Gam[m][a][a] for a in idx), 0))[1][c] for m in idx] for c in idx]     # gam1[c][m]
        w0 = [Si * Bf[m][0] / 2 for m in idx]
        w1 = [[(Si * Bf[m][1][c] - Si * S1[c] * Si * Bf[m][0]) / 2 + gam1[c][m] * sp.eye(d) / 2 for m in idx] for c in idx]   # w1[c][m] = d_c w^m
        c0 = Si * C0
        E = c0 - sum((w1[m][m] for m in idx), sp.zeros(d, d)) - sum((w0[m] * w0[m] for m in idx), sp.zeros(d, d))
        Om = {(m, n): w1[m][n] - w1[n][m] + w0[m] * w0[n] - w0[n] * w0[m] for m in idx for n in idx}
        return dict(S0=S0, S1=S1, laplace=laplace, w=w0, E=E, Om=Om)          # all entries exact Rationals; nothing is simplified numerically


PROBES = [(0, None, ONE)] + [(1, m, Y[m]) for m in idx] + \
         [(2, (m, n), Y[m] * Y[n] * (C(sp.Rational(1, 2)) if m == n else ONE)) for m in idx for n in idx if m <= n]


# ----------------------------------------------------------------- the recursion with a general connection (C3)
def recursion_trace_a2(B, R, setting, J0, X, M, conjugate):
    """Tr a_2(0) by DeWitt's recursion (the functions of rung3_dewitt_bundle.py; flint backend) for -(D^2 + E) with connection Nc
    and endomorphism Ef given as fields; conjugate = False: the original operator, True: conjugated by J(x)."""
    d = J0.shape[0]; XB = B.X
    g, ginv, sqrtg, isqrtg, Dhalf, iDhalf = B.geometry(R)
    G = B.christoffel(g, ginv)
    fld = lambda co: sum(c * XB[0]**m[0] * XB[1]**m[1] * XB[2]**m[2] * XB[3]**m[3] for m, c in co.items() if c != 0) + 0 * XB[0]
    Mf = [[fld({m: M[m][i, j] for m in M}) for j in range(d)] for i in range(d)]
    if setting == 'C':                                   # covector: (D_mu f)_a = d_mu f_a - Gamma^b_{mu a} f_b;  E = g.M as an endomorphism
        Nc = [[[-1 * G[b][mu][a] for b in idx] for a in idx] for mu in idx]
        Ef = B.mmul(g, Mf, 2)
    else:
        Nc = [[[0 * XB[0] for _ in range(d)] for _ in range(d)] for _ in idx]
        Ef = B.mmap(lambda e: B.trunc(e, 2), Mf)
    if conjugate:
        Xf = [[fld({m: X[m][i, j] for m in X}) for j in range(d)] for i in range(d)]
        one = [[(1 if i == j else 0) + 0 * XB[0] for j in range(d)] for i in range(d)]
        J0f = [[J0[i, j] + 0 * XB[0] for j in range(d)] for i in range(d)]
        J0i = J0.inv(); J0if = [[J0i[i, j] + 0 * XB[0] for j in range(d)] for i in range(d)]
        Jf = B.mmul(J0f, B.madd(one, Xf), 4)
        acc = one; term = one
        for _ in range(4):                               # (1 + X)^-1 = sum (-X)^k, X = O(x): exact to x^4
            term = B.mmul(B.mscale(term, -1), Xf, 4); acc = B.madd(acc, term)
        Ji = B.mmul(acc, J0if, 4)
        Nc = [B.madd(B.mmul(B.mmul(Ji, Nc[mu], 4), Jf, 4), B.mmul(Ji, B.mmap(lambda e, mu=mu: B.pd(e, mu), Jf), 4)) for mu in idx]
        Ef = B.mmul(B.mmul(Ji, Ef, 2), Jf, 2)

    def transport(k, RHS, maxdeg, init=None):
        """rung3_dewitt_bundle.transport with the degree-1 part of x.N allowed (N(0) != 0): the same recursion, one more term"""
        xN = B.sum_mat([B.mmap(lambda e, mu=mu: B.trunc(XB[mu] * e, maxdeg + 1), Nc[mu]) for mu in idx])
        xN_h = {m: B.mmap(lambda e, m=m: B.homog(e, m), xN) for m in range(1, maxdeg + 1)}
        P = {}
        for dd in range(maxdeg + 1):
            rhs = B.mmap(lambda e, dd=dd: B.homog(e, dd), RHS) if RHS is not None else B.mzero(d)
            for m, Qm in xN_h.items():
                if dd - m in P:
                    rhs = B.madd(rhs, B.mmul(Qm, P[dd - m], maxdeg), -1)
            P[dd] = init if dd + k == 0 else B.mscale(rhs, sp.Rational(1, dd + k))
        return B.sum_mat([P[dd] for dd in P])

    a0 = transport(0, None, 4, init=B.meye(d))
    f0 = B.mmap(lambda e: B.trunc(Dhalf * e, 4), a0)
    rhs1 = B.madd(B.box_bundle(f0, Nc, ginv, sqrtg, isqrtg, 2), B.mmul(Ef, f0, 2))
    rhs1 = B.mmap(lambda e: B.trunc(iDhalf * e, 2), rhs1)
    a1 = transport(1, rhs1, 2)
    f1 = B.mmap(lambda e: B.trunc(Dhalf * e, 2), a1)
    rhs2 = B.madd(B.box_bundle(f1, Nc, ginv, sqrtg, isqrtg, 0), B.mmul(Ef, f1, 0))
    tr_a2 = sp.expand(sum(B.at0(rhs2[i][i]) for i in range(d)) / 2)
    trE = sum(Ef[i][i] for i in range(d))
    boxTrE = sp.expand(B.at0(sum(B.pd(B.pd(trE, m), m) for m in idx)))
    N0 = [sp.Matrix(d, d, lambda i, j, mu=mu: B.at0(Nc[mu][i][j])) for mu in idx]
    dN = [[sp.Matrix(d, d, lambda i, j, m=m, n=n: B.at0(B.pd(Nc[n][i][j], m))) for n in idx] for m in idx]          # dN[m][n] = d_m N_n (0)
    Om0 = {(m, n): dN[m][n] - dN[n][m] + N0[m] * N0[n] - N0[n] * N0[m] for m in idx for n in idx}
    E0 = sp.Matrix(d, d, lambda i, j: B.at0(Ef[i][j]))
    return tr_a2, boxTrE, (N0, E0, Om0)


def run_case(job):
    """one (setting, background): returns its log lines and verdict (so a pool can run the cases side by side)"""
    setting, ib, R, J0, X, M, mutate = job
    import rung3_dewitt_bundle as B
    import hv1974_gravscalar_M as GM
    B.use_backend('flint')
    d = J0.shape[0]; lines = []
    Ric, Rs, Ric2, Riem2 = B.invariants(R)
    class A_: pass
    args = A_(); args.mutate = mutate
    tb = time.time()
    ref = Direct(R, setting, sp.eye(d), {m: sp.zeros(d, d) for m in X}, M).data()        # J = 1
    new = Direct(R, setting, J0, X, M).data()                                             # the redefined Lagrangian; J not used below this line except in the conjugation formulas
    M0 = M[(0, 0, 0, 0)]
    Om2 = lambda D_: sum((D_['Om'][(m, n)] * D_['Om'][(m, n)]).trace() for m in idx for n in idx)
    s_ok = ref['S0'] == sp.eye(d) and all(x == sp.zeros(d, d) for x in ref['S1'])                 # S = sqrt(g) x fibre metric = 1 to first order
    if setting == 'C':            # Levi-Civita on covectors, elementwise and with its sign: (Omega_mn)_ab = R_abmn
        om_ok = all(ref['Om'][(m, n)][a, b] == R[(a, b, m, n)] for m in idx for n in idx for a in idx for b in idx)
    else:
        om_ok = all(v == sp.zeros(d, d) for v in ref['Om'].values())
    c0 = ref['laplace'] and s_ok and ref['E'] == M0 and all(w == sp.zeros(d, d) for w in ref['w']) and om_ok
    Ji = J0.inv()
    X1 = [X[tuple(1 if i == m else 0 for i in idx)] for m in idx]                         # J^-1 d_m J at the origin = d_m X(0)
    fac = sp.Rational(1, 2) if args.mutate == 'half' else 1
    w_conj = [Ji * ref['w'][m] * J0 + fac * X1[m] for m in idx]
    E_conj = (J0.T if args.mutate == 'esym' else Ji) * ref['E'] * J0
    Om_conj = {k: Ji * v * J0 for k, v in ref['Om'].items()}
    c1w = all(new['w'][m] == w_conj[m] for m in idx)
    c1E = new['E'] == E_conj
    c1O = all(new['Om'][k] == Om_conj[k] for k in Om_conj)
    c1 = new['laplace'] and new['S0'] == (J0.T * ref['S0'] * J0) and c1w and c1E and c1O
    c2 = new['E'].trace() == ref['E'].trace() and (new['E']**2).trace() == (ref['E']**2).trace() and Om2(new) == Om2(ref)
    td = time.time() - tb
    a2_ref, box_ref, (Nr, Er, Or) = recursion_trace_a2(B, R, setting, J0, X, M, False)
    a2_new, box_new, (Nn, En, On) = recursion_trace_a2(B, R, setting, J0, X, M, True)
    same = lambda D_, N_, E_, O_: all(D_['w'][m] == N_[m] for m in idx) and D_['E'] == E_ and all(D_['Om'][k] == O_[k] for k in O_)
    c3in = same(ref, Nr, Er, Or) and same(new, Nn, En, On)      # the fields fed to the recursion carry, at the origin, what the direct method derived
    GBv = Riem2 - 4 * Ric2 + Rs**2
    pred = sp.expand(2 * GM.hv340(d, Om2(ref), ref['E'], Rs, Ric2) + sp.Rational(d, 180) * GBv + box_ref / 6)
    c3 = c3in and sp.expand(a2_new - a2_ref) == 0 and sp.expand(a2_ref - pred) == 0 and sp.expand(box_new - box_ref) == 0
    lines.append(f"{setting}{ib} (d = {d}; R = {Rs}; det J(0) = {J0.det()}):")
    lines.append(f"   C0 J = 1: Laplace type {ref['laplace']}, w(0) = 0 {all(w == sp.zeros(d, d) for w in ref['w'])}, E(0) = M(0) {ref['E'] == M0}, S = 1 to first order {s_ok}, Omega " + ("== Riemann elementwise" if setting == 'C' else "== 0 as matrices") + f" {om_ok}, Tr Omega^2 = {Om2(ref)}" + (f" (-Riem^2 = {-Riem2})" if setting == 'C' else "") + f"  -> {'PASS' if c0 else 'FAIL'}")
    lines.append(f"   C1 direct == conjugated: S {new['S0'] == (J0.T * ref['S0'] * J0)}, Laplace type {new['laplace']}, w {c1w}, E {c1E}, Omega {c1O}  -> {'PASS' if c1 else 'FAIL'}   [w'(0) nonzero: {any(w != sp.zeros(d, d) for w in new['w'])}]")
    lines.append(f"   C2 traces: Tr E {new['E'].trace()} / {ref['E'].trace()}, Tr E^2 {(new['E']**2).trace()} / {(ref['E']**2).trace()}, Tr Omega^2 {Om2(new)} / {Om2(ref)}  -> {'PASS' if c2 else 'FAIL'}")
    lines.append(f"   C3 recursion: its input fields == the direct method's w, E, Omega at the origin (original and conjugated) {c3in}; Tr a_2(0) conjugated {a2_new}, original {a2_ref}, 2(3.40) + (d/180)GB + (1/6)Box Tr E = {pred}  -> {'PASS' if c3 else 'FAIL'}   [direct {td:.0f} s, total {time.time() - tb:.0f} s]")
    return lines, bool(c0 and c1 and c2 and c3)


def version_line():
    return f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}, python-flint {flint.__version__}; {sys.executable}"


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261110)
    ap.add_argument('--backgrounds', type=int, default=2)
    ap.add_argument('--settings', default='A,B,C')
    ap.add_argument('--mutate', default='none', choices=['none', 'half', 'esym'], help='controls: break the conjugation formula (must FAIL)')
    ap.add_argument('--jobs', type=int, default=1, help='cases run side by side (there are settings x backgrounds of them; the covector cases dominate)')
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# conjugation_proving: seed {args.seed}, settings {args.settings}, {args.backgrounds} backgrounds each, mutate {args.mutate}, jobs {args.jobs}")
    say(version_line())
    import rung3_dewitt_bundle as B
    DIM = {'A': 1, 'B': 2, 'C': 4}
    jobs = []
    for setting in args.settings.split(','):
        rng = random.Random(args.seed + 1000 * 'ABC'.index(setting))        # one stream per setting: a subset of settings draws the same instances
        for ib in range(args.backgrounds):
            R = B.riemann_sample(rng); J0, X, M = sample(rng, DIM[setting])
            jobs.append((setting, ib, R, J0, X, M, args.mutate))
    ok = True
    if args.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(min(args.jobs, len(jobs))) as pool:
            it = pool.imap(run_case, jobs, chunksize=1)          # in order; each case is printed as soon as it and those before it are done
            for lines, good in it:
                for ln in lines:
                    say(ln)
                ok &= good
    else:
        for j in jobs:
            lines, good = run_case(j)
            for ln in lines:
                say(ln)
            ok &= good
    say(f"\nstatus={'PASS' if ok else 'FAIL'}")
    say(f"END conjugation_proving {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
