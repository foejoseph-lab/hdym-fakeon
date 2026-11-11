#!/usr/bin/env python3
"""
Rung 3, step 3b (the improved tensor), the question of 10 Nov: sqrt(6) or sqrt(12)?

The rung 3 note (4 Nov, from the author's overnight s1 in hdym) records the field redefinition that removes a R phi^2 at a = 1/12 as
phi-tilde = sqrt(6) artanh(phi/sqrt(6)), F = 1 - phi^2/6.  A hand calculation on 10 Nov with 't Hooft-Veltman's own normalisation gave
sqrt(12) and F = 1 - phi^2/12.  This script derives it, for a general normalisation c of the curvature term:

    L_J = sqrt(g) ( c f R - 1/2 (d phi)^2 ),      f = 1 - a phi^2 / c      [i.e.  c R - a R phi^2 - 1/2 (d phi)^2;  R standard sign]
    L_E = sqrt(gE) ( c R[gE] - 1/2 K (d phi)^2_gE ),   gE = f g,   K = 1/f + 3 c f'^2 / f^2.

't Hooft-Veltman (7.1): sqrt(g)(-R - 1/2 (d phi)^2 + a R phi^2) with their R = minus the standard one: c = 1.
The common "reduced Planck mass = 1" normalisation sqrt(g)(R/2 - 1/2 (d phi)^2 - xi R phi^2 / 2), xi = 1/6:  c = 1/2, a = 1/12.

What is computed (exact rationals, python-flint polynomials through gravscalar_operator_derive's helpers; no formula for the Weyl
transformation of R is typed in -- R[f g] is computed from the metric f g):
  W1  pointwise, at the origin of a GENERIC polynomial metric (delta + random linear + random quadratic: not normal coordinates) with a
      generic polynomial phi:      L_E - L_J = - d_mu ( 3 c sqrt(g) g^{mu nu} d_nu f ),
      a total derivative, for several (a, c) including a != 1/12.   [The right side and K are my hand derivation of 10 Nov; W1 is
      the test of that derivation.  If W1 holds the two Lagrangians give the same field equations and the map is established.]
  W2  symbolic in phi:  K f^2 = 1 + (a phi^2 / c)(12 a - 1);   so K = 1/f^2 exactly at a = 1/12, for ANY c;
      the canonical field  phi-tilde = integral d phi / f = sqrt(12 c) artanh( phi / sqrt(12 c) ),   F = f = 1 - phi^2/(12 c).
PREDICTION (rule 25; in this file as first written, before the first run): W1 holds for every (a, c) tried; hence
  c = 1 ('t Hooft-Veltman): sqrt(12), F = 1 - phi^2/12;      c = 1/2: sqrt(6), F = 1 - phi^2/6.
  The sqrt(6) of the note is then the c = 1/2 normalisation, not 't Hooft-Veltman's.  (I could not check which normalisation the
  author's s1 used: it is in neither repository as fetched on 10 Nov.)

Run:  python3 weyl_map_check.py --log results/weyl_map_check.log        (seconds)

RECORD (rule 11; sandbox 10 Nov 18:40): status=PASS, 1 s.  W1 on ten cases, of which the two with a = 0 are 0 = 0: eight effective.
FRESH REVIEW (10 Nov, separate agent): no error found.  Mutations on a scratch copy -- K = 1/f alone, 3/2 for 3 in K, the sign or
the coefficient of the divergence term, g^-1 dropped from it, sqrt(gE) = f instead of f^2, the 1/f(0) on R dropped -- each fail on
all eight cases; sqrt(g) dropped from the divergence fails on seven; cubic terms added to g and phi change nothing (third derivatives
cannot enter).  Limits it asked to have stated:
  * W1 holds for ANY f, so it does not test f = 1 - a phi^2/c: the c in sqrt(12 c) rests on the one-line identification
    c f R = c R - a R phi^2 and on the 3c in K (tested at c = 1/2 and 3/2).
  * c = 1 for 't Hooft-Veltman is read from the rung 3 note's summary of (7.1), not from the scan; their (6.3), R_mn = -1/2 D phi D phi,
    has magnitude 1/(2c) and agrees.  That the note's sqrt(6) came from c = 1/2 is CONSISTENT WITH the algebra, no more: s1 cannot
    be consulted and a plain slip there is as possible.
  * The map needs f > 0, i.e. |phi| < sqrt(12 c); the script has m = 0 (their (7.2)-(7.7) carry a mass term, which becomes V/f^2);
    and the step from the map to a statement about the divergence assumes the paper's Jacobian-one remark (p. 70), itself open.
"""
import sys, os, time, argparse, random
import sympy as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gravscalar_operator_derive as D
from gravscalar_operator_derive import idx, XR, ZERO, ONE, C, d, psum, at0_rational, CTX

out = None


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def Tx(p, xdeg):
    return D.T(p, xdeg, 0, 0)


def rpoly(rng, lo, hi, size=3):
    """a random polynomial in x with all monomials of degree lo..hi"""
    import itertools
    tot = ZERO
    for m in itertools.product(range(hi + 1), repeat=4):
        if lo <= sum(m) <= hi:
            tot = tot + C(sp.Rational(rng.randint(-size, size), rng.randint(1, 3))) * XR[0]**m[0] * XR[1]**m[1] * XR[2]**m[2] * XR[3]**m[3]
    return tot


def mat_mul(A, B, xdeg):
    return [[Tx(psum(A[i][k] * B[k][j] for k in idx), xdeg) for j in idx] for i in idx]


def inverse(g, xdeg):
    """series inverse of a polynomial metric with g(0) = delta"""
    A = [[g[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]
    inv = [[ONE if i == j else ZERO for j in idx] for i in idx]; term = [row[:] for row in inv]
    for _ in range(xdeg + 2):
        term = mat_mul([[-t for t in row] for row in term], A, xdeg)
        inv = [[inv[i][j] + term[i][j] for j in idx] for i in idx]
    return inv


def sqrt_det(g, xdeg):
    A = [[g[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]
    tr = lambda M: psum(M[i][i] for i in idx)
    A2 = mat_mul(A, A, xdeg); A3 = mat_mul(A2, A, xdeg)
    ld = Tx(tr(A) - A2_half(tr(A2)) + tr(A3) * C(sp.Rational(1, 3)), xdeg)      # log det to third order in A (A = O(x))
    return Tx(ONE + ld * C(sp.Rational(1, 2)) + ld * ld * C(sp.Rational(1, 8)) + ld * ld * ld * C(sp.Rational(1, 48)), xdeg)


def A2_half(p):
    return p * C(sp.Rational(1, 2))


def ricci_scalar_at_0(g):
    """R(0) of a polynomial metric with g(0) = delta, standard sign (positive on spheres)"""
    ginv = inverse(g, 2)
    dg = [[[d(g[a][b], c) for c in idx] for b in idx] for a in idx]
    Gam = [[[Tx(psum(ginv[m][l] * (dg[l][b][a] + dg[l][a][b] - dg[a][b][l]) for l in idx) * C(sp.Rational(1, 2)), 1) for b in idx] for a in idx] for m in idx]
    Ric = [[Tx(psum(d(Gam[a][b][c], a) - d(Gam[a][b][a], c) for a in idx)
               + psum(Gam[a][a][e] * Gam[e][b][c] - Gam[a][c][e] * Gam[e][b][a] for a in idx for e in idx), 0) for c in idx] for b in idx]
    return at0_rational(Tx(psum(ginv[b][c] * Ric[b][c] for b in idx for c in idx), 0))


def w1(rng, a, c):
    # a generic metric with g(0) = delta (a linear change of coordinates makes any metric so at one point; the linear and quadratic
    # jets are unrestricted), a generic phi with phi(0) != 0
    g = [[None] * 4 for _ in idx]
    for i in idx:
        for j in idx:
            if i <= j:
                g[i][j] = g[j][i] = (ONE if i == j else ZERO) + rpoly(rng, 1, 2, 2) * C(sp.Rational(1, 4))
    phi = C(sp.Rational(rng.randint(1, 3), rng.randint(2, 5))) + rpoly(rng, 1, 2)
    f = Tx(ONE - C(sp.Rational(a) / c) * phi * phi, 2)
    f0 = at0_rational(Tx(f, 0)); phi0 = at0_rational(Tx(phi, 0))
    fp0 = -2 * sp.Rational(a) / c * phi0
    K0 = 1 / f0 + 3 * c * fp0**2 / f0**2
    ginv = inverse(g, 2); sg = sqrt_det(g, 2)
    dphi0 = [at0_rational(Tx(d(phi, m), 0)) for m in idx]                     # g^{-1}(0) = delta
    kin0 = sum(x * x for x in dphi0)
    LJ = c * f0 * ricci_scalar_at_0(g) - kin0 / 2                             # sqrt(g)(0) = 1
    # Einstein frame: rescale coordinates so that gE(0) = delta too:  gE = f g;  x = y / sqrt(f0) is avoided by computing R of the
    # metric  f g / f0  (a constant rescaling: R[lambda g] = R[g] / lambda), which has value delta at the origin
    gE_n = [[Tx(f * g[i][j] * C(1 / f0), 2) for j in idx] for i in idx]
    RE = ricci_scalar_at_0(gE_n) / f0
    LE = f0**2 * (c * RE - K0 * (kin0 / f0) / 2)                              # sqrt(gE) = f^2 sqrt(g);  gE^{-1} = g^{-1}/f
    V = [Tx(sg * psum(ginv[m][n] * d(f, n) for n in idx), 1) for m in idx]
    divV = sum(at0_rational(Tx(d(V[m], m), 0)) for m in idx)
    return sp.nsimplify(LE - LJ), -3 * c * divV, K0, f0


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261110)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0:
            print(f"refusing to overwrite the non-empty log {args.log}; to start over:  del {args.log.replace('/', chr(92))}"); sys.exit(1)
        out = open(args.log, 'w')
    say(f"# weyl_map_check: seed {args.seed}")
    say(D.version_line())
    rng = random.Random(args.seed)
    ok = True
    for a, c in ((sp.Rational(1, 12), 1), (sp.Rational(1, 12), sp.Rational(1, 2)), (sp.Rational(1, 5), 1), (sp.Rational(-2, 7), sp.Rational(3, 2)), (0, 1)):
        for rep in range(2):
            lhs, rhs, K0, f0 = w1(rng, a, c)
            good = sp.simplify(lhs - rhs) == 0
            ok &= good
            say(f"W1 a = {a}, c = {c}: L_E - L_J = {lhs};  -d(3 c sqrt(g) g d f) = {rhs};  equal: {good}   [f(0) = {f0}, K(0) = {K0}, K f^2 = {sp.nsimplify(K0 * f0**2)}]")
    phi, a, c = sp.symbols('phi a c', positive=True)
    f = 1 - a * phi**2 / c
    K = 1 / f + 3 * c * sp.diff(f, phi)**2 / f**2
    w2a = sp.simplify(K * f**2 - (1 + (a * phi**2 / c) * (12 * a - 1))) == 0
    Kc = sp.simplify(K.subs(a, sp.Rational(1, 12)) * f.subs(a, sp.Rational(1, 12))**2)
    cand = sp.sqrt(12 * c) * sp.atanh(phi / sp.sqrt(12 * c))
    w2b = sp.simplify(sp.diff(cand, phi) - 1 / f.subs(a, sp.Rational(1, 12))) == 0
    say(f"\nW2 K f^2 = 1 + (a phi^2/c)(12a - 1): {w2a};  at a = 1/12, K f^2 = {Kc} for any c;  d/dphi [sqrt(12c) artanh(phi/sqrt(12c))] = 1/f = sqrt(K): {w2b}")
    ok &= w2a and Kc == 1 and w2b
    say(f"   c = 1   ('t Hooft-Veltman (7.1)):      phi-tilde = sqrt(12) artanh(phi/sqrt(12)),  F = 1 - phi^2/12")
    say(f"   c = 1/2 (R/2, reduced Planck mass 1):  phi-tilde = sqrt(6)  artanh(phi/sqrt(6)),   F = 1 - phi^2/6")
    say(f"\nstatus={'PASS' if ok else 'FAIL'}")
    say(f"END weyl_map_check {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
