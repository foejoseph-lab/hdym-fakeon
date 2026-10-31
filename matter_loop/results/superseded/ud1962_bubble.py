#!/usr/bin/env python3
"""
Rung 1, technique (a): the 1962 computation, computer-aided.

One real scalar (mass m, non-minimal coupling xi) on g = delta + h, Euclidean, flat-space perturbation
theory. The one-loop graviton self-energy from the single scalar bubble (two h-vertices), done by the
Feynman-parameter method (P) with a hard momentum cutoff Lambda on the shifted loop momentum, as a 1962
calculation would have been regulated. The Lambda^4, Lambda^2, log Lambda pieces are separated; the
log Lambda piece at order p^4 is projected onto the two quadratic-curvature structures and compared with
't Hooft-Veltman 1974 eq. (3.40) (Ann. IHP A 20, 69), which says, in units sqrt(g)/(8 pi^2 (n-4)):
    xi = 0   : (1/120) R_{mu nu} R^{mu nu} + (1/240) R^2
    xi = 1/6 : (1/120) (R_{mu nu} R^{mu nu} - R^2/3)
Dictionary cutoff <-> dim reg: int d^nk/(2pi)^n (k^2+Delta)^-2 = (1/16pi^2)(2/(4-n) + ...) in dim reg and
(1/16pi^2)(log(Lambda^2/Delta) - 1) with a cutoff, so (1/16pi^2) log Lambda^2 <-> -1/(8pi^2(n-4)) = -1/eps.
Sign: 't Hooft-Veltman write the Lagrangian as L = -1/2 (d phi)^2 + ... with a Euclidean metric, i.e.
L_theirs = -L_E (ours), and Delta L is added to L_theirs to cancel the divergence; so Delta L_theirs =
+Gamma_div (ours).  Hence if the coefficient of (1/16pi^2) log Lambda^2 in Gamma_2 is c, then
Gamma_div = -c/eps and the 't Hooft-Veltman coefficient is  -c.
Expected therefore:  A = -1/120,  B = -[ (1/4)(xi - 1/6)^2 - 1/360 ] = -(60 xi^2 - 20 xi + 1)/240,
the R^2 coefficient of (3.40) being (1/4)(M - R/6)^2 -> (1/4)(xi - 1/6)^2 R^2 plus -R^2/360 from the
(1/120)(Ric^2 - R^2/3) term.  At xi = 1/6 the square vanishes and A Ric^2 + B R^2 = -(1/120)(Ric^2 - R^2/3),
the Weyl^2 structure modulo Gauss-Bonnet.

Not in this run: the tadpole (seagull) graph. It carries no p^4 and so cannot touch the target; it is
where the Lambda^4 and the bulk of the Lambda^2 terms of the 1962 abstract live. Next step.

Conventions. S = 1/2 int sqrt(g) [ g^{mu nu} d_mu phi d_nu phi + m^2 phi^2 + xi R phi^2 ].
sqrt(g) g^{mu nu} = delta - h_{mu nu} + (1/2) h delta + O(h^2);  sqrt(g) = 1 + h/2;
sqrt(g) R = R1 + O(h^2),  R1 = d_mu d_nu h_{mu nu} - d^2 h   (R = +2 on the unit S^2).
Fourier: f(x) = int_p f(p) e^{ipx}.  First-order vertex (S_1 = 1/2 int phi(k) phi(q) h_{mu nu}(p) V_{mu nu}(k,q;p),
k+q+p = 0):
    V_{mu nu}(k,q;p) = (k_mu q_nu + k_nu q_mu)/2 - (1/2) delta_{mu nu} (k.q - m^2) + xi (p^2 delta_{mu nu} - p_mu p_nu).
Gamma_2 = int_p (1/2) h_{mu nu}(p) Pi_{mu nu rho sigma}(p) h_{rho sigma}(-p),
    Pi = -(1/2) int d^4k/(2pi)^4  V_{mu nu}(-k, k-p; p) V_{rho sigma}(p-k, k; -p) / ((k^2+m^2)((k-p)^2+m^2)).
Quadratic forms to match (momentum space, p -> i p):
    R1(p)        = -p_mu p_nu h_{mu nu} + p^2 h
    R1_{mu nu}(p) = (1/2)( -p_a p_nu h_{a mu} - p_a p_mu h_{a nu} + p_mu p_nu h + p^2 h_{mu nu} )
    int (R1)^2 -> int_p R1(p) R1(-p),  int (R1_{mu nu})^2 -> int_p R1_{mu nu}(p) R1_{mu nu}(-p).

Usage: python3 ud1962_bubble.py [--log results/ud1962_bubble.log]
"""
import sys, time, itertools, argparse
import sympy as sp

t0 = time.time()
ap = argparse.ArgumentParser()
ap.add_argument('--log', default=None)
args = ap.parse_args()
out = open(args.log, 'w') if args.log else None
def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s);
    if out: out.write(s + '\n'); out.flush()

say('# ud1962_bubble: 1962 cutoff bubble, method P, scalar on flat background')

# ---------- symbols ----------
m, xi, Lam, x = sp.symbols('m xi Lambda x', positive=True)
P = sp.symbols('p1:5', real=True)          # external momentum components
L = sp.symbols('l1:5', real=True)          # shifted loop momentum components
H = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{min(i,j)+1}{max(i,j)+1}', real=True))  # symmetric h
p2 = sum(pi**2 for pi in P)
l2 = sum(li**2 for li in L)
delta = sp.eye(4)

def dot(a, b):
    return sum(ai*bi for ai, bi in zip(a, b))

def V(k, q, p):
    """first-order vertex V_{mu nu}(k,q;p) as a 4x4 sympy Matrix"""
    kq = dot(k, q)
    pp = dot(p, p)
    return sp.Matrix(4, 4, lambda i, j:
        sp.Rational(1, 2)*(k[i]*q[j] + k[j]*q[i])
        - sp.Rational(1, 2)*delta[i, j]*(kq - m**2)
        + xi*(pp*delta[i, j] - p[i]*p[j]))

# ---------- Feynman parametrization ----------
# 1/((k^2+m^2)((k-p)^2+m^2)) = int_0^1 dx 1/[(l^2 + Delta)^2],  l = k - x p,  Delta = m^2 + x(1-x) p^2
Delta = m**2 + x*(1 - x)*p2
k = [li + x*pi for li, pi in zip(L, P)]            # k = l + x p
kmp = [ki - pi for ki, pi in zip(k, P)]            # k - p
V1 = V([-ki for ki in k], kmp, P)                  # V(-k, k-p; p)
V2 = V([-ki for ki in kmp], k, [-pi for pi in P])  # V(p-k, k; -p)

# ---------- angular average over l in 4 dimensions ----------
def dfact(n):  # (n-1)!! for even n handled by caller
    r = 1
    for i in range(n - 1, 0, -2):
        r *= i
    return r

def ang_avg_monomial(exps):
    """<l1^a1 l2^a2 l3^a3 l4^a4> over S^3, in units of l^(sum a). Zero if any exponent odd."""
    if any(e % 2 for e in exps):
        return sp.Integer(0)
    A = sum(exps) // 2
    num = 1
    for e in exps:
        num *= dfact(e)            # (e-1)!!
    den = 1
    for j in range(A):
        den *= (4 + 2*j)           # 4*6*8*...
    return sp.Rational(num, den)

lsym = sp.Symbol('ell', positive=True)  # |l|

def ang_avg(expr):
    """replace polynomial in l1..l4 by its angular average, as polynomial in ell=|l|"""
    poly = sp.Poly(sp.expand(expr), *L)
    res = sp.Integer(0)
    for monom, coeff in poly.terms():
        c = ang_avg_monomial(monom)
        if c != 0:
            res += coeff*c*lsym**sum(monom)
    return sp.expand(res)

# ---------- radial integrals with hard cutoff ----------
# int d^4l/(2pi)^4 f(|l|) = (2 pi^2)/(2pi)^4 int_0^Lam dl l^3 f = (1/(8 pi^2)) int_0^Lam dl l^3 f
Dl = sp.Symbol('Delta', positive=True)
def radial(power):
    """(1/8pi^2) int_0^Lam dl l^(3+power) / (l^2+Delta)^2, expanded for large Lam: keep Lam^4, Lam^2, log, const"""
    # with t = l^2, T = Lam^2:  (1/2) int_0^T dt t^(1+k)/(t+Delta)^2, k = power/2, expanded for large T,
    # dropping O(1/T).  Checked numerically below against the exact integral.
    T = Lam**2
    lg = sp.log(Lam**2) - sp.log(Dl)
    if power == 0:
        ser = sp.Rational(1, 2)*(lg - 1)
    elif power == 2:
        ser = sp.Rational(1, 2)*(T - 2*Dl*lg + Dl)
    elif power == 4:
        ser = sp.Rational(1, 2)*(T**2/2 - 2*Dl*T + 3*Dl**2*lg - Dl**2)
    else:
        raise ValueError(power)
    # numerical check of the expansion at Lambda=1e3, Delta=0.7
    exact = sp.integrate(lsym**(3 + power)/(lsym**2 + Dl)**2, (lsym, 0, Lam))
    num_ex = sp.N(exact.subs({Lam: 1000, Dl: sp.Rational(7, 10)}), 30)
    num_se = sp.N(ser.subs({Lam: 1000, Dl: sp.Rational(7, 10)}), 30)
    rel = abs((num_ex - num_se)/num_ex)
    say(f'#   radial k={power//2}: expansion vs exact at Lambda=1e3, Delta=0.7: rel.diff {sp.N(rel, 3)} (expect ~1e-6)')
    return sp.expand(ser/(8*sp.pi**2))

say('# radial integrals (1/8pi^2) int_0^Lambda dl l^(3+2k)/(l^2+Delta)^2, large-Lambda:')
R = {}
for pw in (0, 2, 4):
    R[pw] = radial(pw)
    say(f'#   k={pw//2}: {R[pw]}')

# ---------- assemble Pi ----------
say('# assembling Pi_{mu nu rho sigma}(p): 256 components')
Pi = {}
logLam2 = sp.log(Lam**2)
for (i, j, r, s) in itertools.product(range(4), repeat=4):
    if i > j or r > s:
        continue
    num = ang_avg(V1[i, j]*V2[r, s])          # polynomial in ell
    polyn = sp.Poly(num, lsym)
    val = sp.Integer(0)
    for (pw,), c in polyn.terms():
        val += c*R[pw]
    val = -sp.Rational(1, 2)*val              # the -(1/2) in Pi
    val = val.subs(Dl, Delta)
    Pi[(i, j, r, s)] = sp.expand(val)

# integrate over Feynman parameter x the pieces we need: separate by power of Lambda / log
def split(expr):
    """return dict: 'L4','L2','log','fin' coefficients (log means coefficient of log(Lambda^2))"""
    e = sp.expand(expr)
    lgL = sp.log(Lam)                      # sympy writes log(Lambda^2) as 2*log(Lambda)
    d = {}
    d['L4'] = e.coeff(Lam, 4)
    d['L2'] = e.coeff(Lam, 2)
    rest = sp.expand(e - d['L4']*Lam**4 - d['L2']*Lam**2)
    d['log'] = sp.expand(rest.coeff(lgL, 1)/2)   # coefficient of log(Lambda^2)
    d['fin'] = sp.expand(rest - 2*d['log']*lgL)
    assert not d['fin'].has(Lam), 'leftover Lambda dependence'
    return d

say('# x-integration of the log(Lambda^2) coefficient')
Pi_log = {}
for key, val in Pi.items():
    d = split(val)
    Pi_log[key] = sp.expand(sp.integrate(d['log'], (x, 0, 1)))

def Pi_full(i, j, r, s, table):
    a, b = (i, j) if i <= j else (j, i)
    c, d_ = (r, s) if r <= s else (s, r)
    return table[(a, b, c, d_)]

# quadratic form: Gamma_2,log = (1/2) sum h_{ij} Pi_log_{ij rs} h_{rs}  (coefficient of log Lambda^2)
G2log = sp.Integer(0)
for (i, j, r, s) in itertools.product(range(4), repeat=4):
    G2log += sp.Rational(1, 2)*H[i, j]*Pi_full(i, j, r, s, Pi_log)*H[r, s]
G2log = sp.expand(G2log)

# ---------- curvature quadratic forms ----------
h_tr = sum(H[i, i] for i in range(4))
R1 = sp.expand(-sum(P[i]*P[j]*H[i, j] for i in range(4) for j in range(4)) + p2*h_tr)
R1mn = sp.Matrix(4, 4, lambda i, j: sp.Rational(1, 2)*(
    -sum(P[a]*P[j]*H[a, i] for a in range(4)) - sum(P[a]*P[i]*H[a, j] for a in range(4))
    + P[i]*P[j]*h_tr + p2*H[i, j]))
# h(-p) has the same components in this real-symbol bookkeeping (R1 is even in p), so (R1)^2 -> R1*R1
Q_R2 = sp.expand(R1*R1)
Q_Ric = sp.expand(sum(R1mn[i, j]**2 for i in range(4) for j in range(4)))

# ---------- match the p^4 part ----------
# isolate terms of total degree 4 in p (the p^4 log terms) from G2log
def deg_part(expr, deg):
    poly = sp.Poly(expr, *P)
    return sp.expand(sum(c*sp.prod([Pi_**e for Pi_, e in zip(P, mon)]) for mon, c in poly.terms() if sum(mon) == deg))

G4 = deg_part(G2log, 4)
G2 = deg_part(G2log, 2)
G0 = deg_part(G2log, 0)

A, B = sp.symbols('A B')
# G4 should equal (1/16pi^2) * [ A*Q_Ric + B*Q_R2 ]  (then A,B are the 'tHV coefficients)
resid = sp.expand(16*sp.pi**2*G4 - (A*Q_Ric + B*Q_R2))
eqs = sp.Poly(resid, *P, *[H[i, j] for i in range(4) for j in range(i, 4)]).coeffs()
sol = sp.solve(eqs, [A, B], dict=True)
say('# match of the p^4 log(Lambda^2) part onto (R1_{mu nu})^2 and (R1)^2:')
if not sol:
    say('#   NO SOLUTION: p^4 log part is not a combination of Ric^2 and R^2 (check transversality)')
    status = 'FAIL'
else:
    A_, B_ = sp.factor(sol[0][A]), sp.factor(sol[0][B])
    say(f'#   A (Ric^2) = {A_}')
    say(f'#   B (R^2)   = {B_}')
    # the system was overdetermined; confirm the residual vanishes identically
    chk = sp.expand(resid.subs(sol[0]))
    say(f'#   residual identically zero: {chk == 0}')
    say(f'#   -A, -B are the 1974 Delta-L coefficients (sign: see header)')
    say(f'#   at xi=0  : -A={-A_.subs(xi,0)}  -B={-B_.subs(xi,0)}    [1974 eq. 3.40: 1/120, 1/240]')
    say(f'#   at xi=1/6: -A={-A_.subs(xi,sp.Rational(1,6))}  -B={-B_.subs(xi,sp.Rational(1,6))}    [1974 eq. 3.40: 1/120, -1/360]')
    B_expect = -(sp.Rational(1,4)*(xi - sp.Rational(1,6))**2 - sp.Rational(1,360))
    say(f'#   B as a polynomial in xi vs 1974 -[(1/4)(xi-1/6)^2 - 1/360]: difference = {sp.simplify(B_ - B_expect)}')
    ok = (chk == 0 and sp.simplify(A_ + sp.Rational(1, 120)) == 0 and sp.simplify(B_ - B_expect) == 0)
    status = 'OK' if ok else 'MISMATCH'

# ---------- the lower orders, for the record ----------
say('# lower-order log(Lambda^2) terms (bubble only; tadpole not included, so these are partial):')
say(f'#   p^2 part of Gamma_2,log (coefficient of 1/16pi^2): {sp.factor(sp.expand(16*sp.pi**2*G2))}')
say(f'#   p^0 part of Gamma_2,log (coefficient of 1/16pi^2): {sp.factor(sp.expand(16*sp.pi**2*G0))}')
# sqrt(g) to second order: 1 + h/2 + h^2/8 - h_ij h_ij/4 ; the p^0 log term should be (m^4/4) sqrt(g)^(2)
sqrtg2 = sp.expand(h_tr**2/8 - sum(H[i, j]**2 for i in range(4) for j in range(4))/4)
ratio = sp.simplify(sp.expand(16*sp.pi**2*G0)/sqrtg2)
say(f'#   p^0 log term / sqrt(g)^(2) = {ratio}   [vacuum energy: m^4/4 expected]')
# transversality of the p^4 log part: p_mu Pi_log_{mu nu rho sigma} restricted to degree-4 terms
trans = sp.Integer(0)
for (j, r, s) in itertools.product(range(4), repeat=3):
    tot = sum(P[i]*Pi_full(i, j, r, s, Pi_log) for i in range(4))
    trans += deg_part(sp.expand(tot), 5)**2
say(f'#   transversality of p^4 log part (sum of squares of p.Pi, should be 0): {sp.simplify(trans)}')

# Lambda^4 and Lambda^2 pieces of the bubble (x-integrated), for the record
L4tot = sp.Integer(0); L2tot = sp.Integer(0)
for (i, j, r, s) in itertools.product(range(4), repeat=4):
    d = split(Pi_full(i, j, r, s, Pi))
    L4tot += sp.Rational(1, 2)*H[i, j]*sp.integrate(d['L4'], (x, 0, 1))*H[r, s]
    L2tot += sp.Rational(1, 2)*H[i, j]*sp.integrate(d['L2'], (x, 0, 1))*H[r, s]
say(f'#   bubble Lambda^4 term of Gamma_2 (x 16pi^2): {sp.factor(sp.expand(16*sp.pi**2*L4tot))}')
say(f'#   bubble Lambda^2 term of Gamma_2 (x 16pi^2): {sp.factor(sp.expand(16*sp.pi**2*L2tot))}')

say(f'# END ud1962_bubble elapsed={time.time()-t0:.1f}s status={status}')
if out: out.close()
