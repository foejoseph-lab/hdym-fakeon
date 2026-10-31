#!/usr/bin/env python3
"""
Rung 1, technique (b): proper time (Schwinger 1951 / DeWitt 1965), computer-aided.

Same scalar, same two graphs, same vertices as ud1962_full.py; the propagators are written as Schwinger
integrals, 1/(k^2+m^2) = int_0^inf dt e^{-t(k^2+m^2)}, and the regulator is a lower cutoff t >= 1/Lambda^2 on
proper time. The Gaussian k-integrals are exact; the Feynman parameter reappears as u = x t.

Two outputs:
 1. Gamma_2 split into Lambda^4, Lambda^2, log Lambda^2, fitted to covariant structures as in ud1962_full.py.
    Proper-time regularization respects general covariance, so all three orders should now assemble:
       Lambda^4 : -(1/4) sqrt(g)^(2)                              [from a_0]
       Lambda^2 : (1/2) m^2 sqrt(g)^(2) - (1/2)(1/6 - xi) (sqrt(g) R)^(2)   [from a_2]
       log      : -(1/2) a_4 density                               [unchanged from the cutoff]
    in units (1/16 pi^2), for the quadratic form Q = Gamma_2/(2 Vol) (coefficient of int d^4x F).
 2. The small-t expansion of Tr e^{-tD} at second order in h: tau(t) = -2 t q(t), q the t-integrand of Q,
    tau(t) * 16pi^2 = t^-2 a0^(2) + t^-1 a2^(2) + t^0 a4^(2) + O(t). These are DERIVED here and then compared
    with the covariant heat-kernel coefficients of Vassilevich hep-th/0306138 eqs. 4.26-4.28 expanded to
    second order in h with E = -xi R - m^2, Omega = 0:
       a0 = sqrt(g);  a2 = sqrt(g)(-m^2 + (1/6 - xi) R);
       a4 = sqrt(g)[ (1/2) m^4 + (xi - 1/6) m^2 R + (1/60) Ric^2 + ((60 xi^2 - 20 xi + 1)/120) R^2 ]  mod total
            derivatives and Gauss-Bonnet (Riem^2 -> 4 Ric^2 - R^2 at this order).

Usage: python3 ud1965_propertime.py [--log results/ud1965_propertime.log]
"""
import time, argparse
import sympy as sp
import ud1962_lib as U

t0 = time.time()
ap = argparse.ArgumentParser()
ap.add_argument('--log', default=None)
args = ap.parse_args()
out = open(args.log, 'w') if args.log else None
def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out: out.write(s + '\n'); out.flush()

say('# ud1965_propertime: Schwinger proper time, lower cutoff t >= 1/Lambda^2, bubble + tadpole')
say('# proper-time integrals, large-Lambda forms vs exact incomplete gamma:')
U.check_pt_radials(say)

pi2 = 16*sp.pi**2
xi, m = U.xi, U.m

say('# densities from the metric (second order in h)')
dens = U.second_order_densities()
sqrtg2, A2, sqrtgR2 = dens

# ---------------- 1. the regulated effective action
say('# bubble, proper time')
Pi = U.bubble_Pi(U.pt_two_prop)
Gb = {part: U.quad_form(Pi, part) for part in ('L4', 'L2', 'log')}
say('# tadpole, proper time')
Gt_raw, _ = U.tadpole_raw(U.pt_one_prop(0), U.pt_one_prop(2), dens)
Gt = U.split(Gt_raw)
W = {part: sp.expand(Gb[part] + Gt[part]) for part in ('L4', 'L2', 'log')}

say('# --- Lambda^4 (x16pi^2) ---')
for name, G in (('bubble', Gb['L4']), ('tadpole', Gt['L4']), ('total', W['L4'])):
    sol, res = U.fit(sp.expand(pi2*G), [sqrtg2], ['c'])
    say(f'#   {name}: ' + (f'{sol["c"]} * sqrt(g)^(2), residual {res}' if sol else f'NOT prop. to sqrt(g)^(2): {sp.factor(sp.expand(pi2*G))}'))
say('#   expected total: -1/4 * sqrt(g)^(2)')

say('# --- Lambda^2 (x16pi^2) ---')
for name, G in (('bubble', Gb['L2']), ('tadpole', Gt['L2']), ('total', W['L2'])):
    sol, res = U.fit(sp.expand(pi2*G), [m**2*sqrtg2, sqrtgR2], ['c_m2sqrtg', 'c_sqrtgR'])
    say(f'#   {name}: ' + (f'{sol["c_m2sqrtg"]} * m^2 sqrt(g)^(2) + {sol["c_sqrtgR"]} * (sqrt(g) R)^(2), residual {res}' if sol else 'NOT a combination'))
say('#   expected total: 1/2 * m^2 sqrt(g)^(2) + -(1/6 - xi)/2 * (sqrt(g) R)^(2)')

say('# --- log Lambda^2 (x16pi^2) ---')
G = W['log']
G0 = U.deg_part(pi2*G, 0); G2 = U.deg_part(pi2*G, 2); G4 = U.deg_part(pi2*G, 4)
s0, _ = U.fit(G0, [m**4*sqrtg2], ['c']); s2, _ = U.fit(G2, [m**2*sqrtgR2], ['c']); s4, r4 = U.fit(G4, [U.Q_Ric, U.Q_R2], ['A', 'B'])
say(f'#   total p^0: {s0["c"] if s0 else "no fit"} * m^4 sqrt(g)^(2)        [expected -1/4]')
say(f'#   total p^2: {s2["c"] if s2 else "no fit"} * m^2 (sqrt(g) R)^(2)    [expected -(6 xi - 1)/12]')
say(f'#   total p^4: {s4["A"] if s4 else "no fit"} * Ric^2 + {s4["B"] if s4 else ""} * R^2   [expected -1/120, -(60 xi^2 - 20 xi + 1)/240]')
# constants (gamma_E etc.) differ between regulators; record the finite p^4 constant for the record
ok_div = all([
    U.fit(sp.expand(pi2*W['L4']), [sqrtg2], ['c'])[0] == {'c': -sp.Rational(1, 4)},
    U.fit(sp.expand(pi2*W['L2']), [m**2*sqrtg2, sqrtgR2], ['c_m2sqrtg', 'c_sqrtgR'])[0] == {'c_m2sqrtg': sp.Rational(1, 2), 'c_sqrtgR': sp.factor(-(sp.Rational(1, 6) - xi)/2)},
    s0 == {'c': -sp.Rational(1, 4)},
    s2 is not None and sp.simplify(s2['c'] + (6*xi - 1)/12) == 0,
    s4 is not None and sp.simplify(s4['A'] + sp.Rational(1, 120)) == 0 and sp.simplify(s4['B'] + (60*xi**2 - 20*xi + 1)/240) == 0 and r4 == 0,
])
say(f'#   all three orders covariant and as expected: {ok_div}')

# ---------------- 2. heat-kernel coefficients derived
say('# --- Tr e^{-tD} at second order in h: small-t expansion, derived ---')
t = U.tsym
raw_b = U.bubble_raw(U.pt_two_prop_integrand)
qb = sp.Integer(0)
import itertools
for (i, j, r, s) in itertools.product(range(4), repeat=4):
    qb += sp.Rational(1, 2)*U.H[i, j]*U.sym_get(raw_b, i, j, r, s)*U.H[r, s]
qt, _ = U.tadpole_raw(U.pt_one_prop_integrand(0), U.pt_one_prop_integrand(2), dens)
tau = sp.expand(-2*t*(qb + qt))                      # tau(t) = T_2(t)/(2 Vol), still with x inside
# series in t to O(t^1): write e^{-t Delta} expanded; Delta = m^2 + x(1-x)p^2 is inside raw_b via substitution
tau_ser = sp.series(tau, t, 0, 1).removeO()
tau_ser = sp.expand(tau_ser)
coef = {}
for k, pw in ((0, -2), (2, -1), (4, 0)):
    c = tau_ser.coeff(t, pw)
    c = sp.expand(sp.integrate(sp.expand(c), (U.x, 0, 1)))
    coef[k] = sp.expand(pi2*c)
    say(f'#   a_{k}^(2) derived (x16pi^2): {len(sp.Add.make_args(coef[k]))} terms')

# expected covariant forms at second order
exp0 = sqrtg2
exp2 = sp.expand(-m**2*sqrtg2 + (sp.Rational(1, 6) - xi)*sqrtgR2)
exp4 = sp.expand(sp.Rational(1, 2)*m**4*sqrtg2 + (xi - sp.Rational(1, 6))*m**2*sqrtgR2
                 + sp.Rational(1, 60)*U.Q_Ric + (60*xi**2 - 20*xi + 1)/sp.Integer(120)*U.Q_R2)
ok_hk = True
for k, ex, name in ((0, exp0, 'sqrt(g)'), (2, exp2, 'sqrt(g)(-m^2 + (1/6-xi)R)'), (4, exp4, 'sqrt(g)[m^4/2 + (xi-1/6)m^2 R + Ric^2/60 + (60xi^2-20xi+1)R^2/120]')):
    diff = sp.simplify(sp.expand(coef[k] - ex))
    say(f'#   a_{k}: derived - Vassilevich 4.2{6 + k//2} [{name}] = {diff}')
    ok_hk = ok_hk and diff == 0
# the derived a_4 fitted on its own, as a check that the basis is complete
G4h = U.deg_part(coef[4], 4)
s4h, r4h = U.fit(G4h, [U.Q_Ric, U.Q_R2], ['A', 'B'])
say(f'#   a_4 p^4 part fitted: {s4h["A"] if s4h else "no fit"} Ric^2 + {s4h["B"] if s4h else ""} R^2, residual {r4h}')

status = 'OK' if (ok_div and ok_hk) else 'MISMATCH'
say(f'# END ud1965_propertime elapsed={time.time()-t0:.1f}s status={status}')
if out: out.close()
