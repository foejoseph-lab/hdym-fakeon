#!/usr/bin/env python3
"""
Rung 1 driver: one real scalar on g = delta + h, one loop, second order in h; bubble + tadpole.

  --regulator cutoff      1962: hard momentum cutoff on the (Feynman-parametrised) loop momentum
  --placement shifted     (cutoff only; default) |l| < Lambda on the Feynman-shifted momentum l = k + x p
  --placement unshifted   (cutoff only) |k| < Lambda on the propagator momentum k itself, no Feynman parameter;
                          large-k expansion, exact in Lambda^4, Lambda^2, log Lambda^2 (added 1 Nov 2026)
  --regulator propertime  1965: Schwinger proper time with lower cutoff t >= 1/Lambda^2
  --measure plain         operator K = -d(sqrt g g^{mu nu} d) + sqrt g (m^2 + xi R)   [phi as the integration variable]
  --measure covariant     operator g^(-1/4) K g^(-1/4)   [chi = g^(1/4) phi: DeWitt's covariant measure]
  --heatkernel            (propertime only) derive a_0, a_2, a_4 at second order in h from the small-t expansion
                          of Tr e^{-tD} and compare with Vassilevich hep-th/0306138 eqs. 4.26-4.28.

Each order (Lambda^4, Lambda^2, log Lambda^2) of Q = Gamma_2/(2 Vol) is fitted to the covariant structures:
  Lambda^4 : c sqrt(g)^(2);   Lambda^2 : c1 m^2 sqrt(g)^(2) + c2 (sqrt g R)^(2);
  log      : p^0: c m^4 sqrt(g)^(2);  p^2: c m^2 (sqrt g R)^(2);  p^4: A Ric^2 + B R^2.
Expectations (units 1/16pi^2), from W = 1/2 Tr ln D = -1/2 int dt/t Tr e^{-tD} with the heat-kernel coefficients:
  Lambda^4 -> -1/4 sqrt(g);  Lambda^2 -> -(1/2) a_2 = 1/2 m^2 sqrt g - (1/2)(1/6 - xi) sqrt g R;
  log -> -(1/2) a_4 = -(1/4) m^4 - (1/2)(xi - 1/6) m^2 R - (1/120) Ric^2 - ((60 xi^2 - 20 xi + 1)/240) R^2;
and the p^4 log coefficients are 't Hooft-Veltman 1974 eq. (3.40) with Delta L = -(A Ric^2 + B R^2).

Prediction for --placement unshifted (written 1 Nov 2026 before the run; surface-term formula for shifting a
cut-off integral, int_{|k|<L} [F(k+a) - F(k)] = sum_n (1/n!) a^n . oint n dS d^(n-1)F, with F ~ k^0 at large k):
  Lambda^4 : IDENTICAL to the shifted placement (the order-a surface term is odd and vanishes), hence still
             not proportional to sqrt(g)^(2) -- the hard cutoff fails covariance at Lambda^4 with either placement;
  Lambda^2 : DIFFERS from the shifted placement by a non-covariant term quadratic in p (order a^2 surface term);
  log      : IDENTICAL (the first surface term that could touch it is finite);
so the expected status is LOG-OK_POWER-ORDERS-NONCOVARIANT and the log line equals eq. (3.40) as before.
"""
import time, argparse, itertools
import sympy as sp
import ud1962_lib as U

ap = argparse.ArgumentParser()
ap.add_argument('--regulator', choices=['cutoff', 'propertime'], required=True)
ap.add_argument('--measure', choices=['plain', 'covariant'], default='plain')
ap.add_argument('--placement', choices=['shifted', 'unshifted'], default='shifted')
ap.add_argument('--heatkernel', action='store_true')
ap.add_argument('--log', default=None)
args = ap.parse_args()
t0 = time.time()
out = open(args.log, 'w') if args.log else None
def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out: out.write(s + '\n'); out.flush()

unshifted = args.regulator == 'cutoff' and args.placement == 'unshifted'
name = f'rung1 {args.regulator} {args.measure}' + (' unshifted' if unshifted else '')
say(f'# {name}')
pi2 = 16*sp.pi**2
xi, m = U.xi, U.m
cov = args.measure == 'covariant'

if args.regulator == 'cutoff':
    two, one = U.radial_two_prop, U.radial_one_prop
    U.check_radials(say)
    if unshifted:
        U.check_radial_hard(say)
else:
    two, one = U.pt_two_prop, U.pt_one_prop
    U.check_pt_radials(say)

say('# metric expansion and densities')
ME = U.metric_expansion()
sqrtg2, A2, sqrtgR2 = U.second_order_densities(ME)
say('# tadpole symbol' + (' (covariant measure)' if cov else ''))
sigma2 = U.tadpole_symbol(cov, ME)

vertex = U.V1cov if cov else U.V1
say('# bubble' + (' (cutoff on the unshifted propagator momentum, large-k expansion)' if unshifted else ''))
if unshifted:
    raw = U.bubble_raw_unshifted(vertex)
    Pi = {key: U.split(val) for key, val in raw.items()}
else:
    Pi = U.bubble_Pi(two, vertex)
Gb = {part: U.quad_form(Pi, part) for part in ('L4', 'L2', 'log')}
if unshifted:
    say('# bubble again with the shifted placement, for the difference')
    Pi_sh = U.bubble_Pi(two, vertex)
    Gb_sh = {part: U.quad_form(Pi_sh, part) for part in ('L4', 'L2', 'log')}
    say('# --- unshifted minus shifted placement, bubble (x16pi^2) ---')
    for part in ('L4', 'L2', 'log'):
        d = sp.expand(16*sp.pi**2*(Gb[part] - Gb_sh[part]))
        if d == 0:
            say(f'#   {part}: identical')
        else:
            # the four scalars quadratic in h and in p: fit the difference onto them
            htr = U.h_tr
            b1 = U.p2*htr**2
            b2 = U.p2*sum(U.H[i, j]**2 for i in range(4) for j in range(4))
            b3 = htr*sum(U.P[i]*U.P[j]*U.H[i, j] for i in range(4) for j in range(4))
            b4 = sum(U.P[i]*U.P[j]*U.H[i, k]*U.H[j, k] for i in range(4) for j in range(4) for k in range(4))
            fd, rd = U.fit(d, [b1, b2, b3, b4], ['p2_h2', 'p2_hh', 'h_pph', 'phhp'])
            say(f'#   {part}: differs; fit onto {{p^2 h^2, p^2 h_ij h_ij, h p_i p_j h_ij, p_i h_ik h_jk p_j}} = {fd}, residual {rd}')
say('# tadpole')
Gt = U.split(U.tadpole_from_symbol(sigma2, one(0), one(2)))
W = {part: sp.expand(Gb[part] + Gt[part]) for part in ('L4', 'L2', 'log')}

def report(label, G, basis, names, expect=None):
    sol, res = U.fit(sp.expand(pi2*G), basis, names)
    if sol is None:
        say(f'#   {label}: NOT a combination of {names}; p^0 part = {sp.factor(U.deg_part(pi2*G, 0))}')
        return None
    txt = ' + '.join(f'({sol[n]}) {n}' for n in names)
    say(f'#   {label}: {txt}' + (f'   [expected {expect}]' if expect else ''))
    return sol

ok = True
say('# --- Lambda^4 (x16pi^2) ---')
for lab, G in (('bubble', Gb['L4']), ('tadpole', Gt['L4'])):
    report(lab, G, [sqrtg2], ['sqrtg2'])
s = report('total', W['L4'], [sqrtg2], ['sqrtg2'], '-1/4')
ok &= s is not None and s['sqrtg2'] == -sp.Rational(1, 4)

say('# --- Lambda^2 (x16pi^2) ---')
for lab, G in (('bubble', Gb['L2']), ('tadpole', Gt['L2'])):
    report(lab, G, [m**2*sqrtg2, sqrtgR2], ['m2sqrtg2', 'sqrtgR2'])
s = report('total', W['L2'], [m**2*sqrtg2, sqrtgR2], ['m2sqrtg2', 'sqrtgR2'], '1/2, -(1/6 - xi)/2')
ok &= s is not None and s['m2sqrtg2'] == sp.Rational(1, 2) and sp.simplify(s['sqrtgR2'] + (sp.Rational(1, 6) - xi)/2) == 0

say('# --- log Lambda^2 (x16pi^2) ---')
G = W['log']
s0 = report('total p^0', U.deg_part(G, 0), [m**4*sqrtg2], ['m4sqrtg2'], '-1/4')
s2 = report('total p^2', U.deg_part(G, 2), [m**2*sqrtgR2], ['m2sqrtgR2'], '-(6 xi - 1)/12')
s4 = report('total p^4', U.deg_part(G, 4), [U.Q_Ric, U.Q_R2], ['Ric2', 'R2'], '-1/120, -(60 xi^2 - 20 xi + 1)/240')
ok &= s0 is not None and s0['m4sqrtg2'] == -sp.Rational(1, 4)
ok &= s2 is not None and sp.simplify(s2['m2sqrtgR2'] + (6*xi - 1)/12) == 0
ok &= s4 is not None and sp.simplify(s4['Ric2'] + sp.Rational(1, 120)) == 0 and sp.simplify(s4['R2'] + (60*xi**2 - 20*xi + 1)/240) == 0
say(f'#   all orders covariant and as expected: {bool(ok)}')
log_ok = (s0 is not None and s0['m4sqrtg2'] == -sp.Rational(1, 4)
          and s2 is not None and sp.simplify(s2['m2sqrtgR2'] + (6*xi - 1)/12) == 0
          and s4 is not None and sp.simplify(s4['Ric2'] + sp.Rational(1, 120)) == 0 and sp.simplify(s4['R2'] + (60*xi**2 - 20*xi + 1)/240) == 0)

if args.heatkernel and args.regulator == 'propertime':
    say('# --- heat-kernel coefficients at second order in h, derived from Tr e^{-tD} ---')
    t = U.tsym
    # truncated exponentials: tau = -2 t q(t); need powers t^-2, t^-1, t^0
    def two_trunc(power):
        n = power // 2
        Dl = U.Dl
        ex = sum((-t*Dl)**j/sp.factorial(j) for j in range(0, 4))
        return sp.expand(sp.factorial(n + 1)*ex/(16*sp.pi**2*t**(n + 1)))
    def one_trunc(power):
        n = power // 2
        ex = sum((-t*m**2)**j/sp.factorial(j) for j in range(0, 4))
        return sp.expand(sp.factorial(n + 1)*ex/(16*sp.pi**2*t**(n + 2)))
    raw_b = U.bubble_raw(two_trunc, vertex)
    qb = sp.Integer(0)
    for (i, j, r, s_) in itertools.product(range(4), repeat=4):
        qb += sp.Rational(1, 2)*U.H[i, j]*U.sym_get(raw_b, i, j, r, s_)*U.H[r, s_]
    qt = U.tadpole_from_symbol(sigma2, one_trunc(0), one_trunc(2))
    tau = sp.expand(-2*t*(qb + qt))
    coef = {}
    for k, pw in ((0, -2), (2, -1), (4, 0)):
        c = tau.coeff(t, pw)
        c = sp.expand(sp.integrate(sp.expand(c), (U.x, 0, 1)))
        coef[k] = sp.expand(pi2*c)
    exp0 = sqrtg2
    exp2 = sp.expand(-m**2*sqrtg2 + (sp.Rational(1, 6) - xi)*sqrtgR2)
    exp4 = sp.expand(sp.Rational(1, 2)*m**4*sqrtg2 + (xi - sp.Rational(1, 6))*m**2*sqrtgR2
                     + sp.Rational(1, 60)*U.Q_Ric + (60*xi**2 - 20*xi + 1)/sp.Integer(120)*U.Q_R2)
    okhk = True
    for k, ex, lab in ((0, exp0, 'sqrt g'), (2, exp2, 'sqrt g (-m^2 + (1/6 - xi) R)'),
                       (4, exp4, 'sqrt g [m^4/2 + (xi - 1/6) m^2 R + Ric^2/60 + (60xi^2 - 20xi + 1) R^2/120]')):
        diff = sp.expand(coef[k] - ex)
        if diff == 0:
            say(f'#   a_{k}^(2) derived == Vassilevich eq. 4.2{6 + k//2} [{lab}]')
        else:
            okhk = False
            s_, r_ = U.fit(coef[k], [sqrtg2] if k == 0 else [m**2*sqrtg2, sqrtgR2] if k == 2 else [m**4*sqrtg2, m**2*sqrtgR2, U.Q_Ric, U.Q_R2],
                            ['c'] if k == 0 else ['c1', 'c2'] if k == 2 else ['c1', 'c2', 'A', 'B'])
            say(f'#   a_{k}^(2) derived != table; fit: {s_}; difference from table: {sp.factor(diff) if len(sp.Add.make_args(diff)) < 12 else str(len(sp.Add.make_args(diff))) + " terms"}')
    say(f'#   heat-kernel coefficients reproduced: {okhk}')
    ok &= okhk

status = 'OK' if ok else ('LOG-OK_POWER-ORDERS-NONCOVARIANT' if log_ok else 'MISMATCH')
say(f'# END {name} elapsed={time.time()-t0:.1f}s status={status}')
if out: out.close()
