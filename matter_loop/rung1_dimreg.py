#!/usr/bin/env python3
"""
Rung 1, technique (c): dimensional regularization ('t Hooft-Veltman 1972; applied to gravity 1974).

Same two graphs; the loop integrals are continued to n dimensions and only the pole at n = 4 is kept, in
't Hooft-Veltman's convention 1/eps = 1/(8 pi^2 (n-4)) (Ann. IHP A 20 (1974) 69, eq. 3.3). Power divergences
are absent by construction, so only the log order of (a)/(b) has a counterpart here, and the measure term
(a delta(0)) vanishes identically -- the plain and covariant measures must agree.

Expected: Q_pole = epsinv * [ (1/4) m^4 sqrt(g) + (1/2)(xi - 1/6) m^2 sqrt(g) R + (1/120) Ric^2 + ((1/4)(xi-1/6)^2 - 1/360) R^2 ]
i.e. minus the coefficients of (1/16pi^2) log Lambda^2 found with the cutoff and proper time (dictionary:
(1/16pi^2) log Lambda^2 <-> -1/eps), and, since Delta L_tHV = +Gamma_div in our sign, exactly the
coefficients of 't Hooft-Veltman eq. (3.40) for a real scalar.
"""
import time, argparse
import sympy as sp
import ud1962_lib as U

ap = argparse.ArgumentParser()
ap.add_argument('--measure', choices=['plain', 'covariant'], default='plain')
ap.add_argument('--log', default=None)
args = ap.parse_args()
t0 = time.time()
out = open(args.log, 'w') if args.log else None
def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out: out.write(s + '\n'); out.flush()

xi, m = U.xi, U.m
cov = args.measure == 'covariant'
say(f'# rung1 dimreg {args.measure}: pole parts in units of 1/eps = 1/(8 pi^2 (n-4))')
for pw in (0, 2, 4):
    say(f'#   int d^n l (l^2)^{pw//2}/(l^2+Delta)^2 -> {U.dimreg_two_prop(pw)}')
for pw in (0, 2):
    say(f'#   int d^n k (k^2)^{pw//2}/(k^2+m^2)   -> {U.dimreg_one_prop(pw)}')

ME = U.metric_expansion()
sqrtg2, A2, sqrtgR2 = U.second_order_densities(ME)
sigma2 = U.tadpole_symbol(cov, ME)
vertex = U.V1cov if cov else U.V1

say('# bubble')
raw = U.bubble_raw(U.dimreg_two_prop, vertex)
import itertools
Qb = sp.Integer(0)
for (i, j, r, s) in itertools.product(range(4), repeat=4):
    Qb += sp.Rational(1, 2)*U.H[i, j]*U.sym_get(raw, i, j, r, s)*U.H[r, s]
Qb = sp.expand(sp.integrate(sp.expand(Qb), (U.x, 0, 1)))
say('# tadpole')
Qt = U.tadpole_from_symbol(sigma2, U.dimreg_one_prop(0), U.dimreg_one_prop(2))
Q = sp.expand(Qb + Qt)
assert sp.expand(Q).coeff(U.epsinv, 0) == 0  # pure pole
C = sp.expand(Q.coeff(U.epsinv, 1))

say('# --- coefficient of 1/eps ---')
s0, _ = U.fit(U.deg_part(C, 0), [m**4*sqrtg2], ['m4sqrtg2'])
s2, _ = U.fit(U.deg_part(C, 2), [m**2*sqrtgR2], ['m2sqrtgR2'])
s4, r4 = U.fit(U.deg_part(C, 4), [U.Q_Ric, U.Q_R2], ['Ric2', 'R2'])
say(f'#   p^0: {s0}   [expected 1/4]')
say(f'#   p^2: {s2}   [expected (6 xi - 1)/12]')
say(f'#   p^4: {s4}, residual {r4}   [expected 1/120, (60 xi^2 - 20 xi + 1)/240 = 1974 eq. 3.40]')
ok = (s0 == {'m4sqrtg2': sp.Rational(1, 4)} and s2 is not None and sp.simplify(s2['m2sqrtgR2'] - (6*xi - 1)/12) == 0
      and s4 is not None and s4['Ric2'] == sp.Rational(1, 120) and sp.simplify(s4['R2'] - (60*xi**2 - 20*xi + 1)/240) == 0 and r4 == 0)
say(f'#   matches 1974 eq. 3.40 and the log orders of (a), (b): {ok}')
say(f'# END rung1 dimreg {args.measure} elapsed={time.time()-t0:.1f}s status={"OK" if ok else "MISMATCH"}')
if out: out.close()
