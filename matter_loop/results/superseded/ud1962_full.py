#!/usr/bin/env python3
"""
Rung 1, technique (a), complete: bubble + tadpole, hard cutoff, all three orders of the 1962 abstract.

For each order (Lambda^4, Lambda^2, log Lambda^2) the sum of the two graphs is fitted to the covariant
structures it should form at second order in h:
    Lambda^4 :  c * sqrt(g)^(2)                                        (cosmological term)
    Lambda^2 :  c1 * m^2 sqrt(g)^(2) + c2 * (sqrt(g) R)^(2)            (cosmological + Newton)
    log      :  p^0: c * m^4 sqrt(g)^(2);  p^2: c * m^2 (sqrt(g) R)^(2);  p^4: A Ric^2 + B R^2
A non-zero residual at Lambda^4 or Lambda^2 is the expected failure of a hard momentum cutoff to respect
general covariance (and the shift-dependence of power-divergent integrals); it is recorded, not hidden.
The log order is regulator-independent and is the one compared with 't Hooft-Veltman 1974 eq. (3.40).

Usage: python3 ud1962_full.py [--log results/ud1962_full.log]
"""
import sys, time, argparse
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

say('# ud1962_full: 1962 cutoff, bubble + tadpole, scalar on flat background, all orders')
say('# radial integrals, large-Lambda forms vs exact:')
U.check_radials(say)

say('# bubble: 256 components, three orders, x-integrated')
Pi = U.bubble_Pi()
Gb = {part: U.quad_form(Pi, part) for part in ('L4', 'L2', 'log')}

say('# tadpole: second-order densities from the metric, then the seagull trace')
Gt, (sqrtg2, A2, sqrtgR2) = U.tadpole_form()
say(f'#   sqrt(g)^(2)        = {sp.factor(sqrtg2)}')
say(f'#   check vs 1/8 h^2 - 1/4 h_ij h_ij: {sp.simplify(sqrtg2 - (U.h_tr**2/8 - sum(U.H[i,j]**2 for i in range(4) for j in range(4))/4)) == 0}')
say(f'#   (sqrt(g) g^mm)^(2) trace = {sp.factor(A2.trace())}')
say(f'#   (sqrt(g) R)^(2)    = {sp.expand(sqrtgR2)}')
# a check on (sqrt g R)^(2): its p^2 structure is even in p and h-bilinear; compare with the known form via
# the identity  (sqrt g R)^(2) = -(1/4)(d h_ab)^2-type terms; we just verify it is transverse-compatible by
# checking that it vanishes for a pure gauge h_ab = p_a e_b + p_b e_a  (linearized diffeomorphism), as the
# total derivative structure requires at this order when combined... (it does not vanish alone; skip).

W = {}
for part in ('L4', 'L2', 'log'):
    W[part] = sp.expand(Gb[part] + Gt[part])

pi2 = 16*sp.pi**2
# ---------------- Lambda^4
say('# --- Lambda^4 ---')
for name, G in (('bubble', Gb['L4']), ('tadpole', Gt['L4']), ('total', W['L4'])):
    sol, res = U.fit(sp.expand(pi2*G), [sqrtg2], ['c'])
    if sol is None:
        say(f'#   {name}: not proportional to sqrt(g)^(2); form (x16pi^2) = {sp.factor(sp.expand(pi2*G))}')
    else:
        say(f'#   {name}: (x16pi^2) = {sol["c"]} * sqrt(g)^(2),  residual {sp.simplify(res)}')

# ---------------- Lambda^2
say('# --- Lambda^2 ---')
for name, G in (('bubble', Gb['L2']), ('tadpole', Gt['L2']), ('total', W['L2'])):
    sol, res = U.fit(sp.expand(pi2*G), [m2s := U.m**2*sqrtg2, sqrtgR2], ['c_m2sqrtg', 'c_sqrtgR'])
    if sol is None:
        # fit the p^0 and p^2 parts separately to see which fails
        G0 = U.deg_part(pi2*G, 0); G2 = U.deg_part(pi2*G, 2)
        s0, r0 = U.fit(G0, [U.m**2*sqrtg2], ['c'])
        s2, r2 = U.fit(G2, [sqrtgR2], ['c'])
        say(f'#   {name}: NOT a combination of m^2 sqrt(g)^(2) and (sqrt(g) R)^(2)')
        say(f'#      p^0 part: {"c=" + str(s0["c"]) if s0 else "no fit"}; p^2 part: {"c=" + str(s2["c"]) if s2 else "no fit"}')
        if not s2:
            say(f'#      p^2 part (x16pi^2) = {sp.factor(G2)}')
    else:
        say(f'#   {name}: (x16pi^2) = {sol["c_m2sqrtg"]} * m^2 sqrt(g)^(2) + {sol["c_sqrtgR"]} * (sqrt(g) R)^(2),  residual {sp.simplify(res)}')

# ---------------- log Lambda^2
say('# --- log Lambda^2 (coefficients of (1/16pi^2) log Lambda^2) ---')
for name, G in (('bubble', Gb['log']), ('tadpole', Gt['log']), ('total', W['log'])):
    G0 = U.deg_part(pi2*G, 0); G2 = U.deg_part(pi2*G, 2); G4 = U.deg_part(pi2*G, 4)
    s0, r0 = U.fit(G0, [U.m**4*sqrtg2], ['c'])
    s2, r2 = U.fit(G2, [U.m**2*sqrtgR2], ['c'])
    s4, r4 = U.fit(G4, [U.Q_Ric, U.Q_R2], ['A', 'B'])
    say(f'#   {name}:')
    say(f'#      p^0: {("= " + str(s0["c"]) + " * m^4 sqrt(g)^(2)") if s0 else "no fit: " + str(sp.factor(G0))}')
    say(f'#      p^2: {("= " + str(s2["c"]) + " * m^2 (sqrt(g) R)^(2)") if s2 else "no fit: " + str(sp.factor(G2))}')
    say(f'#      p^4: {("= " + str(s4["A"]) + " * Ric^2 + " + str(s4["B"]) + " * R^2") if s4 else "no fit"}')

# ---------------- verdicts
say('# --- verdicts ---')
G4 = U.deg_part(pi2*W['log'], 4)
s4, r4 = U.fit(G4, [U.Q_Ric, U.Q_R2], ['A', 'B'])
B_expect = -(sp.Rational(1, 4)*(U.xi - sp.Rational(1, 6))**2 - sp.Rational(1, 360))
ok4 = s4 is not None and sp.simplify(s4['A'] + sp.Rational(1, 120)) == 0 and sp.simplify(s4['B'] - B_expect) == 0 and r4 == 0
say(f'#   p^4 log vs 1974 eq. 3.40 (Delta L coefficients are -A, -B): {"OK" if ok4 else "MISMATCH"}')
G0 = U.deg_part(pi2*W['log'], 0)
s0, r0 = U.fit(G0, [U.m**4*sqrtg2], ['c'])
# vacuum energy  E = 1/2 int_k ln(k^2+m^2):  dE/dm^2 = 1/2 J0 -> log part of E = -(m^4/4)(1/16pi^2) log Lambda^2
say(f'#   p^0 log total = {s0["c"] if s0 else "no fit"} * m^4 sqrt(g)^(2)   [vacuum energy 1/2 Tr ln: -1/4 expected]')
G2 = U.deg_part(pi2*W['log'], 2)
s2, r2 = U.fit(G2, [U.m**2*sqrtgR2], ['c'])
say(f'#   p^2 log total = {s2["c"] if s2 else "no fit"} * m^2 (sqrt(g) R)^(2)   [a_2 = (1/6 - xi) R: should vanish at xi = 1/6]')
s4L2, r4L2 = U.fit(U.deg_part(pi2*W['L2'], 2), [sqrtgR2], ['c'])
say(f'#   Lambda^2 p^2 total covariant (prop. to (sqrt g R)^(2))? {"yes, c=" + str(s4L2["c"]) if s4L2 else "NO - non-covariant remainder (hard cutoff)"}')
sL4, rL4 = U.fit(sp.expand(pi2*W['L4']), [sqrtg2], ['c'])
say(f'#   Lambda^4 total covariant (prop. to sqrt g^(2))? {"yes, c=" + str(sL4["c"]) if sL4 else "NO - non-covariant remainder (hard cutoff)"}')

status = 'OK' if ok4 else 'MISMATCH'
say(f'# END ud1962_full elapsed={time.time()-t0:.1f}s status={status}')
if out: out.close()
