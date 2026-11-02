#!/usr/bin/env python3
"""
2 Nov 2026, the check that located the k_perp^2 slip (results note section 9).  Run from matter_loop/.  Compares, at p = 0 in dim reg
with the plain measure, the code's bubble and tadpole finite parts with an evaluation EXACT IN n of the same two pieces:
    tadpole = int [1/2 sqrt(g)^(2) K_0 - (h/4) h_{mu nu} k_mu k_nu + 1/2 (hh)_{mu nu} k_mu k_nu]/K_0,
    bubble  = -1/4 int [(h/2) K_0 - h_{mu nu} k_mu k_nu]^2 / K_0^2,     K_0 = k^2 + m^2 (n-dimensional k^2),
with <k_mu k_nu> = delta k^2/n, <k k k k> = (delta delta + ..) k^4/(n(n+2)), and the total shown to be -(1/2) Gamma(-n/2) m^n (4 pi)^{-n/2}
sqrt(g)^(2) = sqrt(g)^(2) x the flat vacuum energy, i.e. the MS-bar constant -3/8 (units 1/16 pi^2).  At the time of the check the
code's tadpole came out 1/(512 pi^2) h^2 - 1/(256 pi^2) h_ij h_ij short of the exact one: the k_perp^2 of the flat trace term.
Kept with the first (wrong) log, rung1_massive_plain_sandbox_4d_tadpole_symbol.log.  The call U.tadpole_dimreg_full(sigma2) below is
the pre-fix signature; after the fix the code agrees with the exact evaluation.
"""
import sys; sys.path.insert(0, ".."); import sympy as sp, itertools, ud1962_lib as U
m, xi = U.m, U.xi
ME = U.metric_expansion(); sqrtg2, A2, sqrtgR2 = U.second_order_densities(ME)
sigma2 = U.tadpole_symbol(False, ME)
Qt = U.tadpole_dimreg_full(sigma2)   # pre-fix call; now defaults to the plain measure WITH the k_perp^2 term
Qt0 = U.deg_part(Qt, 0)
print('tadpole p0:', sp.factor(sp.expand(Qt0.coeff(U.epsinv,1))), '| finite:', U.fit(sp.expand(Qt0 - Qt0.coeff(U.epsinv,1)*U.epsinv), [m**4*U.h_tr**2, m**4*sum(U.H[i,j]**2 for i in range(4) for j in range(4))], ['h2','hh']))
tab = U.bubble_massive('dimreg', False, order=0)
Qb = sum(sp.Rational(1,2)*U.H[i,j]*U.sym_get(tab,i,j,r,s)*U.H[r,s] for (i,j,r,s) in itertools.product(range(4), repeat=4))
Qb = sp.expand(Qb.subs(U.Ssym, U.p2)); Qb0 = U.deg_part(Qb, 0)
print('bubble p0 pole:', U.fit(sp.expand(Qb0.coeff(U.epsinv,1)), [m**4*U.h_tr**2, m**4*sum(U.H[i,j]**2 for i in range(4) for j in range(4))], ['h2','hh']))
print('bubble p0 fin:', U.fit(sp.expand(Qb0 - Qb0.coeff(U.epsinv,1)*U.epsinv), [m**4*U.h_tr**2, m**4*sum(U.H[i,j]**2 for i in range(4) for j in range(4))], ['h2','hh']))
# analytic, exact in n: tadpole = int [1/2 sqrtg2 K0 - (h/4) hkk + 1/2 hhkk]/K0 ; bubble = -1/4 int K1^2/K0^2
n = 4 + U.epsd
X = U._eps_split(U.dimreg_one_prop_full(2))  # int k^2/K0 = (A,B)
J0 = U._eps_split(U.dimreg_one_prop_full(0))
K4 = U._eps_split(U.dimreg_two_prop_full(4)); K4 = (K4[0].subs(U.Dl, m**2), K4[1].subs(U.Dl, m**2))
K2 = U._eps_split(U.dimreg_two_prop_full(2)); K2 = (K2[0].subs(U.Dl, m**2), K2[1].subs(U.Dl, m**2))
K0 = U._eps_split(U.dimreg_two_prop_full(0)); K0 = (K0[0].subs(U.Dl, m**2), K0[1].subs(U.Dl, m**2))
h2 = U.h_tr**2; hh = sum(U.H[i,j]**2 for i in range(4) for j in range(4))
def comb(c, pair):  # c rational in epsd times (A/eps + B) -> finite part
    c0 = c.subs(U.epsd,0); c1 = sp.diff(c, U.epsd).subs(U.epsd,0)
    return sp.expand(c0*pair[1] + c1*pair[0]), sp.expand(c0*pair[0])
# tadpole finite: 1/2 sqrtg2 (J2 + m^2 J0) - (h/4) h_mumu <k^2>/n ... careful: hkk -> h_tr k^2/n ; hhkk -> hh k^2/n
t_fin = sp.expand(sp.Rational(1,2)*sqrtg2*(comb(sp.Integer(1), X)[0] + m**2*comb(sp.Integer(1), J0)[0]) - U.h_tr*U.h_tr/4*comb(1/n, X)[0] + sp.Rational(1,2)*hh*comb(1/n, X)[0])
print('analytic tadpole finite:', U.fit(t_fin, [m**4*h2, m**4*hh], ['h2','hh']))
# bubble: -1/4 int [(h/2)K0 - hkk]^2/K0^2 = -1/4 [ h^2/4 int 1 - h hkk... ] : -1/4 [ (h^2/4)(int 1) - (h/1) h_tr <k^2>/n int 1/K0 ... ] do exactly:
# (h/2 K0 - hkk)^2/K0^2 = h^2/4 - h hkk/K0 + (hkk)^2/K0^2 ; <hkk> = h_tr k^2/n ; <(hkk)^2> = (h2+2hh) k^4/(n(n+2))
b_fin = sp.expand(-sp.Rational(1,4)*(-U.h_tr*U.h_tr*comb(1/n, X)[0] + (h2+2*hh)*comb(1/(n*(n+2)), K4)[0]))
print('analytic bubble finite:', U.fit(b_fin, [m**4*h2, m**4*hh], ['h2','hh']))
tot = sp.expand(16*sp.pi**2*(t_fin + b_fin))
print('analytic total finite (x16pi^2):', U.fit(tot, [m**4*sqrtg2], ['sqrtg2']))
