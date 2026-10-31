"""
The input check removed from rung1_dewitt.py before its first logged run (31 Oct 2026), reconstructed from the
session; the debug run's output was not saved.  It asked that nabla_gamma nabla_delta R_{abcd}, symmetrised in
the derivative pair (gamma, delta) only, vanish at the origin for every component.  That object cannot vanish
for generic curvature (differentiated second Bianchi identity + Ricci identity), so the check failed on random
samples and passed on S^4 and S^2 x S^2.  It was replaced by computing Box R(0).  Kept per the repository rule
that superseded code stays in results/superseded/.  Context: Rf = Riemann of the ansatz metric as polynomials
through degree 2, R = input Riemann at the origin, G = Christoffels, at0 = evaluate at the origin.
"""
# dG0 = {(l, a, b, c): at0(sp.diff(G[l][a][b], X[c])) for l in idx for a in idx for b in idx for c in idx}
# def nabla2(k, gam, dl):
#     a, b, c, d = k
#     v = sp.diff(Rf[k], X[gam], X[dl])
#     v -= sum(dG0[(l, dl, a, gam)] * R[(l, b, c, d)] + dG0[(l, dl, b, gam)] * R[(a, l, c, d)]
#              + dG0[(l, dl, c, gam)] * R[(a, b, l, d)] + dG0[(l, dl, d, gam)] * R[(a, b, c, l)] for l in idx)
#     return at0(v)
# ok = True
# for k in Rf:
#     for gam in idx:
#         for dl in idx:
#             if gam <= dl and sp.expand(nabla2(k, gam, dl) + nabla2(k, dl, gam)) != 0:
#                 ok = False
# rep['sym_nabla2_riemann_at_0'] = ok
