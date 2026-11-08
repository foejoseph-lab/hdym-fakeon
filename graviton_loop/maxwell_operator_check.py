#!/usr/bin/env python3
"""Validation of graviton_operator_derive.py's machinery on Maxwell in Feynman gauge: L = sqrt(g)[F^2/4 + (D.A)^2/2].
Expected: minimal with s = -1 (EL = -D^2 A ...), E = -Ric exactly, connection term on constant vectors +(1/3) Ric.  ~10 s.
Run:  python3 maxwell_operator_check.py > results/maxwell_operator_check.log
"""
import sympy as sp, itertools, random, sys, time
sys.argv=['x']
import os
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'graviton_operator_derive.py')).read().split('def main():')[0]
g={}; exec(src,g); globals().update({k:v for k,v in g.items() if not k.startswith('__')})
rng=random.Random(20261108); R=riemann_sample(rng); Ric,Rs,Ric2,Riem2=invariants(R)
geo=Geometry(R)
Gb,gbinv,sqrtgb=geo.Gb,geo.gbinv,geo.sqrtgb
def L_maxwell(A, gauge_sign):
    # A: list of 4 ring elements (A_mu lower); L = sqrt(g)[ 1/4 F_{mu nu} F^{mu nu} + gauge_sign * 1/2 (D.A)^2 ]
    T=trunc
    dA=[[T(A[n].diff(XR[m]),3) for n in idx] for m in idx]       # d_m A_n
    F=[[dA[m][n]-dA[n][m] for n in idx] for m in idx]
    F2=T(sum((gbinv[m][a]*gbinv[n][b]*F[m][n]*F[a][b] for m in idx for n in idx for a in idx for b in idx),RING(0)),2)
    DA=[[T(dA[m][n]-sum((Gb[l][m][n]*A[l] for l in idx),RING(0)),3) for n in idx] for m in idx]   # D_m A_n
    divA=T(sum((gbinv[m][n]*DA[m][n] for m in idx for n in idx),RING(0)),3)
    return T(sqrtgb*(F2/4 + gauge_sign*divA**2/2),2)
def probes_vec():
    P=[]
    for a in idx:
        E=[RING(1 if i==a else 0) for i in idx]
        P.append(('0',(a,),E))
        for m in idx: P.append(('1',(a,m),[E[i]*XR[m] for i in idx]))
        for m in idx:
            for n in idx:
                if m<=n:
                    w=QQ(1,2) if m==n else QQ(1)
                    P.append(('2',(a,m,n),[E[i]*XR[m]*XR[n]*w for i in idx]))
    return P
def el_vec(base, gauge_sign, probes):
    G={}
    for kind,key,p in probes:
        A=[base[i]+DL*p[i] for i in idx]
        G[(kind,key)]=coeff_dl(L_maxwell(A,gauge_sign))
    EL=[]
    for a in idx:
        G0=G[('0',(a,))]; G1=[G[('1',(a,m))] for m in idx]; G2={(m,n):G[('2',(a,m,n))] for m in idx for n in idx if m<=n}
        expr=G0
        for m in idx: expr-=(-XR[m]*G0+G1[m]).diff(XR[m])
        for m in idx:
            for n in idx:
                mm,nn=min(m,n),max(m,n); w=QQ(1) if m==n else QQ(1,2)
                F2=XR[m]*XR[n]*G0/2-(XR[n]*G1[m]+XR[m]*G1[n])/2+w*G2[(mm,nn)]
                expr+=F2.diff(XR[m]).diff(XR[n])
        EL.append(to_sympy(at0(expr)))
    return EL
probes=probes_vec()
t=time.time()
# kinetic probe
base=[XR[0]**2/2 if i==1 else RING(0) for i in idx]
print("EL[x0^2/2 e_1] =", el_vec(base,1,probes), " (minimal: s e_1)")
base=[XR[0]*XR[2] if i==1 else RING(0) for i in idx]
print("EL[x0 x2 e_1] =", el_vec(base,1,probes), " (want 0)")
# constant fields -> O e_a = s (D^2 e_a + E e_a); D^2 on constant vector at 0: -(d_m Gamma^c_{m a}) e_c
dGam0={(m,c,a): to_sympy(at0(Gb[c][m][a].diff(XR[m]))) for m in idx for c in idx for a in idx}
O=sp.zeros(4,4)
for a in idx:
    base=[RING(1 if i==a else 0) for i in idx]
    col=el_vec(base,1,probes)
    for i in idx: O[i,a]=col[i]
D2c=sp.Matrix(4,4,lambda c,a: -sum(dGam0[(m,c,a)] for m in idx))
RicM=sp.Matrix(4,4,lambda a,b: Ric[(a,b)])
print("time", round(time.time()-t), "s")
print("D2 on constants:", D2c.tolist(), " vs (1/3)Ric:", (D2c - RicM/3)==sp.zeros(4,4), " vs -(1/3)Ric:", (D2c + RicM/3)==sp.zeros(4,4))
for s in (1,2,-1,-2):
    E=(O/s - D2c)
    print(f"s={s}: E - Ric = {(E-RicM).applyfunc(sp.simplify).tolist()}  E + Ric = {(E+RicM).applyfunc(sp.simplify).tolist()}")
