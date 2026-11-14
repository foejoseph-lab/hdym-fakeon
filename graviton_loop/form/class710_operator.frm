* class710_operator.frm
* Rung 3, step 3b, second method for route 2: FORM.  STAGE F1 -- THE CLASS (7.10): K(phi), c A SYMBOL.
* PREDICTIONS WRITTEN 13 Nov 2026 17:52 ET (clock read by the command that wrote them), before any F1 code existed and
* BEFORE any log of the author's R8 or main runs of class710_operator_derive.py had been seen.  The author held his yes until
* those runs were in (pushed 19:31 ET, all three status=PASS) and gave it at 19:36 ET.  They stand as first written (rule 11).
* Working names: graviton_loop/form/class710_operator.frm and graviton_loop/form_class710_check.py.
* F0 stands behind it: the author's run of 13 Nov 17:47 ET, status=PASS, the same program hash as the sandbox's (c937891be4b0).
*
* WHAT CHANGES FROM F0.
*   The Lagrangian is sqrt(g)(R - 1/2 K(Phi) g^{ab} d Phi d Phi), Phi = phi-bg + x p, with K(Phi) = K + x p K1 + 1/2 x^2 p^2 K2;
*   K, K1, K2 are the function and its derivatives at phi-bg(x), so D_m K = K1 v_m, D_m K1 = K2 v_m (a third derivative must
*   never appear).  The gauge condition is C_m = D^n h_mn - 1/2 D_m h - c G p v_m with G = K (the image family; G = 1, the
*   literal gauge, by a switch, reported and not predicted).
*   FORM writes, besides the two equations and the ghost, the covariant derivative D_n B^m of the first-derivative coefficient
*   (free indices n, m), which the conversion to omega, E, Omega needs.  The reader does the 11 x 11 matrix algebra on random
*   backgrounds:  S from the second-derivative part,  omega^m = 1/2 S^-1 B^m,  E = S^-1 C - D_m omega^m - omega.omega,
*   Omega_mn = D_m omega_n - D_n omega_m + [omega_m, omega_n] + the metric's curvature on h,  then (3.40) by
*   hv1974_gravscalar_M.hv340, the ghost, and the R8 fit -- the same downstream functions as class710_operator_derive.py.
*   Typed into the reader, not derived: D_m S = K1 v_m on the scalar entry (immediate once F1.1 holds), and the curvature of the
*   metric connection on symmetric tensors.
*
* PREDICTIONS (rule 25).  u, Wt, k1, beta, delta as in class710_operator_derive.py's header, with delta = 1/2 k1 beta.
*   F1.1  Second-derivative terms: P Box h and K Box p, for every c.
*   F1.2  First-derivative terms: at c = 1 one only, K1 v.D p in the scalar equation.  For general c the rest is graviton-scalar
*         mixing proportional to (1 - c) K v.
*   F1.3  At c = 1, on three off-shell and three on-shell backgrounds with random K jets and with a = 1/12: E equals the
*         prediction R3 of class710_operator_derive.py (the same target the jet probes are compared with), 11 x 11; Tr Omega^2
*         equals step 3's.
*   F1.4  K2 appears in FORM's algebraic part and cancels in E (it enters C and D.omega with opposite signs).  No K3 anywhere.
*   F1.5  The ghost: Box eta + (Ric - c K v v) eta.
*   F1.6  At c = 1 the total by (3.40) obeys R6 pointwise and, on shell, R7: (203/320) K^2 (u^2)^2 after the Gauss-Bonnet and
*         Ricci-identity multiples.
*   F1.7  [R8 BY FORM -- can fail]  For c in {-1, 0, 1/2, 2, 3} as for c = 1, on six on-shell backgrounds each (a = 1/12):
*         alpha = 203/320, b = 0, zero residual on the two backgrounds beyond the four fitted.
*         alpha(c) is a polynomial of degree at most four in c (E is quadratic in c, the total quadratic in E; the ghost's part
*         quadratic).  Six equal values of a quartic make it constant.  So if F1.7 holds, alpha = 203/320 for EVERY c in this
*         family -- given FORM's operator and the degree count, and no longer a matter of sampling.
*         g, the Gauss-Bonnet coefficient: reported for each c; I expect 7/12 and do not know it.
*   F1.8  [THE TWO METHODS ON THE NEW CASE]  On the jet-probe runs' own backgrounds (regenerated from the same seed by the same
*         code), FORM's pointwise total equals the number in the author's logs: at c = 0 and c = 2 (five backgrounds each, the
*         R8 runs) and at c = 1 (the twelve of the main run).  The logs are read by the program; nothing is typed.
*         This is the comparison at c != 1 that F0 could not make, and the only independent test of E, omega, Omega there.
*   F1.9  [controls, must FAIL]  G = 1 in place of K at c = 1 breaks F1.2 and F1.3;  the commutator's sign reversed breaks
*         F1.3;  D_n B^m left out of the conversion breaks F1.3 whenever K1 != 0.
*   F1.10 [a question]  whether any of the four things F0 could not test (order of the two by-parts derivatives, the index
*         order in DT, the second-order term of GG, the unused Leibniz lines) becomes testable.  I expect not: K changes the
*         scalar sector only.  To be answered by mutation, not assumed.
*
* WHAT A PASS WOULD AND WOULD NOT SHOW.
*   Would: two methods of different kinds agreeing on the section 7 class -- exact numbers at a point from the Lagrangian by
*   jet probes, and a covariant symbolic derivation -- at c = 1 through a common target and at c = 0, 2 number against number;
*   and gauge-coefficient independence of the on-shell number for the whole c-family.
*   Would not: the literal gauge; anything shared downstream (hv340, the ghost's -2, the reductions, the on-shell substitution,
*   the R8 argument that alpha is the invariant); the local Jacobian; gauge families other than this one.
*   As the review of class710 said: on shell the K-theory's total at c equals the canonical theory's at c, so F1.7 is step 3's
*   gauge-coefficient independence, shown on the section 7 theory.
*
* COST.  FORM: seconds.  The reader: seconds to a minute (11 x 11 exact matrix algebra on about fifty backgrounds).
*
* WHAT HAPPENED (sandbox, FORM 4.3.1, 13 Nov; this block written 13 Nov 19:42 ET, clock read by the command that wrote it).
*   The program is gravscalar_operator.frm with K: the Lagrangian's K(Phi), the gauge factor GF, three Leibniz lines for
*   K, K1, K2, and a last module that writes D_jn B^jm.  0.00 s; 39 + 13 + 7 + 8 + 10 terms.
*   Every prediction F1.1-F1.9 held on the first run of the reader (form_class710_check.py, 44 s); see its header for the
*   numbers, for what F1.8 compared with the author's logs and checkpoints, and for the fresh review.
*   F1.10 answered by mutation: none of the four things F0 could not test becomes testable here (the reader's header).
*   The last module does not differentiate curvature: it relies on the first-derivative coefficients containing only v, K, K1
*   (true of the output; the reader refuses an output where it is not).
*   Options:  form -d TAG=literal -d LITERAL ...   the gauge condition's scalar term without K (reported; a control at c = 1);
*             form -d TAG=flip -d FLIP ...         the commutator's sign reversed (control).
*
* Run (WSL; from the repository's graviton_loop/form):
*   cd /mnt/c/GitHub/hdym-fakeon/graviton_loop/form
*   ~/form/sources/form class710_operator.frm
*   ~/form/sources/form -d TAG=literal -d LITERAL class710_operator.frm
*   ~/form/sources/form -d TAG=flip -d FLIP class710_operator.frm
* Each writes ../results/form_class710_<TAG>_{Eh,Ep,Gh,DBh,DBp}.out (overwritten; the evidence is the reader's log).
*
#-
#ifndef `TAG'
#define TAG "main"
#endif
Off statistics;
Dimension 4;
Symbols x, e, y, c, p, dp, pp, K, K1, K2, K3, KK, GF;
Indices a, b, m, n, r, s, i1, i2, i3, i4, i5, i6, i7, k1, k2, mu1, ia, ib, jm, jn;
Tensors h(symmetric), dh(symmetric), W(symmetric), Ric(symmetric), hh(symmetric);
Tensors H1, H2, H3, dH1, dH2, p1, p2, dp1, v, R, DW, DRic, S2, pS2, eta, E1, E2, SE2;
CFunctions GI, GG, DG, TT, DT, DD1, DD2, RF, CC, SG;

* ---------------------------------------------------------------------------------------------------------------------------
* the Lagrangian  sqrt(g)(R[g] - 1/2 K(Phi) g^{ab} d Phi d Phi) - 1/2 C.C,  g = g-bar + x h,  Phi = phi-bg + x p;  the x^2 part
* K(Phi) = K + x p K1 + 1/2 x^2 p^2 K2  (KK below);  the gauge condition's scalar term carries GF = K (image family) or 1 (literal)
* ---------------------------------------------------------------------------------------------------------------------------
Local L = SG*( GI(i1,i2)*RF(i1,i2) - 1/2*KK*GI(i3,i4)*(v(i3) + x*p1(i3))*(v(i4) + x*p1(i4)) )
        - 1/2*x^2*CC(i5)*CC(i5);
sum i1, i2, i3, i4, i5;
id KK = K + x*p*K1 + x^2*p^2*K2/2;

* (i) sqrt(g) and the inverse metric to second order
#ifdef `NOSQ'
id SG = 1 + x/2*h(i6,i6) + x^2*( -1/4*h(i6,i7)*h(i6,i7) );
#else
id SG = 1 + x/2*h(i6,i6) + x^2*( -1/4*h(i6,i7)*h(i6,i7) + 1/8*h(i6,i6)*h(i7,i7) );
#endif
sum i6, i7;
* (ii) the Ricci tensor of g-bar + x h, exactly, in terms of the difference of connections GG and its derivative DG(a,b,m,n) = D_n GG^a_bm
repeat;
  id once RF(b?,m?) = Ric(b,m) + DG(k1,b,m,k1) - DG(k1,b,k1,m) + GG(k1,k1,k2)*GG(k2,b,m) - GG(k1,m,k2)*GG(k2,b,k1);
  sum k1, k2;
endrepeat;
* the gauge condition
repeat;
  id once CC(m?) = H1(m,k1,k1) - 1/2*H1(k1,k1,m) - c*GF*p*v(m);
  sum k1;
endrepeat;
repeat;
  id once GI(a?,b?) = d_(a,b) - x*h(a,b) + x^2*h(a,k1)*h(k1,b);
  sum k1;
endrepeat;
* GG^a_bm = 1/2 g^{ad} x (D_b h_dm + D_m h_db - D_d h_bm) = 1/2 g^{ad} x TT(d,b,m), to second order
repeat;
  id once GG(a?,b?,m?) = x/2*TT(a,b,m) - x^2/2*h(a,k1)*TT(k1,b,m);
  sum k1;
endrepeat;
repeat;
  id once DG(a?,b?,m?,n?) = x/2*DT(a,b,m,n) - x^2/2*( H1(a,k1,n)*TT(k1,b,m) + h(a,k1)*DT(k1,b,m,n) );
  sum k1;
endrepeat;
#ifdef `LITERAL'
id GF = 1;
#else
id GF = K;
#endif
if (count(x,1) > 2) discard;
id TT(r?,b?,m?) = H1(r,m,b) + H1(r,b,m) - H1(b,m,r);
id DT(r?,b?,m?,n?) = H2(r,m,b,n) + H2(r,b,m,n) - H2(b,m,r,n);
if (count(x,1) != 2) discard;
id x = 1;
.sort

* ---------------------------------------------------------------------------------------------------------------------------
* (iv) the Euler-Lagrange expressions: first variation in the direction (dh, dp), derivatives moved off the variation
* ---------------------------------------------------------------------------------------------------------------------------
id h(a?,b?) = h(a,b) + y*dh(a,b);
id H1(a?,b?,m?) = H1(a,b,m) + y*dH1(a,b,m);
id H2(a?,b?,m?,n?) = H2(a,b,m,n) + y*dH2(a,b,m,n);
id p = p + y*dp;
id p1(a?) = p1(a) + y*dp1(a);
if (count(y,1) != 1) discard;
id y = 1;
.sort
* dH2(a,b,m,n) = D_n D_m dh_ab:  int X D_n D_m dh = int dh D_m D_n X  (D_n acts first)
id dH2(a?,b?,m?,n?) = dh(a,b)*DD1(n)*DD2(m);
id dH1(a?,b?,m?) = -dh(a,b)*DD1(m);
id dp1(m?) = -dp*DD1(m);
#do j = 1, 2
if (match(DD`j'(m?)));
  id h(a?,b?) = h(a,b) + e*H1(a,b,mu1);
  id H1(a?,b?,m?) = H1(a,b,m) + e*H2(a,b,m,mu1);
  id H2(a?,b?,m?,n?) = H2(a,b,m,n) + e*H3(a,b,m,n,mu1);
  id p = p + e*p1(mu1);
  id p1(a?) = p1(a) + e*p2(a,mu1);
  id v(a?) = v(a) + e*W(a,mu1);
  id W(a?,b?) = W(a,b) + e*DW(a,b,mu1);
  id Ric(a?,b?) = Ric(a,b) + e*DRic(a,b,mu1);
  id K = K + e*K1*v(mu1);
  id K1 = K1 + e*K2*v(mu1);
  id K2 = K2 + e*K3*v(mu1);
  if (count(e,1) != 1) discard;
  id e = 1;
  id DD`j'(m?) = replace_(mu1,m);
endif;
.sort
#enddo

* (iii) second derivatives: the symmetrised one plus half the commutator,  D_n D_m h_ab = S2(a,b,m,n) + 1/2 [D_n, D_m] h_ab
repeat;
#ifdef `FLIP'
  id once H2(a?,b?,m?,n?) = S2(a,b,m,n) - 1/2*( R(a,k1,n,m)*h(k1,b) + R(b,k1,n,m)*h(a,k1) );
#else
  id once H2(a?,b?,m?,n?) = S2(a,b,m,n) + 1/2*( R(a,k1,n,m)*h(k1,b) + R(b,k1,n,m)*h(a,k1) );
#endif
  sum k1;
endrepeat;
id p2(a?,b?) = pS2(a,b);
.sort
Local Eh = L;
Local Ep = L;
.sort
Skip L;
if (expression(Eh));
  id dp = 0;
  id dh(a?,b?) = d_(a,ia)*d_(b,ib);
endif;
if (expression(Ep));
  id dh(a?,b?) = 0;
  id dp = 1;
endif;
.sort
Skip L, Eh, Ep;

* ---------------------------------------------------------------------------------------------------------------------------
* the ghost: the gauge condition on a gauge transformation,  dh_ab = D_a eta_b + D_b eta_a,  dp = eta.v
* E1(a,m) = D_m eta_a,  E2(a,m,n) = D_n D_m eta_a
* ---------------------------------------------------------------------------------------------------------------------------
#ifdef `LITERAL'
Local Gh = H1(ia,k1,k1) - 1/2*H1(k1,k1,ia) - c*p*v(ia);
#else
Local Gh = H1(ia,k1,k1) - 1/2*H1(k1,k1,ia) - c*K*p*v(ia);
#endif
sum k1;
id H1(a?,b?,m?) = E2(a,b,m) + E2(b,a,m);
id p = eta(k2)*v(k2);
sum k2;
repeat;
#ifdef `FLIP'
  id once E2(a?,m?,n?) = SE2(a,m,n) - 1/2*R(a,k1,n,m)*eta(k1);
#else
  id once E2(a?,m?,n?) = SE2(a,m,n) + 1/2*R(a,k1,n,m)*eta(k1);
#endif
  sum k1;
endrepeat;
.sort
Skip L, Eh, Ep, Gh;

* ---------------------------------------------------------------------------------------------------------------------------
* D_jn B^jm: the covariant derivative of the coefficient of the first-derivative terms (D_jm h, D_jm p), the fields' slots held
* fixed (hh, pp are not differentiated).  The reader needs it for D.omega and Omega.
* ---------------------------------------------------------------------------------------------------------------------------
Local DBh = Eh;
Local DBp = Ep;
.sort
Skip L, Eh, Ep, Gh;
if ((count(H1,1) == 0) && (count(p1,1) == 0)) discard;
id H1(a?,b?,m?) = hh(a,b)*d_(m,jm);
id p1(m?) = pp*d_(m,jm);
id v(a?) = v(a) + e*W(a,jn);
id W(a?,b?) = W(a,b) + e*DW(a,b,jn);
id Ric(a?,b?) = Ric(a,b) + e*DRic(a,b,jn);
id K = K + e*K1*v(jn);
id K1 = K1 + e*K2*v(jn);
id K2 = K2 + e*K3*v(jn);
if (count(e,1) != 1) discard;
id e = 1;
.sort
#write <../results/form_class710_`TAG'_Eh.out> "%E", Eh
#write <../results/form_class710_`TAG'_Ep.out> "%E", Ep
#write <../results/form_class710_`TAG'_Gh.out> "%E", Gh
#write <../results/form_class710_`TAG'_DBh.out> "%E", DBh
#write <../results/form_class710_`TAG'_DBp.out> "%E", DBp
#message TAG `TAG': five files written to ../results
.end
