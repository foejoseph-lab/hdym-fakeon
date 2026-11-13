* gravscalar_operator.frm
* Rung 3, step 3b, second method for route 2: FORM.  STAGE F0 -- THE KNOWN CASE (step 3, K = 1) FIRST.
* PREDICTIONS WRITTEN 13 Nov 2026 17:11 ET (clock read by the command that wrote them), before any FORM code for the physics
* existed; the author's yes 17:25 ET.  They stand as first written (rule 11); what happened follows them.
* Working names: graviton_loop/form/gravscalar_operator.frm (FORM) and graviton_loop/form_operator_check.py (reads FORM's output).
*
* WHY F0.  The author's rule of 10 Nov: a new method stands beside a proven one on a case both can do before it is trusted on a
* new case.  Step 3's gauge-fixed operator at K = 1 is known three ways already: the print (4.31)-(4.33), the jet-probe derivation
* (gravscalar_operator_derive.py, the author's PASS run), and hv1974_gravscalar_M.build_E.  FORM must reproduce it before it
* touches the class (7.10).
*
* WHAT FORM DOES, AND FROM WHAT.  It is a different kind of method from the Python one: symbolic, covariant, valid on any
* background, with c kept as a SYMBOL; the Python method is exact numbers at one point on random backgrounds.
*   Fed to FORM (each an identity of tensor calculus, not a result of the paper):
*     (i)   g^{ab} = g-bar^{ab} - h^{ab} + h^{ac} h_c^b  and  sqrt(g) = sqrt(g-bar)(1 + 1/2 h - 1/4 h_ab h^ab + 1/8 h^2), to second order;
*     (ii)  the exact difference of connections, delta-Gamma^a_bc = 1/2 g^{ad}(D_b h_dc + D_c h_db - D_d h_bc), and the exact
*           Ricci tensor of g-bar + h in terms of it (D delta-Gamma - D delta-Gamma + delta-Gamma delta-Gamma - ...);
*     (iii) the commutator of two covariant derivatives on a tensor (the Ricci identity), in the curvature convention of
*           gravscalar_operator_derive.py -- the one place a sign convention enters;
*     (iv)  the covariant Euler-Lagrange formula  dL/dF - D_c dL/d(D_c F) + D_c D_d dL/d(D_c D_d F).
*   Derived by FORM: the second-order Lagrangian of sqrt(g)(R - 1/2 (d phi)^2) - 1/2 sqrt(g-bar) C^2 with
*     C_m = D^n h_mn - 1/2 D_m h - c phi d_m phi-bg; its Euler-Lagrange expressions for h_ab and phi; derivatives brought to a
*     fixed order with (iii); the result written to a file as a sum of tensor monomials in Riem, v = d phi-bg, W = DD phi-bg, delta.
*   Left to Python: evaluating FORM's written tensor expression on random backgrounds (plain index sums in exact rationals) and
*     comparing.  FORM is not asked to canonicalise Riemann monomials -- its weak point -- and nothing is typed from one log to another.
*   NOT SETTLED until a sandbox trial: whether (iv) with the reordering (iii) is practical in FORM as described.  If it is not, this
*     header is amended and shown again before the plan changes; the fallback is FORM for (i)-(ii) and the bilinear form only.
*
* PREDICTIONS (rule 25).
*   F0.1  Second-derivative terms.  For every c the h-equation's are P^{ab,cd} Box h_cd (P = 1/2 (gg + gg) - 1/2 g g in the
*         paper's weights) and the scalar's Box phi: the terms D_a D^c h_cb, D_a D_b h and g_ab D.D.h cancel identically
*         between the action and the gauge term.  (That is what the de Donder part of the gauge is for; c does not enter.)
*   F0.2  First-derivative terms.  None at c = 1.  For general c they are the graviton-scalar mixing only and are proportional
*         to (1 - c) v.  [The Python c = 0 smoke of 13 Nov found 60 nonzero first-derivative entries at a = 1/12, 4 at c = 1.]
*   F0.3  The algebraic part at c = 1, evaluated on three random backgrounds (off shell), equals hv1974_gravscalar_M.build_E as
*         an 11 x 11 matrix, exactly.  This is the content of F0.
*   F0.4  The ghost, from C under delta h = D eta + D eta, delta phi = eta.d phi-bg:  Box eta_m + (Ric - c v v)_m^n eta_n, no
*         first-derivative terms, for symbolic c.
*   F0.5  [control, must FAIL F0.3]  the sign of the Riemann term in (iii) flipped;  [control]  the 1/8 h^2 of sqrt(g) dropped.
*   F0.6  [a question, not a prediction]  how many terms and how many seconds.  I expect seconds.
*
* WHAT A PASS WOULD AND WOULD NOT SHOW.
*   Would: that FORM, from the action and four identities, gives the operator the print and the jet probes give -- a symbolic
*   derivation beside a numerical one; and a formula in c, which the Python method can only sample.
*   Would not: anything new about the physics (every number here is known); the traces, (3.40), the Faddeev-Popov -2, the
*   reductions, the on-shell substitution (all still Python's, shared); FORM on the class (7.10) -- that is F1, with its own header.
*   A convention error common to (iii) and to build_E's target would be invisible; the control F0.5 shows the sign is at least seen.
*
* AFTER F0 (not part of this yes): F1, the class (7.10) with K(phi) and symbolic c -- the targets are R1, R3, R5 of
* class710_operator_derive.py and, at general c, a covariant formula for the operator the R8 runs sample at c = 0 and 2.
*
* WHAT HAPPENED (sandbox, FORM 4.3.1, 13 Nov; this block written 13 Nov 17:30 ET, clock read by the command that wrote it).
*   The plan worked as described; nothing in the header had to change.  Run time 0.00 s; 35 + 8 terms for the two equations.
*   F0.1, F0.2, F0.3 held on three random backgrounds on the first run (form_operator_check.py).  See that script's header for
*   the controls and the ghost.
*   HOW IT IS DONE.  h, D h, D D h are separate tensors h, H1(a,b,m) = D_m h_ab, H2(a,b,m,n) = D_n D_m h_ab (derivative
*   indices in the order applied); likewise p, p1, p2 for the scalar and v = D phi-bg, W = D D phi-bg.  Derivatives of a
*   product (Leibniz) by the substitution f -> f + e D f kept to first order in e.  The variation by f -> f + y df kept to first
*   order in y; then the derivatives are moved off df by parts.  Second derivatives are written as the symmetrised one, S2, plus
*   half the commutator; S2, H1, h are then independent data at a point, so the reader can test each order by itself.
*   CONVENTION (iii):  [D_m, D_n] V_a = R(a,e,m,n) V_e,  Ric(b,n) = R(a,b,a,n)  -- gravscalar_operator_derive.py's, whose
*   metric is delta - 1/3 R(mu,a,nu,b) y^a y^b.
*   Options (the reader's controls):  form -d TAG=flip -d FLIP ...   the commutator's sign reversed;
*                                     form -d TAG=nosq -d NOSQ ...   the 1/8 h^2 of sqrt(g) dropped.
*
*   AFTER THE FRESH REVIEW (see form_operator_check.py's header for its table): typed correctly by reading; 17 of 21 mutations
*   caught.  F0 cannot test, because here they change nothing: the order of the two by-parts derivatives ("D_n acts first"
*   below), the derivative-index order in DT's first term, the x^2 term of GG, and the Leibniz lines for H3, DW, DRic.
*   F0.1's "cancel identically" is true once the symmetries of S2 are imposed, which the reader does and FORM does not.
*
* Run (WSL; from the repository's graviton_loop/form):
*   cd /mnt/c/GitHub/hdym-fakeon/graviton_loop/form
*   ~/form/sources/form gravscalar_operator.frm
*   ~/form/sources/form -d TAG=flip -d FLIP gravscalar_operator.frm
*   ~/form/sources/form -d TAG=nosq -d NOSQ gravscalar_operator.frm
* Each writes ../results/form_gravscalar_<TAG>_{Eh,Ep,Gh}.out (overwritten: they are FORM's output, regenerated in a blink;
* the evidence is form_operator_check.py's log, which is refused if it exists).
*
#-
#ifndef `TAG'
#define TAG "main"
#endif
Off statistics;
Dimension 4;
Symbols x, e, y, c, p, dp;
Indices a, b, m, n, r, s, i1, i2, i3, i4, i5, i6, i7, k1, k2, mu1, ia, ib;
Tensors h(symmetric), dh(symmetric), W(symmetric), Ric(symmetric);
Tensors H1, H2, H3, dH1, dH2, p1, p2, dp1, v, R, DW, DRic, S2, pS2, eta, E1, E2, SE2;
CFunctions GI, GG, DG, TT, DT, DD1, DD2, RF, CC, SG;

* ---------------------------------------------------------------------------------------------------------------------------
* the Lagrangian  sqrt(g)(R[g] - 1/2 g^{ab} d Phi d Phi) - 1/2 C.C,  g = g-bar + x h,  Phi = phi-bg + x p;  the x^2 part
* ---------------------------------------------------------------------------------------------------------------------------
Local L = SG*( GI(i1,i2)*RF(i1,i2) - 1/2*GI(i3,i4)*(v(i3) + x*p1(i3))*(v(i4) + x*p1(i4)) )
        - 1/2*x^2*CC(i5)*CC(i5);
sum i1, i2, i3, i4, i5;

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
  id once CC(m?) = H1(m,k1,k1) - 1/2*H1(k1,k1,m) - c*p*v(m);
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
Local Gh = H1(ia,k1,k1) - 1/2*H1(k1,k1,ia) - c*p*v(ia);
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
#write <../results/form_gravscalar_`TAG'_Eh.out> "%E", Eh
#write <../results/form_gravscalar_`TAG'_Ep.out> "%E", Ep
#write <../results/form_gravscalar_`TAG'_Gh.out> "%E", Gh
#message TAG `TAG': wrote ../results/form_gravscalar_`TAG'_Eh.out, _Ep.out, _Gh.out
.end
