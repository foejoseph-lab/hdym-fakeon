#!/usr/bin/env python3
"""
Rung 3, step 3c: Deser & van Nieuwenhuizen 1974, Einstein-Maxwell -- the operator from the action, (3.40), the paper's (4.16a), 137/60.
HEADER ONLY.  NO CODE EXISTS YET.  PREDICTIONS WRITTEN 13 Nov 2026 22:37 ET (clock read by the command that wrote this file),
for the author's yes before any code.  They stand as first written (rule 11); corrections will follow them, dated.

WHAT THE PRINT SAYS (Phys. Rev. Lett. 32, 245 and Phys. Rev. D 10, 401, read in full as images 13 Nov; rung 3 note 4n).
  (2.1)  L = -sqrt(-g)[kappa^-2 R + 1/4 F^2], signature +2, their R minus the standard one (my inference from three places).
  (2.2a) g-bar = g + kappa h,  A-bar = kappa^-1 A + a: kappa drops out of the quadratic part.
  (2.4)  background: G = -1/2 T, D.F = 0, T_mn = F_ma F_n^a - 1/4 g_mn F^2;  (2.5) R = 0, 4 Ric^2 = T^2, Ric.T = -2 Ric^2.
  (3.1)  THE GAUGE IS STATED: -1/2 (D^n h_mn - 1/2 D_m h)^2 - 1/2 (D.a)^2.  Plain de Donder and Lorentz.  No matter term in the
         gravitational condition (where 't Hooft-Veltman's (4.27) has -phi d phi-bg).
  (2.6)-(2.7d), (3.10): the graviton-photon mixing h_ab Q^{ab rs} D_r a_s stays, a first-derivative term.
  (3.12)-(3.13): a vector ghost and a scalar ghost, mixed in one corner only.
  (4.11)-(4.13) the traces and the non-ghost part; (4.15) the ghost; (4.16a) the total off shell; (4.16b) on shell 137/60 Ric^2.
  Printed slips (note 4n): the letter's c2 = -29/10 (the paper: -1/40); (A6)'s 5/14 (from 't Hooft-Veltman's (5.24): 1/8).

SET-UP (this repository's conventions: Euclidean, standard curvature sign, as in steps 2 and 3).
  L = sqrt(g)( R - 1/4 F_mn F^mn ),  metric g + h,  potential A + a.  Background equations: Ric = 1/2 T, R = 0, D.F = 0.
  Terms linear in the curvature change sign against the print (Ric.T, R F^2, Riem.F.F); squares do not.
  Fields: 14 components, h (10, the basis of steps 2-3) and a (4).  Ghosts: 5, eta_m and eta_5.
  Background at the origin: a random algebraic curvature tensor (normal coordinates), F_mn (6 rationals), D_a F_bc (20 rationals,
  the cyclic identity imposed), from a polynomial potential A(x) to third order.  Second derivatives of F cannot reach the origin.
  Operator form -(S D^2 + A^m D_m + B) from the Lagrangian alone by jet probes with first-order data (class710_operator_derive's
  first_order(), operator_data(), checkpoint and 'expand' variant with its cross-check against 'translate'), then S^-1, omega, E,
  Omega, (3.40) with its first-derivative term.  14 x 15 = 210 base fields of degree <= 2.

STAGE A: THE PAPER'S GAUGE (3.1).   "Steps 1-3" are this repository's closed results.
  A0 [validation]  F = 0, DF = 0: the operator is block diagonal -- step 2's h operator and the Feynman-gauge vector of
        maxwell_operator_check (E = -Ric there); no first-derivative terms.
  A1  S = diag(P, delta).  First-derivative terms in the mixed h-a entries ONLY, linear in F, no curvature:
        the cross term of L_2 is, up to one overall sign I have not fixed,
            h^{am} F_m^b (D_a a_b - D_b a_a) - 1/2 h F^{rs} D_r a_s          [hand expansion of -1/4 sqrt(g) F-bar^2; matches (2.7d)]
        so omega is off-diagonal, antisymmetric with respect to S, and NOT a gradient: Tr Omega^2 depends on F and DF.
  A2  E: the h-h block is step 2's E plus terms quadratic in F ((2.7c) and omega^2); the a-a block is -Ric plus terms quadratic
        in F; the mixed blocks are linear in DF.  Block by block I predict only this structure; the content is A3.
  A3 [the traces, as printed, their sign of R]  with F2_mn = F_ma F_n^a, F^4 = tr F2^2:
        (4.11a)+(4.11b)+(4.11c):  Tr Y^2 = -28 Ric^2 + 7 R^2 - R F^2 - 2 Ric.F2 + 17/2 F^4 - 11/2 (F^2)^2 - 9/2 (D_a F_bc)^2 + (D.F)^2
        (4.12a)+(4.12b)+(4.12c):  Tr X^2 = 7 Ric^2 + 1/18 R^2 + Ric.F2 + 3/4 F^4 + 3/8 (F^2)^2 + 3/4 (D_a F_bc)^2 + 1/2 (D.F)^2
        [(4.12c) once: the two off-diagonal blocks of X are -dQ and -1/2 dQ, and their two products add to (dQ)^2]
        [the sum of the printed lines is mine; whether the run's Tr Omega^2 and Tr (E - R/6)^2 match them term by term depends on
         total derivatives the paper may have dropped inside (4.11)-(4.12).  If they do not match but A4 does, that is the answer.]
  A4 [non-ghost, (4.13)]  1/2 [ 1/12 Tr Y^2 + 1/2 Tr X^2 + 14/60 (Ric^2 - R^2/3) ], pointwise, fitted on ten structures
        (Ric^2, R^2, Riem^2; F^4, (F^2)^2; R F^2, Ric.F2, Riem.F.F; (D_a F_bc)^2, (D.F)^2):
            7/10 Ric^2 + 4/15 R^2 + 1/6 Ric.T + 13/24 T^2 + 1/6 (D.F)^2  + g_A Gauss-Bonnet,
        with ZERO coefficient on (D_a F_bc)^2 and on Riem.F.F before any reduction (the paper's -3/8 + 3/8), and zero on
        (F^2)^2 and R F^2 beyond what T contains.  [can fail: the paper may have used (4.3)-(4.4) where I assume it did not.
        Then the pointwise total differs from this line by a multiple of the identity (4.4), and I will say so, not refit quietly.]
  A5 [ghost, derived]  from C_m and C_5 under delta h = D eta + D eta, delta a_m = eta^a F_am + D_m eta_5 ((2.10a), (2.10b')):
        vector E_gh = Ric (standard sign; step 3's with the scalar removed), scalar E = 0, and ONE mixed corner:
        (eta_5 row, eta_m column) = F_am D^m + D^m F_am.  Nothing in the other corner.
        The gauge symmetry of the ungauged operator under these variations, on shell, is checked first (as on 10 Nov).
  A6 [ghost contribution, (4.15)]  -( 1/4 Ric^2 + 7/24 R^2 ) after Gauss-Bonnet, INDEPENDENT of F and DF.
        This one is a proof, not a prediction: a block-triangular operator has the divergence of its diagonal blocks.  The run
        checks the code against it.  The numbers are (5.21) at zero scalar plus one complex minimal scalar.
  A7 [total off shell, (4.16a)]  9/20 Ric^2 - 1/40 R^2 + 1/6 Ric.T + 13/24 T^2 + 1/6 (D.F)^2      [their sign of R; ours: -1/6 Ric.T]
        i.e. c1 = 9/20, c2 = -1/40 (NOT the letter's -29/10), c4 = 13/24, c7 = 1/6, c9 = 1/6, c5 = c6 = c8 = 0.
  A8 [THE NUMBER]  on shell (Ric = 1/2 T, R = 0, D.F = 0), pointwise after Gauss-Bonnet:  (137/60) Ric^2.
  A9 [a question]  the Gauss-Bonnet coefficient g: it is the coefficient of Riem^2, so it cannot depend on F; I expect the F = 0
        value, which the script computes from steps 1-2's own functions.  Reported.
  Controls (must FAIL): omega halved; the vector ghost's Ric with the other sign (the paper remarks on this sign); the photon's
        scalar ghost left out (the paper remarks that it contributes: the total must move by one complex scalar loop).

THE PROTOCOL WHEN THE RUN AND THE PRINT DIFFER (author, 13 Nov 22:35 ET).
  The paper wins: the code or the reading is wrong; full stop until why and by how much is known.
  We win WITH the cause shown on the paper's side (from its own pieces): the printed value goes to the ledger verdict with the
  reason; the supported value is the target afterwards.  We differ and cannot say why: unresolved, stays O.  Print is never deleted.
  Already decided on paper: c2 (the paper's own pieces give -1/40), (A6) (no 14 can arise), the power of lambda in (5.4).

STAGE B: A GAUGE FAMILY THAT IS NOT IN THE PAPER (ours; the author's yes 13 Nov 22:35 ET; outside the protocol above).
  C_m = D^n h_mn - 1/2 D_m h - c F_mn a^n: 't Hooft-Veltman's kind of matter term, carried to the Maxwell case.  c = 0 is the paper.
  B1 [hand algebra, signs unchecked]  at c = 1 the gauge term cancels the two pieces of the mixing that 't Hooft-Veltman's term
        cancels for the scalar,  h^{am} F_m^b D_a a_b  and  -1/2 h F^{rs} D_r a_s,  and leaves  -h^{am} F_m^b D_b a_a.
        So NO value of c removes the first-derivative mixing: unlike the scalar case, nothing clicks into place.
  B2  for c != 0 the ghost is no longer triangular (the vector row gains -c F_mn D^n eta_5 and -c F_mn F^{an} eta_a): it sees F.
  B3 [can fail]  on shell the divergence is the paper's along c, as an integral: (137/60) Ric^2 up to total derivatives.
  B is NOT coded with stage A.  It needs first the list of on-shell invariants and total derivatives for (Weyl, F, DF), as R8
  needed on 13 Nov, and gets its own header block and the author's yes then.  Any coefficient it produces is ours, labelled so,
  and is never set beside the paper's off-shell numbers.

WHAT A PASS OF STAGE A WOULD AND WOULD NOT SHOW.
  Would: 137/60 and the five off-shell coefficients from the Lagrangian alone, in the paper's stated gauge, with the three zeros;
  the first-order machinery on an operator whose connection is not a gradient (its first real test: on 13 Nov omega was a gradient
  in the image gauge); the letter's c2 against the paper's.
  Would not: anything independent of steps 1-3 in the Faddeev-Popov sign, the Gauss-Bonnet reduction, the on-shell substitution;
  the overall sign of the mixing (A -> -A is a symmetry: every trace is even in F); gauge independence (stage B); the Jacobian.
  One method only until FORM has done the same (rule 12): FORM's header is written before its code -- a vector field with its
  own gauge condition and ghost, the Leibniz and commutator pattern of 13 Nov, the reader comparing entry by entry.
  DeWitt's recursion (ledger DN-2) needs the bundle connection to carry omega; not attempted here; it may stumble there.

COST.  No forecast.  210 base fields against route 2's 176, a heavier Lagrangian (two inverse metrics on F-bar^2).  A sandbox
smoke on a few columns first, then a five-minute smoke on the author's machine, and the forecast from that.  Resumable
(<log>.ckpt, keys carrying a hash of the code); logs refused, not overwritten; END line.
"""
raise SystemExit("einstein_maxwell_operator_derive.py: header only (13 Nov 2026). No code yet; see the docstring.")
