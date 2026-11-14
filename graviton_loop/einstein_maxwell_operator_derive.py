#!/usr/bin/env python3
"""
Rung 3, step 3c: Deser & van Nieuwenhuizen 1974, Einstein-Maxwell -- the operator from the action, (3.40), the paper's (4.16a), 137/60.
PREDICTIONS WRITTEN 13 Nov 2026 22:37 ET (clock read by the command that wrote this file),
before any code existed, for the author's yes.  They stand as first written (rule 11); corrections follow them, dated.

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

AS CODED, AND THE RECORD (this block written 13 Nov 2026 23:31 ET, clock read by the command that inserted it; the author's yes to the
header 13 Nov 22:42 ET).  HEADER ONLY no longer: stage A is coded below.  Stage B is not.
  A3 was corrected before the header was shown: (4.12c) counts once, not twice (the line above carries the reason).
  The background potential: A_n = 1/2 F_mn y^m + 1/2 q_n,mr y^m y^r + cubic, so D_r F_mn(0) = q_n,mr - q_m,nr; on shell the four
  conditions D.F = 0 are imposed on q and the Ricci part of the curvature is set to T/2 (with_ricci).  The cubic terms are drawn
  nonzero on every background with F != 0 and are in no prediction.
  A4-A8 are tested POINTWISE against the predicted formula on each background, not by a ten-structure fit: a fit needs twelve
  backgrounds and is only wanted if a prediction fails.  The F = 0 part of every prediction is computed in the run from
  hv1974_gravscalar_M's own functions (nothing typed); the F part is the paper's three numbers.  The ten structures, Tr Omega^2,
  Tr (E + R/6)^2 and Tr E are printed for every background (the raw data are kept).
  A5a, the gauge symmetry of the ungauged operator: on pure-gauge fields with eta linear and quadratic only (eta(0) = 0), because
  the backgrounds satisfy the field equations AT the origin and not their first derivatives (D Ric = 0 there, D T is not).
  The vector ghost is ghost_operator_derive.ghost_operator on the same metric (the author's PASS run of 10 Nov); the scalar ghost's
  row is new code (photon_ghost_row).  A6 is the proof stated above; the code adds the two diagonal blocks.
  SMOKE RECORD (rule 11; sandbox, 2 cores, python 3.13.18, sympy 1.13.1, python-flint 0.9.0).
  0 (four single base fields, by hand): the measured first-derivative columns equal the hand cross term of A1, overall sign +.
  1 (one off-shell background, --mode min --crosscheck 0 --symmetry 0 --jobs 2, 1990 s; started 22:50 ET):
     A4 holds: (3.40) on the 14 fields = steps 1-2 at F = 0 - 1/6 Ric.T + 13/24 T^2 + 1/6 (D.F)^2, pointwise, with NO multiple of
     (D_a F_bc)^2 or Riem.F.F and nothing on (F^2)^2 or R F^2 beyond T.  A5 holds.  A7 holds: the total is (4.16a) + (13/24)
     Gauss-Bonnet -- c2 = -1/40, the paper's, not the letter's.  The three controls break it.
     A1 printed FAILS: its last clause, "S omega antisymmetric", was mis-coded by me -- each off-diagonal h component is one
     variable standing for two entries, so the bilinear form is Wd S with Wd = 2 on those; with Wd the measured term is
     antisymmetric (checked on the predicted term, which the measured one equals).  Corrected in the code after this smoke.
  2 (one on-shell background, with the symmetry test and the cross-check; forecast 50 minutes): STOPPED at the author's word
     (23:31 ET: a test known from the start to take 50 minutes is not a sandbox smoke; planned tests here are 20 minutes at most).
     Two base fields done; nothing learned.  So A0, A5a, A8 (the 137/60 itself), the cross-check of 'expand' against 'translate'
     and --mode all have NOT been run anywhere.
  Per base field on one sandbox core: 'plain' 14 s alone, 31 s with both cores busy; 'expand' 29 s alone, 50 s with both busy.
  NOT REVIEWED (rule 9).  A fresh-agent review was started at 23:28 ET and stopped at the author's word at 23:31 before it
  reported anything.  This file goes out unreviewed and says so; nothing here is to be leaned on until it is reviewed.

Run (from graviton_loop; the author's machine).  First the smoke, about five minutes by my arithmetic from the sandbox
(84 base fields; forecasts for this machine have been wrong in both directions -- it prints a timing line):
      python einstein_maxwell_operator_derive.py --jobs 16 --onshell 1 --offshell 0 --f0 0 --mode min --crosscheck 0 --symmetry 0 --log results/einstein_maxwell_smoke_onshell.log
  It ends status=INCOMPLETE by design if every check holds; the line that matters is A8.
  The full run (2 on-shell + 2 off-shell + 1 F = 0 backgrounds, every base field, the symmetry test, the cross-check) is not to be
  started before the smoke's timing line is read and the review is done.  Resumable: the same command with --resume added.
"""
import sys, os, time, argparse, random, itertools, platform, json, hashlib
import sympy as sp
import flint

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import class710_operator_derive as CL          # the ring, the truncations, first_order(), the checkpoint (the author's PASS runs of 13 Nov)
import gravscalar_operator_derive as D
import ghost_operator_derive as GH
import hv1974_gravscalar_M as GM

idx = range(4)
PAIRS = D.PAIRS
COMPS = PAIRS + [('a', m) for m in idx]          # 10 + 4 = 14
NC = len(COMPS)
CTX, Y, TT, EPS, DL, ZERO, ONE, HALF = CL.CTX, CL.Y, CL.TT, CL.EPS, CL.DL, CL.ZERO, CL.ONE, CL.HALF
C, rat, psum, d, lin_expand, mono = CL.C, CL.rat, CL.psum, CL.d, CL.lin_expand, CL.mono
out = None
SELF = 'einstein_maxwell_operator_derive.py'


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def _code_version():
    h = hashlib.sha1()
    here = os.path.dirname(os.path.abspath(__file__))
    for fn in (SELF, 'class710_operator_derive.py', 'gravscalar_operator_derive.py', 'ghost_operator_derive.py', 'hv1974_gravscalar_M.py'):
        with open(os.path.join(here, fn), 'rb') as fh:
            h.update(fh.read().replace(b'\r\n', b'\n'))
    return h.hexdigest()[:12]


CODE_VERSION = _code_version()


def ypow(m):
    return Y[0]**m[0] * Y[1]**m[1] * Y[2]**m[2] * Y[3]**m[3]


class Background(CL.Background):
    """the metric of class710_operator_derive.Background (no scalar) and a background potential
        A_n(y) = 1/2 F_mn y^m + 1/2 q_n,mr y^m y^r + cubic,    so that  F_mn(0) = F_mn,  d_r F_mn(0) = q_n,mr - q_m,nr  (= D_r F_mn(0))."""
    def __init__(self, R, F, Q, cubA, tdeg, o=0):
        super().__init__(R, sp.zeros(4, 1), sp.zeros(4, 4), {}, [1, 0, 0], tdeg, o)
        A = []
        for n in idx:
            p = psum(C(F[m, n] / 2) * Y[m] for m in idx) + psum(C(Q[n][(min(m, r), max(m, r))] / 2) * Y[m] * Y[r] for m in idx for r in idx)
            p = p + psum(C(c) * ypow(mm) for mm, c in cubA[n].items())
            if tdeg:
                p = p + psum(TT[m] * d(p, m) for m in idx)
            A.append(p)
        self.Fbg = [[d(A[n], m) - d(A[m], n) for n in idx] for m in idx]

    def L2(self, h, a, gauge=True):
        """the eps^2 coefficient of  sqrt(G)(R[G] - 1/4 G^{ma} G^{nb} Fb_mn Fb_ab)  [G = g + eps h, Fb = F-bg + eps (d a - d a)]
           - 1/2 sqrt(g) g^{mn} C_m C_n - 1/2 sqrt(g) (g^{mn} D_m a_n)^2,   C_m = g^{nr} D_r h_mn - 1/2 D_m h    [(3.1); gauge=False: left out]"""
        T = self.T
        G = [[self.gb[i][j] + EPS * h[i][j] for j in idx] for i in idx]
        Ginv = self.series_inverse(G, 3)
        A = [[G[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]; A2 = self.mat_mul(A, A, 2); A3 = self.mat_mul(A2, A, 2)
        tr = lambda M: psum(M[i][i] for i in idx)
        logdet = T(tr(A) - tr(A2) * C(HALF) + tr(A3) * C(sp.Rational(1, 3)), 2)
        l2 = T(logdet * logdet, 2)
        sqrtg = T(ONE + logdet * C(HALF) + l2 * C(sp.Rational(1, 8)) + T(l2 * logdet, 2) * C(sp.Rational(1, 48)), 2)
        dG = [[[T(d(G[i][j], c), 3) for c in idx] for j in idx] for i in idx]
        Gam = [[[T(psum(Ginv[m][l] * (dG[l][j][i] + dG[l][i][j] - dG[i][j][l]) for l in idx) * C(HALF), 3)
                 for j in idx] for i in idx] for m in idx]
        Ric = [[None] * 4 for _ in idx]
        for b in idx:
            for c in idx:
                if c < b:
                    Ric[b][c] = Ric[c][b]; continue
                Ric[b][c] = T(psum(d(Gam[i][b][c], i) - d(Gam[i][b][i], c) for i in idx)
                              + psum(Gam[i][i][e] * Gam[e][b][c] - Gam[i][c][e] * Gam[e][b][i] for i in idx for e in idx), 2)
        Rsc = T(psum(Ginv[b][c] * Ric[b][c] for b in idx for c in idx), 2)
        Fb = [[T(self.Fbg[m][n] + EPS * (d(a[n], m) - d(a[m], n)), 2) for n in idx] for m in idx]
        F1 = [[T(psum(Ginv[m][i] * Fb[m][n] for m in idx), 2) for n in idx] for i in idx]            # F^i_n
        F2 = [[T(psum(F1[i][n] * Ginv[n][j] for n in idx), 2) for j in idx] for i in idx]            # F^ij
        FF = T(psum(F2[i][j] * Fb[i][j] for i in idx for j in idx), 2)
        L = T(sqrtg * (Rsc - FF * C(sp.Rational(1, 4))), 2)
        if gauge:
            Gb, gbinv = self.Gb, self.gbinv
            Dh = [[[T(d(h[i][j], c) - psum(Gb[l][c][i] * h[l][j] + Gb[l][c][j] * h[i][l] for l in idx), 3)
                    for c in idx] for j in idx] for i in idx]
            Dtrh = [T(psum(gbinv[i][j] * Dh[i][j][c] for i in idx for j in idx), 3) for c in idx]
            Cl = [T(psum(gbinv[n][r] * Dh[m][n][r] for n in idx for r in idx) - Dtrh[m] * C(HALF), 3) for m in idx]
            C2 = T(self.sqrtgb * psum(gbinv[m][n] * Cl[m] * Cl[n] for m in idx for n in idx), 2)
            Da = T(psum(gbinv[m][n] * (d(a[n], m) - psum(Gb[l][m][n] * a[l] for l in idx)) for m in idx for n in idx), 3)
            L5 = T(self.sqrtgb * T(Da * Da, 2), 2)
            L = T(L - C(HALF) * EPS * EPS * (C2 + L5), 2)
        return CTX.from_dict({m[:8] + (0,) + m[9:]: c for m, c in L.to_dict().items() if m[8] == 2})


def comp_field(comp, mon):
    h = [[ZERO for _ in idx] for _ in idx]; a = [ZERO for _ in idx]
    if comp[0] == 'a':
        a[comp[1]] = mon
    else:
        i, j = comp
        h[i][j] = mon
        if i != j:
            h[j][i] = mon
    return h, a


def probe_fields():
    P = []
    for comp in COMPS:
        P.append(('0', comp, None, comp_field(comp, ONE)))
        for m in idx:
            P.append(('1', comp, m, comp_field(comp, Y[m])))
        for m in idx:
            for n in idx:
                if m <= n:
                    P.append(('2', comp, (m, n), comp_field(comp, Y[m] * Y[n] * (C(HALF) if m == n else ONE))))
    return P


def el(bg, base, probes, gauge=True):
    """the 14 Euler-Lagrange components on the base field, each as (value at the origin, four first derivatives);
    class710_operator_derive.el with 14 components"""
    bh, ba = base
    T = bg.T
    G = {}
    for kind, comp, key, (ph, pa) in probes:
        h = [[bh[i][j] + DL * ph[i][j] for j in idx] for i in idx]
        a = [ba[i] + DL * pa[i] for i in idx]
        G[(kind, comp, key)] = CL.coeff_dl(bg.L2(h, a, gauge))
    half = C(HALF)
    outv = []
    for comp in COMPS:
        G0 = G[('0', comp, None)]
        G1 = [G[('1', comp, m)] for m in idx]
        expr = G0
        for m in idx:
            expr = expr - d(G1[m] - Y[m] * G0, m)
        for m in idx:
            for n in idx:
                mm, nn = min(m, n), max(m, n)
                w = ONE if m == n else half
                F2 = Y[m] * Y[n] * G0 * half - (Y[n] * G1[m] + Y[m] * G1[n]) * half + w * G[('2', comp, (mm, nn))]
                expr = expr + d(d(F2, m), n)
        v0, g1 = lin_expand(T(expr, 0), 4 if bg.tdeg else 0)
        if comp[0] != 'a' and comp[0] != comp[1]:
            v0 = v0 / 2; g1 = [x / 2 for x in g1]          # the independent variable h_ab (a < b) enters as h_ab + h_ba
        outv.append((v0, g1))
    return outv


def base_from_spec(spec):
    h = [[ZERO for _ in idx] for _ in idx]; a = [ZERO for _ in idx]
    for k, mon, c in spec:
        term = C(c) * ypow(mon)
        if COMPS[k][0] == 'a':
            a[COMPS[k][1]] = a[COMPS[k][1]] + term
        else:
            i, j = COMPS[k]
            h[i][j] = h[i][j] + term
            if i != j:
                h[j][i] = h[j][i] + term
    return h, a


def pure_gauge(bg, eta, eta5):
    """delta h_ab = D_a eta_b + D_b eta_a,  delta a_n = eta^m F_mn + d_n eta_5   ((2.10a), (2.10b') at lowest order)"""
    T = bg.T
    h = [[T(d(eta[j], i) + d(eta[i], j) - psum(bg.Gb[c][i][j] * eta[c] for c in idx) * C(2), 3) for j in idx] for i in idx]
    eup = [T(psum(bg.gbinv[m][b] * eta[b] for b in idx), 3) for m in idx]
    a = [T(psum(eup[m] * bg.Fbg[m][n] for m in idx) + d(eta5, n), 3) for n in idx]
    return h, a


_W = {}
WAYS = CL.WAYS


def _init(Rd, F, Q, cubA):
    _W['args'] = (Rd, F, Q, cubA); _W['bg'] = {}; _W['probes'] = probe_fields()


def _task(job):
    tag, way, spec = job
    t0 = time.time()
    if way not in _W['bg']:
        _W['bg'][way] = Background(*_W['args'], *WAYS[way])
    bg = _W['bg'][way]
    if tag[0] == 'sym':                               # a pure-gauge base field on the UNGAUGED operator
        eta = [ZERO] * 4; eta5 = ZERO
        k, mon, c = spec
        if k == 4:
            eta5 = C(c) * ypow(mon)
        else:
            eta[k] = C(c) * ypow(mon)
        r = el(bg, pure_gauge(bg, eta, eta5), _W['probes'], gauge=False)
    else:
        r = el(bg, base_from_spec(spec), _W['probes'])
    return tag, way, r, time.time() - t0


CROSS = [('lin', 10, 0), ('rad', 12), ('lin', 1, 2), ('rad', 0), ('lin', 7, 1), ('lin', 13, 3)]


def task_list(mode, cols, method, ncross, nsym):
    """'all': every base field (14 x 15 = 210, + 14 radial); 'min': what E, omega, Omega need (84); 'cols': a smoke -- first-order
    data only for the columns `cols`.  The first `ncross` of CROSS also by the other method; `nsym` pure-gauge fields (symmetry)."""
    J = []; spec_of = {}
    for k in range(NC):
        spec_of[('const', k)] = [(k, (0, 0, 0, 0), 1)]
        J.append((('const', k), 'plain'))
        for m in idx:
            spec_of[('lin', k, m)] = [(k, mono(m), 1)]
            J.append((('lin', k, m), method if (mode != 'cols' or k in cols) else 'plain'))
        spec_of[('rad', k)] = [(k, mono(m, m), HALF) for m in idx]
        J.append((('rad', k), method if (mode != 'cols' or k in cols) else 'plain'))
    if mode == 'all':
        for k in range(NC):
            for m in idx:
                for n in idx:
                    if m <= n:
                        spec_of[('quad', k, m, n)] = [(k, mono(m, n), HALF if m == n else 1)]
                        J.append((('quad', k, m, n), 'plain'))
    other = 'translate' if method == 'expand' else 'expand'
    for tag in CROSS[:ncross]:
        if mode != 'cols' or tag[1] in cols:
            J.append((tag, other))
    sym = []
    for k in range(5):                                # eta_0..3 and eta_5: linear and quadratic (eta(0) = 0: only the field equations AT the origin are needed)
        for m in idx:
            sym.append((('sym', k, mono(m)), (k, mono(m), 1)))
        for m in idx:
            for n in idx:
                if m <= n:
                    sym.append((('sym', k, mono(m, n)), (k, mono(m, n), 1)))
    for tag, spec in sym[:nsym] if nsym < len(sym) else sym:
        spec_of[tag] = spec; J.append((tag, 'plain'))
    return [(tag, way, spec_of[tag]) for tag, way in J]


def bgkey(Rd, F, Q, cubA):
    txt = repr((CODE_VERSION, sorted((k, str(x)) for k, x in Rd.items()), [str(x) for x in F], [sorted((k, str(x)) for k, x in qn.items()) for qn in Q],
                [sorted((k, str(x)) for k, x in cn.items()) for cn in cubA]))
    return hashlib.sha1(txt.encode()).hexdigest()[:16]


def derive(Rd, F, Q, cubA, mode, cols, method, ncross, nsym, jobs, label, t00, ck):
    J = task_list(mode, cols, method, ncross, nsym)
    bk = bgkey(Rd, F, Q, cubA)
    keyof = lambda tag, way: f"{bk}|{tag}|{way}"
    way_of = {}
    for tag, way, _ in J:
        way_of.setdefault(tag, way)
    cost = {'translate': 0, 'expand': 1, 'plain': 2}
    J.sort(key=lambda j: cost[j[1]])
    res = {}; times = {}; ntot = len(J); todo = []
    for job in J:
        k = keyof(job[0], job[1])
        if k in ck.done:
            res[(job[0], job[1])] = ck.done[k][0]; times.setdefault(job[1], []).append(ck.done[k][1])
        else:
            todo.append(job)
    if len(todo) < ntot:
        say(f"   {label}: {ntot - len(todo)} of {ntot} base fields read from the checkpoint, {len(todo)} to compute")

    def take(tag, way, r, dt):
        res[(tag, way)] = r; times.setdefault(way, []).append(dt)
        ck.put(keyof(tag, way), r, dt)
        if len(res) % 20 == 0 or len(res) == ntot:
            say(f"   {label}: {len(res)}/{ntot} base fields  [{time.time() - t00:.0f} s]")
    if todo and jobs > 1:
        import multiprocessing as mp
        with mp.Pool(jobs, initializer=_init, initargs=(Rd, F, Q, cubA)) as pool:
            for tag, way, r, dt in pool.imap_unordered(_task, todo, chunksize=1):
                take(tag, way, r, dt)
    elif todo:
        _init(Rd, F, Q, cubA)
        for job in todo:
            take(*_task(job))
    tline = "; ".join(f"{k}: {len(v_)} tasks, mean {sum(v_) / len(v_):.1f} s, sum {sum(v_):.0f} s" for k, v_ in times.items())
    val0 = lambda tag: sp.Matrix([x[0] for x in res[(tag, way_of[tag])]])
    C0 = sp.Matrix.hstack(*[val0(('const', k)) for k in range(NC)])
    B0 = [sp.Matrix.hstack(*[val0(('lin', k, m)) for k in range(NC)]) for m in idx]
    fo = lambda tag, way: CL.first_order(tag, way, res[(tag, way)], C0, B0)
    z = [sp.zeros(NC, 1)] * 4
    grad = lambda tag: fo(tag, way_of[tag])[1] if way_of[tag] != 'plain' else z
    B1 = [[sp.Matrix.hstack(*[grad(('lin', k, m))[c] for k in range(NC)]) for c in idx] for m in idx]       # B1[m][c] = d_c B^m
    S0 = sp.Matrix.hstack(*[val0(('rad', k)) for k in range(NC)]) / 4
    S1 = [sp.Matrix.hstack(*[grad(('rad', k))[c] for k in range(NC)]) / 4 for c in idx]
    laplace = None
    if mode == 'all':
        laplace = all(val0(('quad', k, m, n)) == (S0[:, k] if m == n else sp.zeros(NC, 1)) for k in range(NC) for m in idx for n in idx if m <= n)
    other = 'translate' if method == 'expand' else 'expand'
    cross = {}
    for tag in CROSS[:ncross]:
        if way_of.get(tag) == method and (tag, other) in res:
            cross[tag] = fo(tag, method) == fo(tag, other)
    symres = {tag: [x[0] for x in res[(tag, 'plain')]] for tag in way_of if tag[0] == 'sym'}
    return dict(C0=C0, B0=B0, B1=B1, S0=S0, S1=S1, laplace=laplace, tline=tline, cross=cross, sym=symres)


def operator_data(dd, gam1, mutate='none'):
    """omega, E, Tr Omega^2 at the origin from the measured coefficient fields (class710_operator_derive.operator_data, 14 x 14)"""
    S0, S1, B0, B1, C0 = dd['S0'], dd['S1'], dd['B0'], dd['B1'], dd['C0']
    Si = S0.inv()
    f = HALF if mutate != 'halfw' else sp.Rational(1, 4)
    w0 = [Si * B0[m] * f for m in idx]
    w1 = [[(Si * B1[m][c] - Si * S1[c] * Si * B0[m]) * f + gam1[c][m] * sp.eye(NC) / 2 for m in idx] for c in idx]       # w1[c][m] = d_c w^m
    c0 = Si * C0
    E = c0 - sum((w1[m][m] for m in idx), sp.zeros(NC, NC)) - sum((w0[m] * w0[m] for m in idx), sp.zeros(NC, NC))
    Om = {(m, n): w1[m][n] - w1[n][m] + w0[m] * w0[n] - w0[n] * w0[m] for m in idx for n in idx}
    Om2 = sp.expand(sum((Om[k] * Om[k]).trace() for k in Om))
    return w0, E.applyfunc(sp.expand), Om2


def hcoords(h):
    return [h[i, j] for i, j in PAIRS]


def hbasis(p):
    i, j = PAIRS[p]; M = sp.zeros(4, 4); M[i, j] = 1; M[j, i] = 1
    return M


def predicted_B(F, sigma=1):
    """the first-derivative coefficient B^r at the origin from the cross term (header, A1)
        X = h^{am} F_m^b (d_a a_b - d_b a_a) - 1/2 h F^{rs} d_r a_s :
        (B^r)_{(ab), s} = dX/dh_ab at d_r a_s = 1;      (B^r)_{s, (ab)} = -(d/dy_r) dX/d(d_rho a_s) on h = E_ab y_r, summed over rho"""
    B = [sp.zeros(NC, NC) for _ in idx]
    for r in idx:
        for s in idx:
            M = sp.zeros(4, 4)
            for i in idx:
                for j in idx:
                    M[i, j] = (sp.Integer(i == r) * F[j, s] + sp.Integer(j == r) * F[i, s]) / 2 - (sp.Integer(i == s) * F[j, r] + sp.Integer(j == s) * F[i, r]) / 2 - sp.Integer(i == j) * F[r, s] / 2
            for p in range(10):
                B[r][p, 10 + s] = sigma * M[PAIRS[p][0], PAIRS[p][1]]
        for p in range(10):
            Eb = hbasis(p)
            for s in idx:
                B[r][10 + s, p] = -sigma * (sum(Eb[r, m] * F[m, s] for m in idx) - sum(Eb[s, m] * F[m, r] for m in idx) - Eb.trace() * F[r, s] / 2)
    return B


def photon_ghost_row(bg, F, DdivF):
    """C_5 = g^{mn} D_m (delta a_n) at the origin on every polynomial eta of degree <= 2 (eta_m and eta_5): the scalar ghost's row"""
    T = bg.T

    def C5(eta, eta5):
        _, a = pure_gauge(bg, eta, eta5)
        e = T(psum(bg.gbinv[m][n] * (d(a[n], m) - psum(bg.Gb[l][m][n] * a[l] for l in idx)) for m in idx for n in idx), 0)
        return lin_expand(e, 0)[0]
    ok = True; notes = []
    for k in idx:                                      # on eta_k
        e = lambda mon, c=1: ([C(c) * ypow(mon) if i == k else ZERO for i in idx], ZERO)
        c_const = C5(*e((0, 0, 0, 0)))
        ok &= c_const == DdivF[k]                      # D^m F_km
        for r in idx:
            ok &= C5(*e(mono(r))) == F[k, r]           # F_km D^m on eta_k = y_r
        for m in idx:
            for n in idx:
                if m <= n:
                    ok &= C5(*e(mono(m, n))) == 0
    lap = True
    for mon, want in [((0, 0, 0, 0), 0)] + [(mono(r), 0) for r in idx] + [(mono(m, n), sp.Integer(2 if m == n else 0)) for m in idx for n in idx if m <= n]:
        lap &= C5([ZERO] * 4, ypow(mon)) == want       # Box eta_5: unit Laplacian, nothing else
    return ok, lap


def structures(R, F, DF):
    """the ten pointwise structures of the header (A4) and the derived ones"""
    Ric = sp.Matrix(4, 4, lambda b, e: sum(R[(i, b, i, e)] for i in idx)); Rs = sp.expand(Ric.trace())
    F2 = F * F.T                                       # F2_mn = F_ma F_na
    Fsq = sp.expand(F2.trace())                        # F^2 = F_mn F_mn
    st = {'Ric2': sp.expand((Ric * Ric).trace()), 'R2': Rs**2, 'Riem2': sp.expand(sum(x**2 for x in R.values())),
          'F4': sp.expand((F2 * F2).trace()), 'F22': Fsq**2, 'RF2': sp.expand(Rs * Fsq), 'RicF2': sp.expand((Ric * F2).trace()),
          'RiemFF': sp.expand(sum(R[(i, j, k, l)] * F[i, j] * F[k, l] for i, j, k, l in itertools.product(idx, repeat=4))),
          'DaF2': sp.expand(sum(DF[(r, m, n)]**2 for r, m, n in itertools.product(idx, repeat=3))),
          'DF2': sp.expand(sum(sum(DF[(m, m, n)] for m in idx)**2 for n in idx))}
    st['RicT'] = sp.expand(st['RicF2'] - st['RF2'] / 4); st['T2'] = sp.expand(st['F4'] - st['F22'] / 4)
    st['GB'] = sp.expand(st['Riem2'] - 4 * st['Ric2'] + st['R2'])
    T = F2 - Fsq * sp.eye(4) / 4
    return Ric, Rs, T, st


STKEYS = ['Ric2', 'R2', 'Riem2', 'F4', 'F22', 'RF2', 'RicF2', 'RiemFF', 'DaF2', 'DF2']
TARGET = {'Ric2': sp.Rational(9, 20), 'R2': sp.Rational(-1, 40), 'RicT': sp.Rational(-1, 6), 'T2': sp.Rational(13, 24), 'DF2': sp.Rational(1, 6)}   # (4.16a), standard sign of R
NG_F = {'RicT': sp.Rational(-1, 6), 'T2': sp.Rational(13, 24), 'DF2': sp.Rational(1, 6)}                                                        # the F part of (4.13)


def sample(rng, kind):
    R = D.riemann_sample(rng)
    if kind == 'F = 0':
        F = sp.zeros(4, 4); Q = [{(m, r): sp.Integer(0) for m in idx for r in idx if m <= r} for _ in idx]
        cubA = [{} for _ in idx]
        return R, F, Q, cubA
    F = sp.zeros(4, 4)
    for m in idx:
        for n in idx:
            if m < n:
                F[m, n] = sp.Rational(rng.choice((-3, -2, -1, 1, 2, 3)), rng.randint(1, 3)); F[n, m] = -F[m, n]
    Q = [{(m, r): sp.Rational(rng.randint(-3, 3), rng.randint(1, 3)) for m in idx for r in idx if m <= r} for _ in idx]
    cubA = [{mm: sp.Rational(rng.randint(-2, 2), rng.randint(1, 2)) for mm in itertools.product(range(4), repeat=4) if sum(mm) == 3} for _ in idx]
    if kind == 'on shell':
        qf = lambda n, m, r: Q[n][(min(m, r), max(m, r))]
        for n in idx:                                  # D^m F_mn = sum_m (q_n,mm - q_m,nm) = 0, by shifting q_n,kk for one k != n
            k = (n + 1) % 4
            div = sum(qf(n, m, m) - qf(m, n, m) for m in idx)
            Q[n][(k, k)] = Q[n][(k, k)] - div
        F2 = F * F.T
        R = D.with_ricci(R, (F2 - F2.trace() * sp.eye(4) / 4) / 2)
    return R, F, Q, cubA


def DF_of(Q):
    qf = lambda n, m, r: Q[n][(min(m, r), max(m, r))]
    return {(r, m, n): qf(n, m, r) - qf(m, n, r) for r, m, n in itertools.product(idx, repeat=3)}          # D_r F_mn


def run_background(kind, i, R, F, Q, cubA, args, first, t00, ck):
    DF = DF_of(Q)
    Ric, Rs, Tm, st = structures(R, F, DF)
    divF = [sum(DF[(m, k, m)] for m in idx) for k in idx]                 # D^m F_km
    onsh = (Ric - Tm / 2) == sp.zeros(4, 4) and all(x == 0 for x in divF)
    label = f"{kind} {i}"
    say(f"\n## {label}: R = {Rs}, F^2 = {sp.expand((F * F.T).trace())}, (D.F)^2 = {st['DF2']}, (D_a F_bc)^2 = {st['DaF2']}; field equations hold at 0: {onsh}")
    say("   structures: " + ", ".join(f"{k} = {st[k]}" for k in STKEYS))
    cols = [int(x) for x in args.cols.split(',')] if args.mode == 'cols' else list(range(NC))
    nsym = args.symmetry if (onsh and first.get('sym', True)) else 0
    if nsym:
        first['sym'] = False
    ncross = args.crosscheck if first.get('cross', True) and kind != 'F = 0' else 0
    if ncross:
        first['cross'] = False
    dd = derive(R, F, Q, cubA, args.mode, cols, args.method, ncross, nsym, args.jobs, label, t00, ck)
    say(f"   timing: {dd['tline']}")
    flags = {}
    if dd['cross']:
        say(f"   cross-check, first-order data by 'translate' and by 'expand' on {len(dd['cross'])} base fields: {'identical' if all(dd['cross'].values()) else 'DIFFER: ' + str(dd['cross'])}")
        flags['cross'] = all(dd['cross'].values())
    if dd['sym']:
        bad = [tag for tag, v in dd['sym'].items() if any(x != 0 for x in v)]
        say(f"   A5a gauge symmetry of the UNGAUGED operator on {len(dd['sym'])} pure-gauge fields (eta_m, eta_5; linear and quadratic), on shell: {'all 14 components vanish on every one' if not bad else 'FAILS on ' + str(bad[:6])}")
        flags['A5a'] = not bad
    bg1 = Background(R, F, Q, cubA, 0, 1)
    gam1 = bg1.gam1()
    w0, E, Om2 = operator_data(dd, gam1)
    P10 = D.op_matrix10(lambda h: h / 2 - sp.eye(4) * h.trace() / 4)
    Sw = sp.diag(P10, sp.eye(4))
    hh = lambda M: M[:10, :10]; aa = lambda M: M[10:, 10:]
    # A1
    a1s = dd['S0'] == Sw and all(dd['S1'][c] == sp.zeros(NC, NC) for c in idx)
    diag_zero = all(hh(dd['B0'][m]) == sp.zeros(10, 10) and aa(dd['B0'][m]) == sp.zeros(4, 4) for m in idx)
    sig = [s for s in (1, -1) if all(dd['B0'][m] == predicted_B(F, s)[m] for m in idx)]
    nfirst = sum(1 for m in idx for x in dd['B0'][m] if x != 0)
    Wd = sp.diag(*[1 if (c[0] == 'a' or c[0] == c[1]) else 2 for c in COMPS])     # h_ab (a < b) is one variable standing for two entries: the bilinear form is Wd S
    antis = all((Wd * Sw * w0[m]).T == -(Wd * Sw * w0[m]) for m in idx)
    say(f"   A1  S = diag(P, 1), constant to first order: {a1s}" + (f"; Laplace type on all {NC * 10} quadratic base fields: {dd['laplace']}" if dd['laplace'] is not None else "")
        + f";  first-derivative terms in the mixed entries only: {diag_zero} ({nfirst} nonzero entries);  = the hand cross term with overall sign {sig if sig else 'NEITHER SIGN'};  omega antisymmetric with respect to S (as a bilinear form on the 14 variables): {antis}")
    flags['A1'] = bool(a1s and diag_zero and (sig or kind == 'F = 0') and antis and dd['laplace'] is not False)
    if args.mode == 'cols':
        say("   (smoke, --mode cols: the rest needs every column)")
        return flags, None
    # A0 / A2: the diagonal blocks against steps 2 and 1 at F = 0
    E0, Om2_h0 = GM.build_E(R, Ric, Rs, sp.zeros(4, 1), sp.zeros(4, 4))
    Eh0 = E0[:10, :10].applyfunc(sp.expand)
    if kind == 'F = 0':
        a0 = hh(E) == Eh0 and aa(E) == (-Ric).applyfunc(sp.expand) and E[:10, 10:] == sp.zeros(10, 4) and E[10:, :10] == sp.zeros(4, 10) and nfirst == 0
        say(f"   A0  F = 0: h block = step 2's E, a block = -Ric, no mixing, no first-derivative terms: {a0};  Tr Omega^2 = step 2's + the vector's (-Riem^2): {Om2 == sp.expand(Om2_h0 - st['Riem2'])}")
        flags['A0'] = bool(a0 and Om2 == sp.expand(Om2_h0 - st['Riem2']))
    else:
        mixed_from_DF = None
        say(f"   A2  h block - step 2's E: {'zero' if hh(E) == Eh0 else 'nonzero'} (expected nonzero, quadratic in F);  a block + Ric: {'zero' if aa(E) == (-Ric).applyfunc(sp.expand) else 'nonzero'};  mixed blocks: {'zero' if E[:10, 10:] == sp.zeros(10, 4) and E[10:, :10] == sp.zeros(4, 10) else 'nonzero'} (expected linear in DF)")
    # the traces
    X = E + Rs / 6 * sp.eye(NC)
    trX2 = sp.expand((X * X).trace())
    say(f"   A3  Tr Omega^2 = {Om2};  Tr (E + R/6)^2 = {trX2};  Tr E = {sp.expand(E.trace())}")
    # A4: the non-ghost part, pointwise
    ng = GM.hv340(NC, Om2, E, Rs, st['Ric2'])
    ng0 = sp.expand(GM.hv340(10, Om2_h0, Eh0, Rs, st['Ric2']) + GM.hv340(4, -st['Riem2'], -Ric, Rs, st['Ric2']))
    ng_pred = sp.expand(ng0 + sum(NG_F[k] * st[k] for k in NG_F))
    a4 = ng == ng_pred
    say(f"   A4  (3.40) on the 14 fields, pointwise: {ng};  steps 1-2 at F = 0 ({ng0}) - 1/6 Ric.T + 13/24 T^2 + 1/6 (D.F)^2 = {ng_pred};  equal: {a4}" + ("" if a4 else f"   difference: {sp.expand(ng - ng_pred)}"))
    flags['A4'] = a4
    ctl = GM.hv340(NC, *operator_data(dd, gam1, 'halfw')[:0:-1], Rs, st['Ric2']) != ng_pred if kind != 'F = 0' else None
    # A5: the ghosts
    s_gh, second_ok, first_zero, E_gh, _ = GH.ghost_operator(D.Background(R, sp.zeros(4, 1), sp.zeros(4, 4), {m: sp.Integer(0) for m in itertools.product(range(4), repeat=4) if sum(m) == 3}), 0, 1)
    vec_ok = s_gh == 1 and second_ok and first_zero and E_gh is not None and E_gh.applyfunc(sp.expand) == Ric.applyfunc(sp.expand)
    corner, lap = photon_ghost_row(Background(R, F, Q, cubA, 0, 0), F, divF)
    say(f"   A5  vector ghost: unit Laplacian, no first-derivative terms, E_gh = Ric: {vec_ok};  scalar ghost: unit Laplacian and nothing else: {lap};  the one mixed corner = F_km D^m + D^m F_km on eta_k: {corner}")
    flags['A5'] = bool(vec_ok and lap and corner)
    gvec = GM.hv340(4, -st['Riem2'], Ric, Rs, st['Ric2']); gsc = GM.hv340(1, 0, sp.zeros(1, 1), Rs, st['Ric2'])
    ghost = sp.expand(-2 * (gvec + gsc))
    tot = sp.expand(ng + ghost)
    # A6-A8
    gcoef = None
    f0 = sp.expand(ng0 + ghost)                       # the F = 0 total: must be 9/20 Ric^2 - 1/40 R^2 + g GB
    rem0 = sp.expand(f0 - TARGET['Ric2'] * st['Ric2'] - TARGET['R2'] * st['R2'])
    gcoef = sp.cancel(rem0 / st['GB']) if st['GB'] != 0 else None
    pred = sp.expand(f0 + sum(NG_F[k] * st[k] for k in NG_F))
    a7 = tot == pred
    say(f"   A6  ghost contribution (triangular: the diagonal blocks): {ghost}")
    say(f"   A7  total, pointwise: {tot};  (4.16a) + g Gauss-Bonnet with g = {gcoef}: equal: {a7}" + ("" if a7 else f"   difference: {sp.expand(tot - pred)}"))
    flags['A7'] = a7
    if ctl is not None:
        c2 = sp.expand(ng - 2 * (GM.hv340(4, -st['Riem2'], -Ric, Rs, st['Ric2']) + gsc)) != pred
        c3 = sp.expand(ng - 2 * gvec) != pred
        say(f"   controls (must break A4/A7): omega halved: {'breaks it' if ctl else 'DOES NOT'};  the vector ghost's Ric with the other sign: {'breaks it' if c2 else 'DOES NOT'};  the scalar ghost left out: {'breaks it' if c3 else 'DOES NOT'}")
        flags['controls'] = bool(ctl and c2 and c3)
    if onsh and gcoef is not None:
        rem = sp.expand(tot - gcoef * st['GB'])
        a8 = rem == sp.Rational(137, 60) * st['Ric2']
        say(f"   A8  on shell: total - ({gcoef}) Gauss-Bonnet = {rem};  (137/60) Ric^2 = {sp.Rational(137, 60) * st['Ric2']};  equal: {a8}" + (f"   ratio to Ric^2: {sp.cancel(rem / st['Ric2'])}" if not a8 and st['Ric2'] != 0 else ""))
        flags['A8'] = a8
    return flags, (gcoef, tot, st)


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261114)
    ap.add_argument('--onshell', type=int, default=2)
    ap.add_argument('--offshell', type=int, default=2)
    ap.add_argument('--f0', type=int, default=1, help='backgrounds with F = 0 (the validation A0)')
    ap.add_argument('--mode', choices=('all', 'min', 'cols'), default='all')
    ap.add_argument('--cols', default='10,0', help='with --mode cols: the columns that get first-order data (0-9 the h pairs, 10-13 the potential)')
    ap.add_argument('--method', choices=('expand', 'translate'), default='expand')
    ap.add_argument('--crosscheck', type=int, default=6, help='on the first background with F != 0, this many base fields by the other method too (0-6)')
    ap.add_argument('--symmetry', type=int, default=70, help='on the first on-shell background, this many pure-gauge fields on the ungauged operator (0-70)')
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    t0 = time.time()
    if args.resume and not args.log:
        print("--resume needs --log (the checkpoint is <log>.ckpt)"); sys.exit(1)
    if args.log and os.path.exists(args.log) and os.path.getsize(args.log) > 0 and not args.resume:
        print(f"refusing to overwrite the non-empty log {args.log}; add --resume to continue that run, or to start over:  del {args.log.replace('/', chr(92))}  and  del {args.log.replace('/', chr(92))}.ckpt"); sys.exit(1)
    if args.log and os.path.dirname(args.log):
        os.makedirs(os.path.dirname(args.log), exist_ok=True)
    ck = CL.Checkpoint(args.log + '.ckpt' if args.log else None, args.resume)
    if args.log:
        out = open(args.log, 'a' if args.resume else 'w')
        if args.resume:
            say(f"\n# RESUMED {time.strftime('%Y-%m-%d %H:%M:%S')}: {len(ck.done)} base fields in the checkpoint; everything below is recomputed from it or run afresh")
    say(f"# einstein_maxwell_operator_derive (stage A, the paper's gauge (3.1)): seed {args.seed}, {args.onshell} on-shell + {args.offshell} off-shell + {args.f0} F = 0 backgrounds, mode {args.mode}"
        + (f" (columns {args.cols})" if args.mode == 'cols' else "") + f", method {args.method}, crosscheck {args.crosscheck}, symmetry {args.symmetry}, jobs {args.jobs}")
    say(f"# code version {CODE_VERSION} (hash of this script and its four imports; part of every checkpoint key)")
    say(f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}, python-flint {flint.__version__}; {sys.executable}")
    rng = random.Random(args.seed)
    allflags = []; first = {}; rows = []
    for kind, n in (('on shell', args.onshell), ('off shell', args.offshell), ('F = 0', args.f0)):
        for i in range(n):
            R, F, Q, cubA = sample(rng, kind)
            flags, row = run_background(kind, i, R, F, Q, cubA, args, first, t0, ck)
            allflags.append((kind, flags))
            if row is not None:
                rows.append((kind, row))
            say(f"   -> {kind} {i}: " + ", ".join(f"{k} {'ok' if v else 'FAILS'}" for k, v in flags.items()) + f"   [{time.time() - t0:.0f} s]")
    say("")
    ok = all(all(f.values()) for _, f in allflags)
    gs = {r[0] for _, r in rows if r[0] is not None}
    if rows:
        a9 = len(gs) == 1
        say(f"A9  the Gauss-Bonnet coefficient g on every background: {sorted(gs, key=str)}: {'one value, independent of F' if a9 else 'NOT ONE VALUE'}")
        ok &= a9
    kinds = {k for k, _ in allflags}
    complete = args.mode != 'cols' and {'on shell', 'off shell', 'F = 0'} <= kinds and any('A5a' in f for _, f in allflags) and any('cross' in f for _, f in allflags)
    say("status=" + ("FAIL" if not ok else "PASS" if complete else "INCOMPLETE (PASS needs --mode all or min, an on-shell, an off-shell and an F = 0 background, the symmetry test and the cross-check)"))
    say(f"END einstein_maxwell_operator_derive {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
