#!/usr/bin/env python3
"""
Rung 3, step 3b: 't Hooft-Veltman section 7 -- ROUTE 2, the paper's own: the class (7.10), computed directly.
PREDICTIONS WRITTEN 13 Nov 2026 11:59 ET (clock read by the command that wrote them), before any code existed; shown to the
author, who said yes to them at 12:12 ET.  They stand below as first written (rule 11); corrections follow them.

WHAT THE PRINT SAYS (pp. 88-90, read as images 13 Nov; (7.8)-(7.9) also read by the author, 13 Nov).
  (7.1) L = sqrt(g)(-R - 1/2 (d phi)^2 + a R phi^2 + b R^{mu nu} d phi d phi), b not considered.  (7.8) is (7.1) at b = 0 written
  with f = 1 - a phi^2, plus a mass term -1/2 m^2 phi^2 that (7.1) does not have.  The metric alone is rescaled (7.7), g -> g/f;
  the scalar is NOT redefined.  The result (7.9) is filed under the class (7.10), sqrt(g){-R + 1/2 (d phi)^2 f_1(phi) + f_2(phi)}.
  The gauge is not stated ("not too difficult to see the changes with respect to the treatment of the previous sections").
  On shell or off shell is not stated.  NO RESULT IS PRINTED: only that the counter-Lagrangian is a power series in the background
  field with non-trivial coefficients, that a = 1/12 "is of little help", and that it "seems not to be of any simple form".
  Misprint, probable: (7.9) prints the kinetic term as -(2/2f)(d phi)^2 (two independent readings).  Rescaling (7.8) gives -(1/2f);
  weyl_map_check.py (the author's PASS run) has K = 1/f + 3 f'^2/f^2, i.e. 1/2f.  This script takes K as input and cannot test it.

SET-UP (standard sign, as in step 3b; R_HV = -R).
  L = sqrt(g)( R - 1/2 K(phi) (d phi)^2 ),  m = 0, f_2 = 0.  K is ANY function: the code carries its jets K0, K1, K2 at the
  background value as free rationals.  The theory (7.1) is the special case
        K_a = (1 + (12 a^2 - a) phi^2) / (1 - a phi^2)^2          [= 1/F^2, F = 1 - phi^2/12, at a = 1/12;  = 1 at a = 0],
  so general a is covered by general K, and named values of a are run as cases (0, 1/12, 1/6, 1, one random rational).
  Background: metric in normal coordinates on a random algebraic curvature tensor, phi-bg(x) a polynomial with rational jets.
  Notation at the origin:  u = d phi-bg,  Wt = DD phi-bg + 1/2 k1 u u,  k1 = K1/K0,  beta = tr Wt = Box phi-bg + 1/2 k1 u^2
  (beta = 0 is the scalar field equation of this Lagrangian),  delta = 1/2 k1 beta + 1/4 k1^2 u^2.
  Gauge: -1/2 C^2,  C_mu = D^nu h_mu nu - 1/2 D_mu h - c G(phi-bg) phi d_mu phi-bg, with
        G = K   "image gauge": (4.27) of the canonical-field theory rewritten in phi.  The gauge route 1 assumed by argument.
        G = 1   "literal gauge": (4.27) as printed, untouched.  [run only with the author's yes -- see R8]
  Operator form: -(S D^2 + A^mu D_mu + B) from the Lagrangian alone, by jet probes with the background translated to first order
  (the extension proven in conjugation_proving.py); then multiplied from the left by S^-1 -- the paper's own P^-1 replacement
  (5.3)-(5.4), extended to the scalar entry; then omega from S^-1 A, E, Omega, and (3.40).

PREDICTIONS (rule 25).  "Step 3" means the operators derived in gravscalar_operator_derive.py and ghost_operator_derive.py.
  R0 [validation of the new evaluator]  K0 = 1, K1 = K2 = 0: step 3b's operator and ghost exactly; no first-derivative terms.
  R1 [image gauge, c = 1]  S = diag(P, K0) delta^{mu nu}.  First-derivative terms in ONE entry only, scalar-scalar:
        S^-1 A^mu = k1 u^mu there, so omega_mu = 1/2 k1 u_mu on the scalar, zero mixed, the metric's connection on h.
  R2  Tr Omega^2 = step 3's (the scalar's omega is a gradient).
  R3  E, as an 11 x 11 matrix, on shell AND off shell, with v v -> K0 u u in step 3's h-h block:
        (E phi)_ab = -2 K0 Wt_ab phi,     (E h)_phi = 1/2 beta tr h - Wt^ab h_ab,     E_phi-phi = -K0 u^2 + delta.
      Equivalently E = J^-1 [ E_step3(Riem, v = sqrt(K0) u, W = sqrt(K0) Wt) + delta e_phi-phi ] J,  J = diag(1, sqrt(K0)).
      Everything above is rational in K0, k1; no square root appears in any trace.
  R4  K2 enters nothing at the origin (mutating it changes no output).
  R5 [ghost]  from C under delta h = D eta + D eta, delta phi = eta.d phi-bg:  no first-derivative terms,
        E_gh = Ric - c K0 u u,  no Wt term, no k1.
  R6 [the divergence, pointwise, before the reductions]  (3.40) on (h, phi) minus twice the ghost
        = step 3's value at (Riem, sqrt(K0) u, sqrt(K0) Wt)  +  1/2 delta (R/6 - K0 u^2)  +  1/4 delta^2.
      After the two reductions, off shell, image gauge:
        (1/80) R^2 + (43/120) Ric^2 + 1/2 K0^2 u^4 + (1/12) K0 R u^2 + 1/2 K0 beta^2
              + 1/4 k1 beta (R/6 - K0 u^2) + (1/16) k1^2 beta^2 .
      [caution: the reductions integrate by parts; with phi-dependent coefficients that produces extra terms.  The pointwise
       statement is the prediction; the reduced line is what I expect and will be re-derived, not assumed, when coded.]
  R7 [THE NUMBER]  on shell (Ric = 1/2 K0 u u, R = 1/2 K0 u^2, beta = 0):
        (203/320) K0^2 (u^2)^2       for every K, hence for every a;   = (203/320) K^2 sqrt(g) ((d phi)^2)^2.
      The same expression holds in the original metric of (7.1): sqrt(g)((d phi)^2)^2 is Weyl invariant in four dimensions.
      So on shell a = 1/12 is not special at m = 0.  What a = 1/12 does do: k1 = phi/(3F) and delta = (phi/(6F)) beta are
      rational in phi; for other a the canonical field needs sqrt(1 + (12a^2 - a) phi^2).  OFF shell the k1 terms of R6 are a
      power series in phi with coefficients depending on a -- my candidate for the paper's "non-trivial coefficients", in this
      gauge only; the paper's gauge is unknown, so this cannot be checked against it.
  R8 [a question, needs the author's yes: literal gauge G = 1]  Mixed first-derivative terms appear, proportional to (K0 - 1)
      and to k1; omega is no longer a gradient; E and E_gh = Ric - c u u differ from R3, R5.  Prediction: ON SHELL the total is
      still (203/320) K0^2 u^4.  This IS a test of gauge independence on shell (the image gauge is not: there the operator is
      step 3's by construction).  It also covers the other gauge a 1974 reader might have meant.
  Controls (must FAIL): the scalar's S^-1 left out; omega halved; c = 0 at R1.   [as coded: the first two are applied inside every run, at no cost, and must break R3; c = 0 is a separate run]

WHAT A PASS WOULD AND WOULD NOT SHOW.
  Would: the gauge term and the ghost of the non-canonical theory carried in code (route 1 carried them by argument); a divergence
  computed for the section 7 theory; the translated-probe machinery on a gauge-fixed 11-component operator with a first-derivative
  term; with R8, gauge independence on shell between two gauges.
  Would not: anything about the gauge the paper used; the Faddeev-Popov -2, the two reductions and the on-shell substitution
  (shared with step 3 and with route 1); the local Jacobian (p. 70, ledger row HV-12, open); the mass term; b != 0; the theory
  quantized in the ORIGINAL metric of (7.1) without the Weyl map (a different operator, non-minimal without a new gauge).
  In the image gauge R3-R7 can fail only through the code or through my algebra for delta -- the physics there is route 1's.

SECOND METHOD: FORM (4.3.1; nothing of the physics done in it yet).  Its own header is written before its code.  Plan:
  F0  K = 1 first: the covariant E of step 3 (note 4d's formula), Tr E, Tr E^2, and the five coefficients of (5.22), by index
      contraction -- FORM proven on the known case before it is trusted on the new one (the author's rule of 10 Nov).
  F1  the class (7.10): the gauge-fixed second-order Lagrangian as a covariant formula in K0, k1, u, Wt; then R3, R6, R7 as formulas.
  F2  comparison: the Python run evaluates FORM's written output on its own backgrounds; nothing is typed from one log to the other.
  Not settled until a sandbox trial: how much FORM derives from the action itself and how much is fed to it as expansion formulas.

COST.  No forecast yet.  A sandbox smoke on a subset comes first and the forecast is made from its per-task line
(sandbox single-core seconds x tasks / 3).  A guess only: several times route 1's 54 minutes on the author's machine.

AMENDMENTS AND CORRECTIONS (13 Nov, in the order made; this block written 13 Nov 12:42 ET, clock read by the command that inserted it).
  b = 0 (author's question, 12:05-12:12 ET; hand argument, nothing run).  The paper's reason holds: the b term has four derivatives,
      the order of the one-loop divergence itself, so its loop effects start at six and cannot cancel (6.5); and the shift
      g -> g + b (d phi d phi - 1/2 g (d phi)^2) removes it at first order in b, leaving a multiple of b ((d phi)^2)^2 -- a scalar
      self-interaction, no change to any propagator, no spin-2 pole.  b = 0 as in the paper (the author's decision, 12:12 ET).
  CORRECTION TO delta (found by the first smoke, log file time 12:36 ET).  The set-up line above defines
      delta = 1/2 k1 beta + 1/4 k1^2 u^2.   WRONG: my slip in copying; the hand algebra gives 1/2 k1 Box phi-bg + 1/4 k1^2 u^2,
      which is 1/2 k1 beta, and I wrote beta where Box phi-bg stood.  The correct statement, which the run tests:
            delta = 1/2 k1 beta.
      The first wording was also inconsistent with R7 (it would not vanish on shell).  R6's reduced line and R7's remark on
      a = 1/12 were written with the right delta and stand.  The smoke's derived scalar-scalar entry was -K0 u^2 + 1/2 k1 beta.
  R8 IS A STATEMENT ABOUT THE INTEGRAL (seen while coding).  Two gauges give the same on-shell divergence only up to a total
      derivative, which need not vanish at a point.  In the image gauge the operator is step 3's and the comparison is pointwise;
      for the literal gauge a test needs several on-shell backgrounds and a fit modulo the total-derivative structures (with
      phi-dependent coefficients).  That fit is NOT in this script: --gauge literal runs and prints, and tests nothing about R8.
      So R8 costs more than "a second run"; it still waits for the author's yes.
  HOW THE FIRST-ORDER DATA ARE OBTAINED ('expand' beside 'translate').  The translated probes of conjugation_proving.py cost
      about 213 s per base field on one sandbox core here.  The default is a variant, 'expand': no translation, every truncation
      one order higher in y, the Euler-Lagrange expression of the base field kept to first order in y, and the base field's own
      y-dependence removed (first_order()).  About 32 s per base field.  It is the same information, but it is different code, so
      on the first background of every run --crosscheck base fields (default 6) are derived both ways and must agree exactly; the
      run fails otherwise.  --method translate runs everything the proven way.
  THE GHOST is not re-derived by new code: ghost_operator_derive.ghost_operator (the author's PASS run of 10 Nov) is called on
      the same background with the gauge condition's scalar coefficient at the origin, c K0 (image) or c (literal).  Exact,
      because that term carries no derivative of the fluctuation; but it means R5 tests one number passed to proven code.
  R6's COMPARISON TARGET is step 3's total as a formula on eight structures, fitted inside this script from
      hv1974_gravscalar_M's own functions on 12 random backgrounds (nothing typed), evaluated at the image background by
      homogeneity (v v -> K0 u u, W -> sqrt(K0) Wt enters quadratically).  R3's target is build_E of the same module, block by block.
      R7 is pointwise: the total minus the Gauss-Bonnet and Ricci-identity multiples fitted for step 3 equals (203/320) K0^2 (u^2)^2.

R8 MADE TESTABLE (the author's yes, 13 Nov 13:54 ET; this block written 13 Nov 13:55 ET, clock read by the command that inserted it,
before any of the --r8 code existed).
  The gauge family.  C_mu = D^nu h_mu nu - 1/2 D_mu h - c K(phi-bg) phi d_mu phi-bg, c any rational: the image gauge at c = 1.
  In the canonical field this is step 3's own family, C = ... - c psi d chi, for every K.  It is the test brief 21 suggested
  ("the step 3 scripts at --c != 1"), run on the section 7 theory.  The LITERAL gauge (no K) at K' != 0 is not in this family
  and stays untested: there the total derivatives carry unknown functions of K and one number cannot be isolated.
  The statement.  On shell at a point the total is an O(4) invariant of dimension four built from the Weyl tensor, u and the
  traceless Wt, and there are exactly four: Weyl^2, (u^2)^2, u.Wt.u, Wt.Wt.  Three combinations are total derivatives on shell:
        Gauss-Bonnet;     Wc.Wc - (tr Wc)^2 + Ric v v  (the Ricci identity; v, Wc the canonical field's);     v.Wc.v  (= 1/2 D_a(v^a v^2) on shell).
  So the on-shell total is  alpha (v^2)^2 + b v.Wc.v + C (Wc.Wc + Ric v v) + g GB  with ONE gauge-invariant number, alpha.
  R8a [can fail]  alpha = 203/320 at c = 0 and at c = 2, as at c = 1.  (alpha(c) is a polynomial of degree at most four in c;
        three values do not prove it constant, they make a non-constant one very unlikely.)
  R8b [a question]  g, the Gauss-Bonnet coefficient: I expect 7/12 for every c and do not know it.  b and C are reported.
  R8c  the four structures suffice: with five on-shell backgrounds the fifth has zero residual.
  How: --r8 with --onshell 5 or more on one case; the per-background checks R1-R3, R6, R7 are written for c = 1 and are
  skipped at other c (the ghost check R5 holds for any c and is kept).  The fit routine is first run, inside every --r8 run,
  on step 3's own formula at c = 1 (hv1974_gravscalar_M.build_E on six random on-shell backgrounds), where it must return
  alpha = 203/320, b = 0, C = 1, g = 7/12 -- the new machinery on the known case.
  Would not catch: an error in the family that is the same at every c; anything about the literal gauge.

WHAT EACH CHECK WOULD AND WOULD NOT CATCH.
  R1 the kinetic normalisation K0, K1 in S and in the first-derivative term, and stray mixed first-derivative terms (the gauge).
  R3 the whole of E: the gauge term's K, the conversion from (S, B, C) to E, the connection; it is the content.
  R2 that the scalar's omega is a gradient and the h-bundle's is the metric's.
  R6, R7 add (3.40) and the ghost to R3; given R2, R3, R5 they can fail only through the fitted step 3 formula or homogeneity.
  NOT caught by anything here: an error common to step 3 and this script in the Faddeev-Popov -2, the reductions, the on-shell
  substitution; the local Jacobian; K2 beyond the fact that it is drawn nonzero and absent from every prediction (R4).
  Orders y^3 and beyond of phi-bg and K's third derivative cannot reach the origin.

SMOKE RECORD (rule 11; sandbox, 2 cores, python 3.13.18, sympy 1.13.1, python-flint 0.9.0).
  1 (log file time 12:36 ET; free K jets (2, -1, 2), one off-shell background, --mode cols --cols 10 --crosscheck 2, 904 s):
     cross-check of 'expand' against 'translate' on two base fields: identical.  R1 held.  R3 FAILED on the scalar-scalar entry:
     derived -1645/144, predicted -265/24; the difference is exactly 1/4 k1^2 u^2 -- the slip in delta corrected above.
  2 (same background, --mode min, 1138 s; after the correction): R1, R2, R3, R5, R6 hold.
  3 (a = 1/12 at phi0 = 1/2, K0 = 2304/2209, one on-shell background, --mode min, 1431 s): R1, R2, R3, R5, R6 hold, and
     R7: total - (7/12) Gauss-Bonnet - (W.W - (tr W)^2 + Ric v v) = 216199872/24398405 = (203/320) K0^2 (u^2)^2.
  4 (the script as it stands, with the in-run controls; a = 1/6, K0 = 600/529, one off-shell background, --mode cols --cols 10,1,
     724 s): R1 holds, both controls break R3, R3 holds on the two columns tested.
  5 (--r8 --c 0, a = 1/12, the on-shell background of smoke 3, --mode min; killed after 13 base fields and restarted with
     --resume: the 13 were read from the checkpoint, a restart without --resume was refused; 1148 s after the resume): runs
     through; 60 nonzero first-derivative entries (4 at c = 1); the ghost is Ric; the pointwise total differs from c = 1's on
     the same background, as it may.  One background cannot test R8.  The R8 fit and its report were exercised with synthetic
     rows from step 3's formula (scratch harness, not committed): alpha = 203/320 recovered, and a total spoiled by v^4/7 FAILS.
  Per base field on one sandbox core: 'expand' 36-45 s, 'plain' 18-25 s, 'translate' 213 s.
  NOT run in the sandbox: --mode all (the 110 quadratic base fields), cases a = 0, 1, rand, --c 2, --gauge literal,
  --method translate as a whole, more than one background per run, and therefore no R8 fit on derived totals.

FRESH REVIEW (13 Nov afternoon, separate agent, no context; scratch copies, repository untouched; about 20 minutes).
No error found in the physics, in 'expand', or in the conversion to omega, E, Omega.  It did not make any complete real run.
Must-fix, all folded in:
  * the control "scalar S^-1 left out" is the identity at K0 = 1, so the main run would have ended status=FAIL on case a = 0;
    it is now skipped there.
  * the cross-check was weakest where needed: the first background was a = 0, where the first-derivative term and d S vanish
    and the 'rad' correction is untested.  The case order now starts with a = 1/12; the R8 runs now use --crosscheck 2.
  * the checkpoint key had no code version: after an edit, --resume would have reused old lines.  A hash of this file and its
    three imports is now part of the key.  A truncated last line no longer swallows the next record.
  * a run at c != 1 or in the literal gauge without --r8 printed status=PASS having tested no prediction; now INCOMPLETE.
What it ran: 'expand' against 'translate' on two base fields outside the CROSS list, at c = 0 and in the literal gauge at
c = 2: identical, and different once the corrections of first_order() are removed.  The plain operator at c = 2 with
K = (1 + s/2)^2 against gravscalar_operator_derive.el_at_origin on the canonical field with the mapped base field, six base
fields: equal -- the K-weighted c-family IS step 3's c-family in the canonical field, at the level of the operator.  A spawn
pool with the real _init/_task.  Every late code path (--mode all, six cases, translate, literal, --r8) through the real main()
with a synthetic operator: no crash.
What the checks do not carry (its mutation table):
  * R6 and R7 are consequences of R2, R3, R5 and the fitted step 3 formula; they catch a wrong argument to (3.40) or a wrong
    ghost coefficient, nothing in E.  R5 compares the ghost with the coefficient it was given.
  * At c = 1 NOTHING sees the index order of d_c B^m or the commutator term of Omega (omega is proportional to one matrix, and
    Tr Omega^2 is blind to a sign).  Both are right by reading; at c != 1 they matter and only R8a would see them.
  * Case a = 0 sees no K-dependent mutation.  R4 rests on K2 being drawn nonzero and R3 passing.
  * In --mode min (all R8 runs) Laplace type is not tested; in --mode all only at the origin.
  * At c != 1 the per-background checks are S and the ghost only: E, omega, Omega there are validated by R8a-R8d and nothing else.
Stated too strongly above, corrected here:
  * R8 is not a new gauge test "of the section 7 theory": on shell delta = 0, so the K-theory's total at c equals the
    canonical theory's at c exactly.  R8 at a = 1/12 is step 3's c-independence test plus a test of this code's handling of K.
  * R8d, added on the reviewer's argument: b = 0 for every c, by the parity (chi, psi) -> (-chi, -psi).
  * Taken on trust by the reviewer and by me: that the on-shell divergence is c-independent up to total derivatives -- that is
    the thing R8a tests, not something established -- and that (v^2)^2 is not itself a total derivative on shell.

Run (from graviton_loop; the author's machine; three runs, in this order, each resumable):
      python3 class710_operator_derive.py --jobs 16 --r8 --c 0 --cases 1/12 --onshell 5 --offshell 0 --mode min --crosscheck 2 --log results/class710_r8_c0.log
      python3 class710_operator_derive.py --jobs 16 --r8 --c 2 --cases 1/12 --onshell 5 --offshell 0 --mode min --crosscheck 2 --log results/class710_r8_c2.log
      python3 class710_operator_derive.py --jobs 16 --log results/class710_operator_derive.log
  The first two are R8 (five on-shell backgrounds each, what E, omega, Omega need); the third is R0-R7: six cases x (1 on-shell +
  1 off-shell) backgrounds, every base field, six of them cross-checked by translation on the first background.
  Forecast from the smokes (sandbox core-seconds / 3): R8 about 15 minutes per background, 1 1/4 hours per run; the main run about
  30 minutes per background, 6 hours.  Forecasts for this machine have been wrong in both directions; each background prints a
  timing line.
  RESUME.  Every finished base field is written at once to <log>.ckpt.  If a run stops for any reason, give the same command
  again with --resume added: nothing finished is recomputed, and the log is appended to.  Without --resume a non-empty log or
  checkpoint is refused, never overwritten.  Give a resumed run the same options; a line of the checkpoint is used only for
  the background, gauge, c and code version it was computed for.
  The c = 0 run doubles as the control for R1: it must report more than 4 nonzero first-derivative entries and
  "E == the c = 1 prediction: False".
"""
import sys, os, time, argparse, random, itertools, platform, json, hashlib
import sympy as sp
import flint

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gravscalar_operator_derive as D          # sympy-level helpers, and its Background for the ghost script
import ghost_operator_derive as GH
import hv1974_gravscalar_M as GM

N = 4
idx = range(N)
PAIRS = D.PAIRS
COMPS = D.COMPS
out = None
CTX = flint.fmpq_mpoly_ctx.get(('y1', 'y2', 'y3', 'y4', 't1', 't2', 't3', 't4', 'eps', 'dl'), 'degrevlex')
GENS = CTX.gens()
Y = GENS[:4]; TT = GENS[4:8]; EPS, DL = GENS[8], GENS[9]
ZERO = CTX.from_dict({})
ONE = CTX.constant(flint.fmpq(1))
HALF = sp.Rational(1, 2)


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out:
        out.write(s + '\n'); out.flush()


def q(r):
    r = sp.Rational(r)
    return flint.fmpq(int(r.p), int(r.q))


def C(r):
    return CTX.constant(q(r))


def rat(c):
    return sp.Rational(int(c.p), int(c.q))


def psum(it):
    tot = ZERO
    for e in it:
        tot = tot + e
    return tot


def d(p, m):
    return p.derivative(m)


def make_T(tdeg, o):
    def T(p, ydeg):
        """keep monomials with y-degree <= ydeg + o, t-degree <= tdeg, eps-degree <= 2, dl-degree <= 1"""
        yd = ydeg + o
        return CTX.from_dict({m: c for m, c in p.to_dict().items()
                              if m[0] + m[1] + m[2] + m[3] <= yd and m[4] + m[5] + m[6] + m[7] <= tdeg and m[8] <= 2 and m[9] <= 1})
    return T


def lin_expand(p, lo):
    """value at the origin and the four first derivatives, in t (lo = 4: the y = 0 part of p) or in y (lo = 0: the t = 0 part)"""
    v = sp.Integer(0); g = [sp.Integer(0)] * 4
    other = 4 - lo
    for m, c in p.to_dict().items():
        if sum(m[other:other + 4]) or m[8] or m[9]:
            continue
        st = sum(m[lo:lo + 4])
        if st == 0:
            v = rat(c)
        elif st == 1:
            g[list(m[lo:lo + 4]).index(1)] = rat(c)
    return v, g


class Background:
    """metric and scalar background; K(phi) = K0 + K1 s + 1/2 K2 s^2, s = phi - phi-bg(origin).  Three ways to run it:
      (tdeg, o) = (0, 0)  'plain':      everything to the orders of gravscalar_operator_derive.py; the operator at the origin only.
      (tdeg, o) = (1, 0)  'translate':  g(y + t), phi-bg(y + t) exact to first order in t; the operator at the point t (the
                                        extension proven in conjugation_proving.py).
      (tdeg, o) = (0, 1)  'expand':     no translation; every truncation one order higher in y, and the Euler-Lagrange expression
                                        kept to first order in y.  The same information as 'translate', several times cheaper."""
    def __init__(self, R, v, W, cub, Kj, tdeg, o=0):
        assert not (tdeg and o)
        self.tdeg = tdeg; self.o = o; T = self.T = make_T(tdeg, o)
        self.Kj = [sp.Rational(x) for x in Kj]

        def shift(p):
            if tdeg == 0:
                return p
            return p + psum(TT[m] * d(p, m) for m in idx)
        g = [[None] * 4 for _ in idx]
        for mu in idx:
            for nu in idx:
                q2 = psum(C(-sp.Rational(1, 3) * R[(mu, a, nu, b)]) * Y[a] * Y[b] for a in idx for b in idx)
                q4 = psum(C(sp.Rational(2, 45) * R[(a, mu, b, l)] * R[(l, c, nu, e)]) * Y[a] * Y[b] * Y[c] * Y[e]
                          for a in idx for b in idx for c in idx for e in idx for l in idx)
                g[mu][nu] = shift((ONE if mu == nu else ZERO) + q2 + q4)
        self.gb = g
        self.gbinv = self.series_inverse(g, 3)
        dgb = [[[d(g[a][b], c) for c in idx] for b in idx] for a in idx]
        self.Gb = [[[T(psum(self.gbinv[m][l] * (dgb[l][b][a] + dgb[l][a][b] - dgb[a][b][l]) for l in idx) * C(HALF), 3)
                     for b in idx] for a in idx] for m in idx]                    # Gb[m][a][b] = Gamma^m_{ab}
        A = [[g[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]; A2 = self.mat_mul(A, A, 2)
        logdet = T(psum(A[i][i] for i in idx) - psum(A2[i][i] for i in idx) * C(HALF), 2)
        self.sqrtgb = T(ONE + logdet * C(HALF) + logdet * logdet * C(sp.Rational(1, 8)), 2)
        phi = psum(C(v[a]) * Y[a] for a in idx) + psum(C(W[a, b] / 2) * Y[a] * Y[b] for a in idx for b in idx)
        phi = phi + psum(C(c) * Y[0]**m[0] * Y[1]**m[1] * Y[2]**m[2] * Y[3]**m[3] for m, c in cub.items())
        self.phit = shift(phi)                                                    # phi-bg(origin) = 0: only its derivatives and K's jets enter
        self.dphit = [d(self.phit, a) for a in idx]
        self.Kbg = self.K_of(self.phit, 3)

    def K_of(self, s, ydeg):
        K0, K1, K2 = self.Kj
        return self.T(C(K0) + C(K1) * s + C(K2 / 2) * s * s, ydeg)

    def mat_mul(self, A, B, ydeg):
        return [[self.T(psum(A[i][k] * B[k][j] for k in idx), ydeg) for j in idx] for i in idx]

    def series_inverse(self, g, ydeg):
        A = [[g[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]
        eye = [[ONE if i == j else ZERO for j in idx] for i in idx]
        inv = [row[:] for row in eye]; term = [row[:] for row in eye]
        for _ in range(6):
            term = self.mat_mul([[-t for t in row] for row in term], A, ydeg)
            inv = [[inv[i][j] + term[i][j] for j in idx] for i in idx]
        return inv

    def L2(self, h, ph, c_gauge, gauge_K):
        """the eps^2 coefficient of  sqrt(G)(R[G] - 1/2 K(Phi) G^{mn} d Phi d Phi) - 1/2 sqrt(g) g^{mn} C_m C_n,  G = g + eps h,
        Phi = phi-bg + eps ph,  C_m = g^{nr} D_r h_mn - 1/2 D_m h - c Gf ph d_m phi-bg,  Gf = K(phi-bg) (image gauge) or 1; to y^2"""
        T = self.T
        G = [[self.gb[a][b] + EPS * h[a][b] for b in idx] for a in idx]
        Ginv = self.series_inverse(G, 3)
        A = [[G[i][j] - (ONE if i == j else ZERO) for j in idx] for i in idx]; A2 = self.mat_mul(A, A, 2); A3 = self.mat_mul(A2, A, 2)
        tr = lambda M: psum(M[i][i] for i in idx)
        logdet = T(tr(A) - tr(A2) * C(HALF) + tr(A3) * C(sp.Rational(1, 3)), 2)
        l2 = T(logdet * logdet, 2)
        sqrtg = T(ONE + logdet * C(HALF) + l2 * C(sp.Rational(1, 8)) + T(l2 * logdet, 2) * C(sp.Rational(1, 48)), 2)
        dG = [[[T(d(G[a][b], c), 3) for c in idx] for b in idx] for a in idx]
        Gam = [[[T(psum(Ginv[m][l] * (dG[l][b][a] + dG[l][a][b] - dG[a][b][l]) for l in idx) * C(HALF), 3)
                 for b in idx] for a in idx] for m in idx]
        Ric = [[None] * 4 for _ in idx]
        for b in idx:
            for c in idx:
                if c < b:
                    Ric[b][c] = Ric[c][b]; continue
                Ric[b][c] = T(psum(d(Gam[a][b][c], a) - d(Gam[a][b][a], c) for a in idx)
                              + psum(Gam[a][a][e] * Gam[e][b][c] - Gam[a][c][e] * Gam[e][b][a] for a in idx for e in idx), 2)
        Rsc = T(psum(Ginv[b][c] * Ric[b][c] for b in idx for c in idx), 2)
        dPhi = [self.dphit[a] + EPS * d(ph, a) for a in idx]
        kin = T(psum(Ginv[m][n] * dPhi[m] * dPhi[n] for m in idx for n in idx), 2)
        Kfull = self.K_of(self.phit + EPS * ph, 2)
        L = T(sqrtg * (Rsc - T(Kfull * kin, 2) * C(HALF)), 2)
        Gb, gbinv = self.Gb, self.gbinv
        Dh = [[[T(d(h[a][b], c) - psum(Gb[l][c][a] * h[l][b] + Gb[l][c][b] * h[a][l] for l in idx), 3)
                for c in idx] for b in idx] for a in idx]
        Dtrh = [T(psum(gbinv[a][b] * Dh[a][b][c] for a in idx for b in idx), 3) for c in idx]
        Gf = self.Kbg if gauge_K else ONE
        Cl = [T(psum(gbinv[n][r] * Dh[m][n][r] for n in idx for r in idx) - Dtrh[m] * C(HALF) - C(c_gauge) * T(Gf * ph, 3) * self.dphit[m], 3) for m in idx]
        C2 = T(self.sqrtgb * psum(gbinv[m][n] * Cl[m] * Cl[n] for m in idx for n in idx), 2)
        L = T(L - C(HALF) * EPS * EPS * C2, 2)
        return CTX.from_dict({m[:8] + (0,) + m[9:]: c for m, c in L.to_dict().items() if m[8] == 2})

    def gam1(self):
        """d_c (g^{ab} Gamma^m_{ab}) at the origin: gam1[c][m]  (g^{ab} = delta^{ab} + second order, Gamma first order)"""
        assert self.tdeg == 1 or self.o == 1
        e = [lin_expand(self.T(psum(self.Gb[m][a][a] for a in idx), 0), 4 if self.tdeg else 0)[1] for m in idx]
        return [[e[m][c] for m in idx] for c in idx]


def comp_field(comp, mono):
    h = [[ZERO for _ in idx] for _ in idx]; ph = ZERO
    if comp == 's':
        ph = mono
    else:
        a, b = comp
        h[a][b] = mono
        if a != b:
            h[b][a] = mono
    return h, ph


def probe_fields():
    P = []
    for comp in COMPS:
        P.append(('0', comp, None, comp_field(comp, ONE)))
        for m in idx:
            P.append(('1', comp, m, comp_field(comp, Y[m])))
        for m in idx:
            for n in idx:
                if m <= n:
                    w = C(HALF) if m == n else ONE
                    P.append(('2', comp, (m, n), comp_field(comp, Y[m] * Y[n] * w)))
    return P


def coeff_dl(p):
    return CTX.from_dict({m[:9] + (0,): c for m, c in p.to_dict().items() if m[9] == 1})


def el(bg, base, c_gauge, gauge_K, probes):
    """the 11 Euler-Lagrange components on the base field, each as (value at the origin, four first derivatives): in t at
    y = 0 ('translate': the operator at the point t on a base field that is a monomial in y) or in y ('expand': the
    Euler-Lagrange expression of the base field as a field, to first order; the caller removes the base field's own y-dependence)"""
    bh, bp = base
    T = bg.T
    G = {}
    for kind, comp, key, (ph_h, ph_p) in probes:
        h = [[bh[a][b] + DL * ph_h[a][b] for b in idx] for a in idx]
        G[(kind, comp, key)] = coeff_dl(bg.L2(h, bp + DL * ph_p, c_gauge, gauge_K))
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
        if comp != 's' and comp[0] != comp[1]:
            v0 = v0 / 2; g1 = [x / 2 for x in g1]          # the independent variable h_ab (a < b) enters as h_ab + h_ba
        outv.append((v0, g1))
    return outv


def base_from_spec(spec):
    h = [[ZERO for _ in idx] for _ in idx]; ph = ZERO
    for k, mon, c in spec:
        term = C(c) * Y[0]**mon[0] * Y[1]**mon[1] * Y[2]**mon[2] * Y[3]**mon[3]
        if COMPS[k] == 's':
            ph = ph + term
        else:
            a, b = COMPS[k]
            h[a][b] = h[a][b] + term
            if a != b:
                h[b][a] = h[b][a] + term
    return h, ph


_W = {}
WAYS = {'plain': (0, 0), 'translate': (1, 0), 'expand': (0, 1)}


def _init(Rd, v, W, cub, Kj, c_gauge, gauge_K):
    _W['args'] = (Rd, v, W, cub, Kj); _W['bg'] = {}; _W['probes'] = probe_fields(); _W['c'] = c_gauge; _W['gK'] = gauge_K


def _task(job):
    tag, way, spec = job
    t0 = time.time()
    if way not in _W['bg']:
        _W['bg'][way] = Background(*_W['args'], *WAYS[way])
    r = el(_W['bg'][way], base_from_spec(spec), _W['c'], _W['gK'], _W['probes'])
    return tag, way, r, time.time() - t0


def mono(*ms):
    e = [0, 0, 0, 0]
    for m in ms:
        e[m] += 1
    return tuple(e)


CROSS = [('lin', 10, 0), ('rad', 10), ('lin', 1, 2), ('rad', 0), ('lin', 7, 1), ('lin', 4, 3)]


def task_list(mode, cols, method, ncross):
    """mode 'all': every base field (165 + the 11 radial ones); 'min': what E, omega, Omega need (66); 'cols': a smoke -- the
    linear base fields to first order only for the columns `cols` of E.  `method` is how the first-order data are obtained;
    the first `ncross` base fields of CROSS are run by the other method as well."""
    J = []
    spec_of = {}
    for k in range(11):
        spec_of[('const', k)] = [(k, (0, 0, 0, 0), 1)]
        J.append((('const', k), 'plain'))
        for m in idx:
            spec_of[('lin', k, m)] = [(k, mono(m), 1)]
            J.append((('lin', k, m), method if (mode != 'cols' or k in cols) else 'plain'))
        spec_of[('rad', k)] = [(k, mono(m, m), HALF) for m in idx]
        J.append((('rad', k), method))
    if mode == 'all':
        for k in range(11):
            for m in idx:
                for n in idx:
                    if m <= n:
                        spec_of[('quad', k, m, n)] = [(k, mono(m, n), HALF if m == n else 1)]
                        J.append((('quad', k, m, n), 'plain'))
    other = 'translate' if method == 'expand' else 'expand'
    for tag in CROSS[:ncross]:
        if mode != 'cols' or tag[1] in cols:
            J.append((tag, other))
    return [(tag, way, spec_of[tag]) for tag, way in J]


def first_order(tag, way, r, C0, B0):
    """value and first derivatives of the coefficient field a base field measures (K^mm summed for 'rad', B^m for 'lin'), as
    11-vectors: for 'translate' as returned; for 'expand' after removing the base field's own y-dependence:
        O[y_m e_k](y)       = B^m(y) e_k + y_m C(y) e_k                 ->  d_c B^m = gradient_c - delta_cm C(0)
        O[1/2 y.y e_k](y)   = sum_m K^mm(y) e_k + y_m B^m(y) e_k + ...  ->  d_c sum K^mm = gradient_c - B^c(0)"""
    v0 = sp.Matrix([x[0] for x in r]); g = [sp.Matrix([x[1][c] for x in r]) for c in idx]
    if way == 'expand':
        k = tag[1]
        if tag[0] == 'lin':
            g[tag[2]] = g[tag[2]] - C0[:, k]
        elif tag[0] == 'rad':
            g = [g[c] - B0[c][:, k] for c in idx]
    return v0, g


def _code_version():
    """a hash of this file and the three it imports: a checkpoint line written by different code is never reused"""
    h = hashlib.sha1()
    here = os.path.dirname(os.path.abspath(__file__))
    for fn in ('class710_operator_derive.py', 'gravscalar_operator_derive.py', 'ghost_operator_derive.py', 'hv1974_gravscalar_M.py'):
        with open(os.path.join(here, fn), 'rb') as fh:
            h.update(fh.read().replace(b'\r\n', b'\n'))
    return h.hexdigest()[:12]


CODE_VERSION = _code_version()


class Checkpoint:
    """every finished base field is appended to <log>.ckpt as one JSON line (exact rationals as strings) and flushed; a run
    restarted with --resume reads the file and computes only what is missing.  The key names the background by a hash of all
    its data, the gauge and c, so a line can never be used for a different computation."""
    def __init__(self, path, resume):
        self.path = path; self.done = {}; self.f = None
        if path is None:
            return
        if os.path.exists(path) and os.path.getsize(path) > 0:
            if not resume:
                print(f"refusing to overwrite the non-empty checkpoint {path}; add --resume to continue that run, or to start over:  del {path.replace('/', chr(92))}"); sys.exit(1)
            with open(path) as fh:
                for line in fh:
                    try:
                        rec = json.loads(line)
                    except ValueError:
                        continue                      # a line cut off by a crash: that base field is simply computed again
                    self.done[rec['k']] = ([(sp.Rational(a), [sp.Rational(x) for x in g]) for a, g in rec['r']], rec['dt'])
        ends_clean = True
        if os.path.exists(path) and os.path.getsize(path) > 0:
            with open(path, 'rb') as fh:
                fh.seek(-1, 2); ends_clean = fh.read(1) == b'\n'
        self.f = open(path, 'a')
        if not ends_clean:
            self.f.write('\n')                       # a line cut off by a crash: the next record must not be glued to it

    @staticmethod
    def bgkey(Rd, v, W, cub, Kj, c_gauge, gauge_K):
        txt = repr((CODE_VERSION, sorted((k, str(x)) for k, x in Rd.items()), [str(x) for x in v], [str(x) for x in W], sorted((k, str(x)) for k, x in cub.items()),
                    [str(x) for x in Kj], str(c_gauge), bool(gauge_K)))
        return hashlib.sha1(txt.encode()).hexdigest()[:16]

    def put(self, key, r, dt):
        if self.f:
            self.f.write(json.dumps({'k': key, 'r': [[str(a), [str(x) for x in g]] for a, g in r], 'dt': round(dt, 1)}) + '\n'); self.f.flush(); os.fsync(self.f.fileno())


def derive(Rd, v, W, cub, Kj, c_gauge, gauge_K, mode, cols, method, ncross, jobs, label, t00, ck):
    J = task_list(mode, cols, method, ncross)
    bk = Checkpoint.bgkey(Rd, v, W, cub, Kj, c_gauge, gauge_K)
    keyof = lambda tag, way: f"{bk}|{tag}|{way}"
    way_of = {}
    for tag, way, _ in J:                                  # the first entry for a base field is the one used; a later one is the cross-check
        way_of.setdefault(tag, way)
    cost = {'translate': 0, 'expand': 1, 'plain': 2}
    J.sort(key=lambda j: cost[j[1]])                       # the slow tasks first
    res = {}; times = {}; ntot = len(J)
    todo = []
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
        with mp.Pool(jobs, initializer=_init, initargs=(Rd, v, W, cub, Kj, c_gauge, gauge_K)) as pool:
            for tag, way, r, dt in pool.imap_unordered(_task, todo, chunksize=1):
                take(tag, way, r, dt)
    elif todo:
        _init(Rd, v, W, cub, Kj, c_gauge, gauge_K)
        for job in todo:
            take(*_task(job))
    tline = "; ".join(f"{k}: {len(v_)} tasks, mean {sum(v_) / len(v_):.1f} s, sum {sum(v_):.0f} s" for k, v_ in times.items())
    val0 = lambda tag: sp.Matrix([x[0] for x in res[(tag, way_of[tag])]])
    C0 = sp.Matrix.hstack(*[val0(('const', k)) for k in range(11)])
    B0 = [sp.Matrix.hstack(*[val0(('lin', k, m)) for k in range(11)]) for m in idx]
    fo = lambda tag, way: first_order(tag, way, res[(tag, way)], C0, B0)
    z = [sp.zeros(11, 1)] * 4
    grad = lambda tag: fo(tag, way_of[tag])[1] if way_of[tag] != 'plain' else z
    B1 = [[sp.Matrix.hstack(*[grad(('lin', k, m))[c] for k in range(11)]) for c in idx] for m in idx]       # B1[m][c] = d_c B^m
    S0 = sp.Matrix.hstack(*[val0(('rad', k)) for k in range(11)]) / 4
    S1 = [sp.Matrix.hstack(*[grad(('rad', k))[c] for k in range(11)]) / 4 for c in idx]
    laplace = None
    if mode == 'all':
        laplace = all(val0(('quad', k, m, n)) == (S0[:, k] if m == n else sp.zeros(11, 1)) for k in range(11) for m in idx for n in idx if m <= n)
    other = 'translate' if method == 'expand' else 'expand'
    cross = {}
    for tag in CROSS[:ncross]:
        if way_of[tag] == method and (tag, other) in res:
            cross[tag] = fo(tag, method) == fo(tag, other)
    return dict(C0=C0, B0=B0, B1=B1, S0=S0, S1=S1, laplace=laplace, tline=tline, cross=cross)


def operator_data(dd, gam1, mutate):
    """omega, E, Omega at the origin from the measured coefficient fields (the formulas of conjugation_proving.py's Direct.data)"""
    S0, S1, B0, B1, C0 = dd['S0'], dd['S1'], dd['B0'], dd['B1'], dd['C0']
    Si = S0.inv()
    if mutate == 'nosinv':
        S0m = S0.copy(); S0m[10, 10] = 1; Si = S0m.inv()
    f = HALF if mutate != 'halfw' else sp.Rational(1, 4)
    w0 = [Si * B0[m] * f for m in idx]
    w1 = [[(Si * B1[m][c] - Si * S1[c] * Si * B0[m]) * f + gam1[c][m] * sp.eye(11) / 2 for m in idx] for c in idx]       # w1[c][m] = d_c w^m
    c0 = Si * C0
    E = c0 - sum((w1[m][m] for m in idx), sp.zeros(11, 11)) - sum((w0[m] * w0[m] for m in idx), sp.zeros(11, 11))
    Om = {(m, n): w1[m][n] - w1[n][m] + w0[m] * w0[n] - w0[n] * w0[m] for m in idx for n in idx}
    Om2 = sp.expand(sum((Om[k] * Om[k]).trace() for k in Om))
    return w0, E, Om2


def K_jets_of_a(a, phi0):
    """K_a = 1/f + 3 f'^2/f^2, f = 1 - a phi^2  (weyl_map_check.py's K at c = 1): its value and two derivatives at phi0"""
    p = sp.Symbol('p'); f = 1 - a * p**2
    K = 1 / f + 3 * sp.diff(f, p)**2 / f**2
    return [sp.Rational(sp.cancel(sp.diff(K, p, n).subs(p, phi0))) for n in range(3)]


def step3_structure(seed, nsamp=12):
    """the pre-reduction coefficients of step 3's total on the eight structures, fitted here from hv1974_gravscalar_M's own
    functions (nothing typed), and the on-shell number after its reduce()"""
    rng = random.Random(seed); rows = []
    for _ in range(nsamp):
        R = GM.riemann_sample(rng); v, W = GM.scalar_sample(rng)
        Ric, Rs, inv = GM.invariants(R, v, W)
        E, Om2 = GM.build_E(R, Ric, Rs, v, W)
        tot = sp.expand(GM.hv340(11, Om2, E, Rs, inv['Ric2']) - 2 * GM.hv340(4, -inv['Riem2'], Ric - v * v.T, Rs, inv['Ric2']))
        rows.append((tot, inv))
    c, ok = GM.fit(rows)
    red = GM.reduce(c)
    on = sp.expand(red['R2'] / 4 + red['Ric2'] / 4 + red['Rv2'] / 2 + red['Ricvv'] / 2 + red['v4'])
    return c, red, on, ok


R8KEYS = ['v4', 'vWv', 'ricid', 'gb']


def fit_r8(rows):
    """on-shell totals -> (alpha, b, C, g) on the four structures from the first four rows; residuals of the others"""
    cs = sp.symbols('r0:4')
    eqs = [sum(c * st[k] for c, k in zip(cs, R8KEYS)) - val for val, st in rows]
    sol = sp.solve(eqs[:4], cs, dict=True)
    if not sol or len(sol[0]) < 4:
        return None, None
    sol = sol[0]
    return {k: sol[c] for c, k in zip(cs, R8KEYS)}, [sp.expand(e.subs(sol)) for e in eqs[4:]]


def r8_selftest(seed):
    """the fit on the known case: step 3's formula (c = 1) on six random on-shell backgrounds of the canonical theory"""
    rng = random.Random(seed); rows = []
    for _ in range(6):
        R = GM.riemann_sample(rng); v, W = GM.scalar_sample(rng)
        W = W - W.trace() * sp.eye(4) / 4
        R = D.with_ricci(R, v * v.T / 2)
        Ric, Rs, inv = GM.invariants(R, v, W)
        E, Om2 = GM.build_E(R, Ric, Rs, v, W)
        tot = sp.expand(GM.hv340(11, Om2, E, Rs, inv['Ric2']) - 2 * GM.hv340(4, -inv['Riem2'], Ric - v * v.T, Rs, inv['Ric2']))
        rows.append((tot, {'v4': inv['v4'], 'vWv': (v.T * W * v)[0], 'ricid': inv['WW'] - inv['trW2'] + inv['Ricvv'], 'gb': inv['Riem2'] - 4 * inv['Ric2'] + inv['R2']}))
    return fit_r8(rows)


def run_background(case, kind, R, u, W, cub, Kj, args, struct, ncross, t00, ck):
    K0, K1, K2 = Kj
    k1 = K1 / K0
    u2 = (u.T * u)[0]
    Wt = W + k1 * u * u.T / 2
    beta = Wt.trace()
    delta = sp.expand(k1 * beta / 2)                # the header's first wording had an extra 1/4 k1^2 u^2: see CORRECTION in the docstring
    Ric, Rs, inv = GM.invariants(R, u, Wt)
    onsh = (Ric - K0 * u * u.T / 2) == sp.zeros(4, 4) and beta == 0
    label = f"{case} {kind}"
    say(f"\n## {label}: K jets ({K0}, {K1}, {K2}), k1 = {k1}; R = {Rs}, u^2 = {u2}, beta = {beta}, delta = {delta}; field equations hold at 0: {onsh}")
    gauge_K = args.gauge == 'image'
    cols = [int(x) for x in args.cols.split(',')] if args.mode == 'cols' else list(range(11))
    dd = derive(R, u, W, cub, Kj, args.c, gauge_K, args.mode, cols, args.method, ncross, args.jobs, label, t00, ck)
    say(f"   timing: {dd['tline']}")
    if dd['cross']:
        say(f"   cross-check, first-order data by 'translate' and by 'expand' on {len(dd['cross'])} base fields: {'identical' if all(dd['cross'].values()) else 'DIFFER: ' + str(dd['cross'])}")
        flags_cross = all(dd['cross'].values())
    else:
        flags_cross = None
    gam1 = Background(R, u, W, cub, Kj, 0, 1).gam1()
    w0, E, Om2 = operator_data(dd, gam1, args.mutate)
    Pbig = sp.diag(D.op_matrix10(lambda h: h / 2 - sp.eye(4) * h.trace() / 4), sp.Matrix([[1]]))
    ess = sp.zeros(11, 11); ess[10, 10] = 1
    Sw = Pbig.copy(); Sw[10, 10] = K0
    flags = {}
    if flags_cross is not None:
        flags['cross'] = flags_cross
    pred_ok = gauge_K and args.c == 1                 # R1-R3, R6, R7 are statements about the image gauge at c = 1
    # R1
    r1s = dd['S0'] == Sw and all(dd['S1'][c] == K1 * u[c] * ess for c in idx)
    nfirst = sum(1 for m in idx for x in dd['B0'][m] if x != 0)
    r1w = all(w0[m] == k1 * u[m] / 2 * ess for m in idx)
    say(f"   R1  S = diag(P, K0) with d_c S = K1 u_c on the scalar only: {r1s}" + (f"; Laplace type on all 110 quadratic base fields: {dd['laplace']}" if dd['laplace'] is not None else "")
        + f";  nonzero first-derivative entries of the raw operator: {nfirst};  omega_m = 1/2 k1 u_m on the scalar only: {r1w}")
    flags['R1'] = (r1s and r1w and dd['laplace'] is not False) if pred_ok else (r1s and dd['laplace'] is not False)       # S does not depend on the gauge
    # predicted E: step 3's printed-and-derived E (hv1974_gravscalar_M.build_E), by homogeneity in v and W, plus delta
    E_u, Om2_step3 = GM.build_E(R, Ric, Rs, u, Wt)
    E_0, _ = GM.build_E(R, Ric, Rs, sp.zeros(4, 1), sp.zeros(4, 4))
    Ep = sp.zeros(11, 11)
    Ep[:10, :10] = E_0[:10, :10] + K0 * (E_u[:10, :10] - E_0[:10, :10])
    Ep[:10, 10] = K0 * E_u[:10, 10]
    Ep[10, :10] = E_u[10, :10]
    Ep[10, 10] = K0 * E_u[10, 10] + delta
    Ed = E.applyfunc(sp.expand); Ep = Ep.applyfunc(sp.expand)
    colsok = {k: Ed[:, k] == Ep[:, k] for k in cols}
    r3 = all(colsok.values())
    if args.mutate == 'none' and pred_ok:             # the two controls that only touch the conversion to E cost nothing: done in every run
        mus = ('nosinv', 'halfw') if K0 != 1 else ('halfw',)          # at K0 = 1 leaving the scalar's S^-1 out changes nothing (review, 13 Nov)
        ctl = {mu: any(operator_data(dd, gam1, mu)[1].applyfunc(sp.expand)[:, k] != Ep[:, k] for k in cols) for mu in mus}
        say("   controls (must break R3): scalar S^-1 left out: " + ("not applicable at K0 = 1" if K0 == 1 else 'breaks it' if ctl['nosinv'] else 'DOES NOT') + f";  omega halved: {'breaks it' if ctl['halfw'] else 'DOES NOT'}")
        flags['controls'] = all(ctl.values())
    if args.mode == 'cols':
        say(f"   R3  E == prediction, columns {cols} only (smoke): {colsok}" + ("" if pred_ok else "   [not a prediction at this c / gauge]"))
        say(f"       scalar-scalar entry: derived {Ed[10, 10]}, predicted -K0 u^2 + delta = {Ep[10, 10]}" if 10 in cols else "")
        if pred_ok:
            flags['R3'] = r3
        return flags, None
    r2 = Om2 == Om2_step3
    if not pred_ok:
        say(f"   (c = {args.c}, gauge {args.gauge}: R2, R3, R6, R7 are not predictions here; for the record: Tr Omega^2 == step 3's: {r2}; E == the c = 1 prediction: {r3})")
    else:
        say(f"   R2  Tr Omega^2 = step 3's ({Om2_step3}): {r2}")
    if pred_ok:
        say(f"   R3  E == prediction (11 x 11): {r3}   [h-h block {Ed[:10, :10] == Ep[:10, :10]}; (E phi)_ab = -2 K0 Wt: {Ed[:10, 10] == Ep[:10, 10]}; (E h)_phi: {Ed[10, :10] == Ep[10, :10]}; E_phiphi = -K0 u^2 + delta: {Ed[10, 10] == Ep[10, 10]} (derived {Ed[10, 10]})]")
        flags['R2'] = r2; flags['R3'] = r3
    # R5: the ghost, by ghost_operator_derive.py on the same background; the gauge condition's scalar coefficient at the origin is c K0 (image) or c (literal)
    c_eff = args.c * K0 if gauge_K else args.c
    s_gh, second_ok, first_zero, E_gh, sset = GH.ghost_operator(D.Background(R, u, W, cub), c_eff, 1)
    r5 = s_gh == 1 and second_ok and first_zero and E_gh is not None and E_gh.applyfunc(sp.expand) == (Ric - c_eff * u * u.T)
    say(f"   R5  ghost: unit second-derivative part {s_gh == 1 and second_ok}, no first-derivative terms {first_zero}, E_gh == Ric - ({c_eff}) u u: {r5}")
    flags['R5'] = r5
    # R6, R7: (3.40)
    hphi = GM.hv340(11, Om2, Ed, Rs, inv['Ric2'])
    ghost = -2 * GM.hv340(4, -inv['Riem2'], E_gh, Rs, inv['Ric2'])
    tot = sp.expand(hphi + ghost)
    c_, red, on, _ = struct
    img = dict(inv); img['Rv2'] = K0 * inv['Rv2']; img['Ricvv'] = K0 * inv['Ricvv']; img['v4'] = K0**2 * inv['v4']; img['trW2'] = K0 * inv['trW2']; img['WW'] = K0 * inv['WW']
    pred = sp.expand(sum(c_[k] * img[k] for k in GM.KEYS) + delta * (Rs / 6 - K0 * u2) / 2 + delta**2 / 4)
    r6 = tot == pred
    row = None
    if onsh:                                          # the four on-shell structures of R8, in rational form (b is reported per K0 u.Wt.u)
        row = (tot, {'v4': img['v4'], 'vWv': K0 * (u.T * Wt * u)[0], 'ricid': img['WW'] - img['trW2'] + img['Ricvv'], 'gb': inv['Riem2'] - 4 * inv['Ric2'] + inv['R2']})
    if not pred_ok:
        say(f"   (3.40) total, pointwise: {tot}")
        return flags, row
    say(f"   R6  (3.40) total, pointwise: {tot};  step 3's eight-structure formula at the image background + 1/2 delta (R/6 - K0 u^2) + 1/4 delta^2: {pred};  equal: {r6}")
    flags['R6'] = r6
    if onsh:
        gb = inv['Riem2'] - 4 * inv['Ric2'] + inv['R2']
        ricid = img['WW'] - img['trW2'] + img['Ricvv']
        rem = sp.expand(tot - c_['Riem2'] * gb - c_['WW'] * ricid)
        r7 = rem == on * K0**2 * u2**2
        say(f"   R7  on shell: total - ({c_['Riem2']}) Gauss-Bonnet - ({c_['WW']}) (W.W - (tr W)^2 + Ric v v) = {rem};  ({on}) K0^2 (u^2)^2 = {sp.expand(on * K0**2 * u2**2)};  equal: {r7}")
        flags['R7'] = r7
    return flags, row


def main():
    global out
    ap = argparse.ArgumentParser()
    ap.add_argument('--seed', type=int, default=20261113)
    ap.add_argument('--cases', default='1/12,free,0,1/6,1,rand', help="values of a (K = K_a at phi0), 'rand' (a random rational a), 'free' (random K jets)")
    ap.add_argument('--phi0', default='1/2')
    ap.add_argument('--onshell', type=int, default=1)
    ap.add_argument('--offshell', type=int, default=1)
    ap.add_argument('--c', default='1', help='coefficient of the scalar term of the gauge condition (a rational; 1 is the image gauge)')
    ap.add_argument('--r8', action='store_true', help='the gauge-independence test: fit the on-shell totals of this run (needs --onshell 5 or more on one case)')
    ap.add_argument('--resume', action='store_true', help='continue a run from its checkpoint file (<log>.ckpt); the log is appended to')
    ap.add_argument('--gauge', choices=('image', 'literal'), default='image')
    ap.add_argument('--mode', choices=('all', 'min', 'cols'), default='all')
    ap.add_argument('--cols', default='10', help='with --mode cols: the columns of E to test (0-9 the h pairs, 10 the scalar)')
    ap.add_argument('--method', choices=('expand', 'translate'), default='expand', help='how the first-order data of the coefficient fields are obtained')
    ap.add_argument('--crosscheck', type=int, default=6, help='on the first background, this many base fields are run by the other method too (0-6)')
    ap.add_argument('--mutate', choices=('none', 'nosinv', 'halfw'), default='none')
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--log', default=None)
    args = ap.parse_args()
    args.c = sp.Rational(args.c)
    t0 = time.time()
    if args.resume and not args.log:
        print("--resume needs --log (the checkpoint is <log>.ckpt)"); sys.exit(1)
    if args.log:
        if os.path.exists(args.log) and os.path.getsize(args.log) > 0 and not args.resume:
            print(f"refusing to overwrite the non-empty log {args.log}; add --resume to continue that run, or to start over:  del {args.log.replace('/', chr(92))}  and  del {args.log.replace('/', chr(92))}.ckpt"); sys.exit(1)
    if args.log and os.path.dirname(args.log):
        os.makedirs(os.path.dirname(args.log), exist_ok=True)
    ck = Checkpoint(args.log + '.ckpt' if args.log else None, args.resume)
    if args.log:
        out = open(args.log, 'a' if args.resume else 'w')
        if args.resume:
            say(f"\n# RESUMED {time.strftime('%Y-%m-%d %H:%M:%S')}: {len(ck.done)} base fields in the checkpoint; everything below is recomputed from it or run afresh")
    say(f"# class710_operator_derive: seed {args.seed}, cases {args.cases}, phi0 = {args.phi0}, {args.onshell} on-shell + {args.offshell} off-shell backgrounds per case, c = {args.c}, gauge {args.gauge}, mode {args.mode}, method {args.method}, crosscheck {args.crosscheck}"
        + (f" (columns {args.cols})" if args.mode == 'cols' else "") + f", mutate {args.mutate}, jobs {args.jobs}")
    say(f"# code version {CODE_VERSION} (hash of this script and its three imports; part of every checkpoint key)")
    say(f"# versions: python {platform.python_version()} ({platform.system()} {platform.machine()}), sympy {sp.__version__}, python-flint {flint.__version__}; {sys.executable}")
    struct = step3_structure(args.seed + 7)
    say(f"# step 3's total on the eight structures (fitted here from hv1974_gravscalar_M, residuals zero: {struct[3]}): {struct[0]}")
    say(f"#   after its two reductions: { {k: v for k, v in struct[1].items() if v != 0} };  on shell: {struct[2]} (v^2)^2   [{time.time() - t0:.0f} s]")
    rng = random.Random(args.seed)
    phi0 = sp.Rational(args.phi0)
    r8rows = []
    if args.r8:
        cst, rst = r8_selftest(args.seed + 11)
        st_ok = cst is not None and cst['v4'] == struct[2] and cst['vWv'] == 0 and cst['ricid'] == struct[0]['WW'] and cst['gb'] == struct[0]['Riem2'] and all(x == 0 for x in rst)
        say(f"# R8 fit on the known case (step 3's formula, c = 1, six on-shell backgrounds): {cst}, residuals {rst}: {'as it must be' if st_ok else 'WRONG -- the fit routine is not to be trusted'}")
        if not st_ok:
            say("status=FAIL"); say(f"END class710_operator_derive {time.time() - t0:.0f} s"); sys.exit(1)
    allflags = []; first = True; predictions_apply = args.gauge == 'image' and args.c == 1 and args.mutate == 'none'
    for case in args.cases.split(','):
        if case == 'free':
            Kj = [sp.Rational(rng.randint(1, 5), rng.randint(1, 3)), sp.Rational(rng.choice((-3, -2, -1, 1, 2, 3)), rng.randint(1, 3)), sp.Rational(rng.randint(-3, 3), rng.randint(1, 3))]
            name = 'free K jets'
        else:
            a = sp.Rational(rng.choice((-2, -1, 1, 2, 3)), rng.randint(2, 7)) if case == 'rand' else sp.Rational(case)
            Kj = K_jets_of_a(a, phi0); name = f"a = {a}"
        for kind, n in (('on shell', args.onshell), ('off shell', args.offshell)):
            for i in range(n):
                R = D.riemann_sample(rng); u, W, cub = D.scalar_sample(rng)
                if kind == 'on shell':
                    k1 = Kj[1] / Kj[0]
                    W = W - (W.trace() + k1 * (u.T * u)[0] / 2) * sp.eye(4) / 4
                    R = D.with_ricci(R, Kj[0] * u * u.T / 2)
                flags, row = run_background(name, f"{kind} {i}", R, u, W, cub, Kj, args, struct, args.crosscheck if first else 0, t0, ck); first = False
                allflags.append((name, kind, flags))
                if row is not None:
                    r8rows.append(row)
                say(f"   -> {name}, {kind} {i}: " + ", ".join(f"{k} {'ok' if v else 'FAILS'}" for k, v in flags.items()) + f"   [{time.time() - t0:.0f} s]")
    ok = all(all(f.values()) for _, _, f in allflags)
    complete = args.mode != 'cols' and any(k == 'on shell' for _, k, _ in allflags) and any(k == 'off shell' for _, k, _ in allflags)
    say("")
    if args.r8:
        complete = args.mode != 'cols' and len(r8rows) >= 5 and len(args.cases.split(',')) == 1
        if len(r8rows) >= 4 and len(args.cases.split(',')) == 1:
            cf, resid = fit_r8(r8rows)
            if cf is None:
                say("R8: the four structures are degenerate on these backgrounds; no fit"); ok = False
            else:
                r8a = cf['v4'] == struct[2]; r8c = all(x == 0 for x in resid)
                say(f"R8a  c = {args.c}: alpha = {cf['v4']}  (203/320 predicted; step 3's own on-shell number computed above: {struct[2]}): {'PASS' if r8a else 'FAIL'}")
                say(f"R8b  Gauss-Bonnet coefficient g = {cf['gb']} (7/12 at c = 1);  b = {cf['vWv']} per K0 u.Wt.u,  C = {cf['ricid']}")
                say(f"R8c  residuals of the {len(resid)} background(s) beyond the four fitted: {resid}: {'PASS' if (r8c and resid) else 'FAIL' if resid else 'none to test'}")
                r8d = cf['vWv'] == 0
                say(f"R8d  b = 0 (the canonical theory and this gauge family are even under (chi, psi) -> (-chi, -psi); v.Wc.v is odd): {'PASS' if r8d else 'FAIL'}")
                ok &= r8a and r8c and r8d
        else:
            say("R8: needs one case and at least four on-shell backgrounds (five for the residual test)")
    elif not predictions_apply:
        complete = False
        say(f"predictions R1-R7 are written for the image gauge with c = 1 and no mutation; this run is not that: {'every check held' if ok else 'checks fail'}")
    say("status=" + ("FAIL" if not ok else "PASS" if complete else "INCOMPLETE (PASS needs --mode all or min, and an on-shell and an off-shell background; with --r8, one case and five on-shell backgrounds; a run at c != 1, in the literal gauge or with --mutate tests no prediction without --r8)"))
    say(f"END class710_operator_derive {time.time() - t0:.0f} s")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
