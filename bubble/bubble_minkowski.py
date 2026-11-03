#!/usr/bin/env python3
"""
Rung 2, (v) the control: the complex-ghost-pair bubble with the loop energy k0 integrated ALONG THE REAL AXIS -- the naive 1969
computation, the one written before Lee-Wick's contour (Nucl. Phys. B 9, 209) and the one that 2017 (Aglietti-Anselmi 1612.06510)
names as the "Minkowski formulation".  Rule 28 (brief 16 amendment, 2 Nov 23:46): nothing from 1612.06510 or later was fetched or
read for this header; 1612.06510's ABSTRACT was read on 31 Oct (chronology, 3b) and is quoted below as memory, marked.  Everything
marked [derived] was obtained from the pole positions and the large-k expansion (sympy, 2 Nov 23:50-23:58, before this file was
written) -- 1969 mathematics.

Model and conventions as bubble_frame.py: a1 = i, a2 = -i (M = 1), metric (+,-,-,-),
    B(p) = int d^4k/(2 pi)^4 1/((k^2 - a1)((k-p)^2 - a2)),   D = 4,
rest frame p = (p0, 0), boosted frame p = (p0, p_vec) with p0^2 = s + p_vec^2, s = p^2.  The spatial integral is cut off on the
sphere |k_vec| < Lambda IN THE FRAME OF THE COMPUTATION (the cutoff is part of the formulation under test).

The k0 integral.  For real p0 the four poles k0 = +-w1, p0 +- w2 (w_i = sqrt(k_vec^2 + a_i), Re > 0) are all off the real axis for
every real k_vec: Im w1 > 0, Im w2 < 0.  The real-axis integral closes in the upper half plane and picks the two poles with Im > 0;
no prescription is needed, nothing is ambiguous, and no pole ever crosses the axis as p0 varies (real p0).  The script selects the
poles by the sign of Im, generically, and checks the residue sum against a direct quadrature of the real-axis integral (rule 6).

[derived] D1.  For (a1, a2) = (i, -i) the UHP poles are +w1 and p0 - w2, and the residue sum is
        h_{+-} = -(i/2) (w1 - w2) / (w1 w2 [p0^2 - (w1 - w2)^2])
   -- with NO positive-energy term (w1 + w2)/(...) at all: it is the mixed term of one +i eps and one -i eps propagator
   (2001.01942 eq. 2.3, bubble_lib.h_fw_mixed, which is written for the (-i, +i) ordering: the sign of (w1 - w2) flips with the
   ordering) with complex masses.  The real-axis contour treats a = +i as an anti-Feynman
   (-i eps) propagator and a = -i as a Feynman one.  (-+) is the same with w1 <-> w2.  For (-i, -i): h_{--} = the standard
   (Lee-Wick contour) value of bubble_lib.h_bubble; for (i, i): h_{++} = MINUS that value (both propagators anti-Feynman: the
   complex conjugate structure).
[derived] D2.  Large-k tail in the rest frame (u = 1/k):  h_{+-} = u^3 (p0^2 - u^2)/(2 p0^4) + O(u^7), identical for (-+).  Hence
        B_Lambda(p0) = ln(Lambda) / (4 pi^2 p0^2) + C(p0) + O(Lambda^-2),
   a LOGARITHMIC divergence with a NONLOCAL coefficient 1/p0^2 (not a polynomial in p: no counterterm absorbs it) which is REAL in B
   (the ordinary local log of a bubble is imaginary in B, +i/(8 pi^2) ln Lambda: the two sit in different parts; in 1612.06510's
   language the ordinary one is anti-Hermitian and this one is not).  The Euclidean subtraction of bubble_frame does not remove it: B_Lambda(s) - B_Lambda(s_ref) ~ (1/4 pi^2)(1/p0^2 -
   1/p0_ref^2) ln Lambda.  The next term, -u^5/(2 p0^4), converges; the u^5 term odd in (a1 - a2) vanishes for a = +-i since
   a1^2 - a2^2 = 0.
[derived] D3.  In the boosted frame the same tail holds with p0^2 = s + p_vec^2, so at fixed s the coefficient of ln Lambda is
   1/(4 pi^2 (s + p_vec^2)): the divergence itself is frame dependent.  (Nakanishi's 1971 crack, bubble_frame.py, was FINITE,
   linear in |p|, and in the imaginary part only; this one is in the divergent part and in the real part.)
[derived] D4.  The symmetric four-pair sum (++) + (--) + (+-) + (-+): the two mixed terms have the SAME tail (not opposite), so
   their divergences add; (++) + (--) = std(-i,-i) - std(i,i), finite and covariant.  No cancellation.
[memory, tainted] M1.  1612.06510's abstract (read 31 Oct, chronology 3b): the Minkowski formulation of higher-derivative
   theories has "nonlocal, non-Hermitian divergences" at one loop and is inconsistent.  D2 is the concrete form this run
   expects; whether 1612.06510's formulation and cutoff are these is to be read AFTER the run (rule 28) and logged beside this.
[obvious to a 1969 reader] O1.  Since the integrand is covariant but the integral is only conditionally convergent (the tail is
   1/k^3 against k^2 dk), the order of integration and the shape of the cutoff matter and Lorentz invariance is not guaranteed;
   this is Nakanishi's objection one step earlier.
[obvious] O2.  No singularity in s on the real axis: nothing happens at s_min(|p|) or at 2 M^2, B_Lambda is smooth in s.

Predictions scored by the run (a = +-i, s_ref = -1 for the subtraction, Lambda in {10, 30, 100, 300, 1000}, dps 15):
  P1. The local slope c_j = [B(Lambda_{j+1}) - B(Lambda_j)] / ln(Lambda_{j+1}/Lambda_j) of the (+-) term tends to
      1/(4 pi^2 p0^2) with the FRAME'S p0^2 (rest: s; boost: s + p^2), relative error < 10^-3 at the last pair, in both frames,
      at s = 1.0 and 2.5, |p| = 0.8.  [from D2, D3]
  P2. The subtracted difference B_Lambda(s) - B_Lambda(s_ref) still grows like ln Lambda with slope (1/4 pi^2)(1/p0^2 - 1/p0_ref^2),
      p0_ref^2 = s_ref + p^2 in the frame: not UV finite.  [from D2]
  P3. R(Lambda) = B_Lambda - ln(Lambda)/(4 pi^2 p0^2) converges (|R(1000) - R(300)| < 10^-4) to a constant C that DIFFERS between
      the rest and boosted frames at the same s (reported, |C_boost - C_rest| > 10^-4).  [from O1; a guess as to size]
  P4. The real-axis (--) term equals bubble_lib.B_momentum(kmax = Lambda) for (-i, -i) to 10^-12, and the (++) term equals MINUS
      bubble_lib.B_momentum for (i, i).  [from D1]
  P5. Residue sum = direct real-axis k0 quadrature at three sample k_vec, both frames, to 10^-12.  [rule 6]
  P6. Smooth in s: on s in {1.5, 1.7, 1.9, 2.1} at Lambda = 100, |p| = 0.8, the second difference of B is small and of one sign
      (no kink at s_min(0.8) = 1.707 or at 2); reported, not scored.  [from O2]
Status PASS if P1, P2, P4, P5 hold; FAIL otherwise.  The nonlocal divergence is the EXPECTED outcome of this run; a UV-finite
subtracted difference would be the surprise and would say D2 is wrong.

RECORD OF THE SMOKE (rule 11; 3 Nov 00:10 ET, sandbox, dps 12, Lambda 10/30/100, s = 2.5; log results/smoke_sandbox_2026-11-03/
minkowski_smoke1.log).  The header above is kept as written.
  - P5, P4 PASS as written (residue sum = direct quadrature to 4e-13; four sample points, the header said three; (--) = Lee-Wick to 0; (++) = -Lee-Wick to 0).
  - Rest frame: slope 0.0101304 vs predicted 1/(4 pi^2 s) = 0.0101321, rel err 1.7e-4: D2 holds.
  - D3 IS WRONG.  Boosted slope 0.0101312 + 1e-6 i -> the SAME 1/(4 pi^2 s), not 1/(4 pi^2 (s + p^2)).  The tail expansion of D2
    dropped the O(1) term k_hat.p_vec in w1 - w2 in the boosted frame: w1 - w2 = x + (a1 - a2 - p^2 + x^2)/(2k) + ..., x = k_hat.p_vec,
    the leading term is odd in k_hat and averages out, and the 1/k term's angular average is (1/2) int_{-1}^{1} dc (A + p^2 c^2)/(A -
    p^2 c^2)^2 = 1/(A - p^2) = 1/s with A = p0^2 (the integrand is d/dc of c/(A - p^2 c^2)).  The coefficient of the nonlocal
    divergence is LORENTZ INVARIANT, 1/(4 pi^2 s); the frame dependence is in the finite part C (P3), smoke: C_rest = -0.000479,
    C_boost = -0.00127 + 0.00098 i with the corrected coefficient subtracted.  Prediction P1 is re-scored for the full run as:
    slope -> 1/(4 pi^2 s) in BOTH frames (same tolerance); its first wording stays above.
  - P2 "FAIL" was a script error with a small finding inside: the subtraction point s_ref = -1 (Euclidean p0 = i, as in bubble_frame)
    is not available to the real-axis formulation -- the pole selection by sign(Im) is not the continuation of the real-p0 one
    (three poles move to the upper half plane) and B_ref grew like Lambda.  The real-axis integral has no Euclidean region to subtract
    at; subtraction is now at a REAL point s_ref = 0.5 (p0 > |p| in both frames).  P2 re-scored as: slope of B(s) - B(s_ref) ->
    (1/4 pi^2)(1/s - 1/s_ref) in both frames, rel err < 10^-3.
  - Timing: the 150 s per boosted point of the first smoke was the ill-defined Euclidean B_ref; with the real s_ref a boosted integral
    is 15 s at dps 12 and 28 s at dps 15, Lambda = 300 (sandbox).  B_ref is computed once per frame, Lambda stops at 300, P6 at
    Lambda = 30.  Expected 5-10 min (rule 16 holds).  Second smoke (dps 12, Lambda 10/30/100): status=PASS as re-scored,
    results/smoke_sandbox_2026-11-03/minkowski_smoke2.log.

AFTER THE READ (rule 28, 3 Nov 17:50 ET; the author's full run, results/minkowski.log, status=PASS, compared with the sandbox smoke
before the paper was opened: boost numbers identical, rest s = 2.5 differing at 1e-6 -- the window bug in B_rest, found by the review).  1612.06510 read 3 Nov 17:45 (EPJC 77 (2017) 84; a fetch-summary of secs. 1-3, 5 and
app. A, not yet the full text).  Their formulation IS this one: loop energy along the real axis, residues closing in the upper half
plane, hard cutoff on the space momenta, rest frame.  Their propagator is the physical pole TIMES the pair, S = 1/(p^2 - m^2 + i eps) *
M^4/((p^2)^2 + M^4), poles at p^2 = +-i M^2 (our a = +-i), partial fractions 1/x - (1/2)/(x - i) - (1/2)/(x + i) at m = 0.
  D1  their power counting, "a propagator calculated on the pole of another falls off half less rapidly": the same statement.  MATCH.
  D2  not in the paper for the bare pair.  Their D = 4 trivial-numerator coefficient is c_00 = 0: the physical pole's cross terms
      cancel the nonlocal log.  P7 below checks that with our building blocks (prediction from the D2 tail with the (-+) selection:
      (0, +i) has tail i(a1 - a2)/(4 p0^2) = +1/4 per unit, (+i, -i) has 1/2; weights 2(-1/2) and 1/2: 1/4 * (-1) + 1/2 * 1/2 = 0).
  D3  WRONG as written (smoke record); they take p = (0, p_s) and "find (3.2) again": the divergent part is a function of p^2 only.
      Our corrected 1/(4 pi^2 s) in both frames AGREES with them.
  P3  the frame dependence of the FINITE part is not discussed in the paper (they state Sigma "is expected to be a function of p^2
      only" and check it for the divergent part).  Candidate "what the era did not write"; to be confirmed on the full text.
  D4  the (++)/(--) sign flip and the cancellation of the local log in the four-pair sum: not in the paper in that form; their
      "the usual divergences are anti-Hermitian, these are not" is the same fact seen from the counterterm side.
  M1  accurate and shallower than the derived lines.
  "No Euclidean region to subtract at": their sec. 5, "the analytic continuation to Euclidean space is not possible".  MATCH.
  Their D = 6 headline, Sigma = -M^4/(2 (4 pi)^3) [M^2/(p^2)^2 - i/p^2] ln(Lambda^2/M^2) + ..., is the same mechanism one dimension up
  (two nonlocal structures there, one here).  Reading note for the ledger (author, 17:50): the large-k_s "pinch", Omega ~ p_s - i eta
  and Omega-bar ~ p_s + i eta closing on the real k0 axis from both sides, is a conjugate pair squeezing the contour, not a pole
  collision -- the same misreading hdym's Krein pairs invited first.

  P7 (added 3 Nov 17:55, author's yes; scored).  The weighted six-pair sum of real-axis bubbles with S's partial-fraction weights
      (1, -1/2, -1/2; the massless physical pole as a1 = -i*10^-10, Feynman) has log-slope |c| < 10^-5 at s = 2.5, rest frame,
      from Lambda = 30 to 300 (c_00 = 0 of 1612.06510 eq. 3.4-3.6), while the bare pair's is 1/(4 pi^2 s).  Sandbox 3 Nov 17:48 (dps 12):
      total -3.6e-7 + 9.1e-7 i (the finite-Lambda remainder, O(Lambda^-2) by D2), bare pair 0.0101311; m = 0 only -- the m != 0
      entries of their table have m-dependent weights and are not checked here.

Run (from bubble/):   python3 -B -u bubble_minkowski.py --json results/minkowski.json > results/minkowski.log
Smoke (sandbox):      python3 -B -u bubble_minkowski.py --smoke
"""
import sys, time, json, argparse
import mpmath as mp
import bubble_lib as L

ap = argparse.ArgumentParser()
ap.add_argument('--dps', type=int, default=15)
ap.add_argument('--json', default=None)
ap.add_argument('--smoke', action='store_true', help='fewer cutoffs, one s, lower dps: minutes, not a result')
ap.add_argument('--lams', default='10,30,100,300')
ap.add_argument('--svals', default='1.0,2.5')
ap.add_argument('--p', type=float, default=0.8)
args = ap.parse_args()
if args.smoke:
    args.dps = min(args.dps, 12); args.lams = '10,30,100'; args.svals = '2.5'
mp.mp.dps = args.dps
t0 = time.time()
I = L.I; PI = L.PI
a_plus, a_minus = mp.mpc(0, 1), mp.mpc(0, -1)
s_ref = mp.mpf('0.5')   # REAL subtraction point (smoke record: the Euclidean point is not available to the real-axis integral)
p = mp.mpf(args.p)
lams = [mp.mpf(x) for x in args.lams.split(',')]
svals = [mp.mpf(x) for x in args.svals.split(',')]
out = {'dps': args.dps, 'p': float(p), 'lams': [float(x) for x in lams], 'svals': [float(x) for x in svals], 'checks': {}}
print(f"# bubble_minkowski: real-axis k0, sphere cutoff in the frame; a = +-i, D = 4, s_ref = {s_ref}, |p| = {p}, dps = {args.dps}, "
      f"Lambda = {[float(x) for x in lams]}, s = {[float(x) for x in svals]}")


def h_real(w1, w2, p0):
    """int dk0/(2 pi) 1/((k0^2 - w1^2)((k0 - p0)^2 - w2^2)) along the REAL axis = i * (sum of residues in the upper half plane).
    Poles and residues: +w1: 1/(2 w1 ((w1-p0)^2 - w2^2)); -w1: -1/(2 w1 ((w1+p0)^2 - w2^2));
    p0+w2: 1/(2 w2 ((p0+w2)^2 - w1^2)); p0-w2: -1/(2 w2 ((p0-w2)^2 - w1^2)).  Selected by sign(Im) -- generic, no hand formula."""
    poles = [(w1, 1 / (2 * w1 * ((w1 - p0) ** 2 - w2 ** 2))),
             (-w1, -1 / (2 * w1 * ((w1 + p0) ** 2 - w2 ** 2))),
             (p0 + w2, 1 / (2 * w2 * ((p0 + w2) ** 2 - w1 ** 2))),
             (p0 - w2, -1 / (2 * w2 * ((p0 - w2) ** 2 - w1 ** 2)))]
    tot = mp.mpc(0)
    for z, r in poles:
        if mp.im(z) == 0:
            raise ValueError("pole on the real k0 axis")
        if mp.im(z) > 0:
            tot += r
    return I * tot


def h_direct(w1, w2, p0):
    """the same integral by quadrature along the real k0 axis (rule 6 check)"""
    f = lambda k0: 1 / ((k0 * k0 - w1 * w1) * ((k0 - p0) ** 2 - w2 * w2))
    pts = sorted(set([-mp.inf, mp.mpf(0), mp.re(p0), mp.re(w1), -mp.re(w1), mp.re(p0 + w2), mp.re(p0 - w2), mp.inf]))
    return mp.quad(f, pts) / (2 * PI)


def omegas_rest(k, a1, a2):
    return mp.sqrt(k * k + a1), mp.sqrt(k * k + a2)


def omegas_boost(kpar, kperp, a1, a2):
    return mp.sqrt(kpar * kpar + kperp * kperp + a1), mp.sqrt((kpar - p) ** 2 + kperp * kperp + a2)


def B_rest(p0, a1, a2, lam):
    """int_{|k|<lam} d^3k/(2 pi)^3 h_real, rest frame.  Breakpoints where Re(w1 + w2) = p0 or Re(w1 - w2) = +-p0 could make the
    denominators small (only the (++)/(--) terms have a near-zero of (w1+w2)^2 - p0^2; the mixed terms' p0^2 - (w1-w2)^2 is
    bounded away from zero in the rest frame since w1 - w2 is purely imaginary there)."""
    f = lambda k: k * k * h_real(*omegas_rest(k, a1, a2), p0)
    pts = [mp.mpf(0), mp.mpf(1)]
    # windows only for the terms whose denominator is (w1+w2)^2 - p0^2, i.e. Im a1 and Im a2 of one sign ((++), (--), (0,-i)); the
    # mixed terms have p0^2 - (w1-w2)^2, bounded away from zero.  (Review of 3 Nov, must-fix: applied to the mixed term too, the
    # forced 1e-3 windows at a real k_t made the rest-frame s = 2.5 values dps-dependent at 7e-7 -- C_rest(2.5) off by 1.5e-3 relative
    # in the author's first log; slopes and every PASS/FAIL unaffected, the error being Lambda-independent.)
    kt2 = (p0 ** 4 - 2 * p0 ** 2 * (a1 + a2) + (a1 - a2) ** 2) / (4 * p0 ** 2)
    kt = mp.sqrt(kt2)
    if mp.sign(mp.im(a1)) == mp.sign(mp.im(a2)) and mp.re(kt) > 0 and 0 < abs(mp.im(kt)) < 1:
        r, im = mp.re(kt), abs(mp.im(kt))
        pts += [r - 4 * im, r - im, r, r + im, r + 4 * im]
    pts += [lam]
    pts = sorted(set(x for x in pts if 0 <= x <= lam))
    return mp.quad(f, pts) / (2 * PI ** 2)


def B_boost(p0, a1, a2, lam):
    """int_{|k|<lam} d^3k/(2 pi)^3 h_real, boosted frame, cylindrical (kpar, kperp), sphere cutoff kpar^2 + kperp^2 < lam^2.
    For the mixed terms |Re(w1 - w2)| <= ~|p| < p0 so p0^2 - (w1 - w2)^2 stays away from zero: no windows needed."""
    def inner(kpar):
        top = mp.sqrt(lam * lam - kpar * kpar)
        f = lambda kperp: kperp * h_real(*omegas_boost(kpar, kperp, a1, a2), p0)
        pts = sorted(set(x for x in [mp.mpf(0), mp.mpf(1) / 2, mp.mpf(2), top] if 0 <= x <= top))
        return mp.quad(f, pts)
    pts = sorted(set([-lam, -mp.mpf(2), mp.mpf(0), p / 2, p, p + 2, lam]))
    return mp.quad(inner, pts) / (4 * PI ** 2)


def slope(vals):
    return [(vals[j + 1] - vals[j]) / mp.log(lams[j + 1] / lams[j]) for j in range(len(vals) - 1)]


def fmt(z):
    return mp.nstr(z, 12)


ok = True
# ---------------------------------------------------------------- P5: residue sum vs direct real-axis quadrature
print("\n## P5: residue sum vs direct real-axis k0 quadrature (rule 6)")
worst = mp.mpf(0)
for (w1, w2, p0, lab) in [(*omegas_rest(mp.mpf('0.7'), a_plus, a_minus), mp.sqrt(mp.mpf('2.5')), 'rest k=0.7 (+-)'),
                          (*omegas_rest(mp.mpf('3.0'), a_plus, a_plus), mp.sqrt(mp.mpf('2.5')), 'rest k=3.0 (++)'),
                          (*omegas_boost(mp.mpf('0.4'), mp.mpf('1.1'), a_plus, a_minus), mp.sqrt(mp.mpf('2.5') + p * p), 'boost (0.4,1.1) (+-)'),
                          (*omegas_boost(mp.mpf('-0.9'), mp.mpf('0.2'), a_minus, a_plus), mp.sqrt(mp.mpf('1.0') + p * p), 'boost (-0.9,0.2) (-+)')]:
    hr, hd = h_real(w1, w2, p0), h_direct(w1, w2, p0)
    d = abs(hr - hd); worst = max(worst, d)
    print(f"  {lab:24s} residues {fmt(hr)}  direct {fmt(hd)}  |diff| {mp.nstr(d, 3)}")
p5 = worst < mp.mpf(10) ** (-12 if args.dps >= 15 else -(args.dps - 3))
print(f"  P5 {'PASS' if p5 else 'FAIL'} (worst {mp.nstr(worst, 3)})"); ok &= p5
out['checks']['P5'] = {'worst': float(worst), 'pass': bool(p5)}

# ---------------------------------------------------------------- P4: (--) and (++) against the Lee-Wick contour value
print("\n## P4: real-axis (--) vs bubble_lib.B_momentum (Lee-Wick contour) with the same cutoff; (++) vs MINUS it  [rest frame]")
s4 = mp.mpf('2.5'); p0 = mp.sqrt(s4); lam4 = mp.mpf(30)
mm_real = B_rest(p0, a_minus, a_minus, lam4); mm_lw = L.B_momentum(4, p0, a_minus, a_minus, kmax=lam4)
pp_real = B_rest(p0, a_plus, a_plus, lam4); pp_lw = L.B_momentum(4, p0, a_plus, a_plus, kmax=lam4)
d_mm, d_pp = abs(mm_real - mm_lw), abs(pp_real + pp_lw)
print(f"  (--): real-axis {fmt(mm_real)}  LW {fmt(mm_lw)}  |diff| {mp.nstr(d_mm, 3)}")
print(f"  (++): real-axis {fmt(pp_real)}  LW {fmt(pp_lw)}  |real + LW| {mp.nstr(d_pp, 3)}   (|real - LW| = {mp.nstr(abs(pp_real - pp_lw), 3)})")
p4 = max(d_mm, d_pp) < mp.mpf(10) ** (-10 if args.dps >= 15 else -(args.dps - 4))
print(f"  P4 {'PASS' if p4 else 'FAIL'}"); ok &= p4
out['checks']['P4'] = {'mm_real': [float(mp.re(mm_real)), float(mp.im(mm_real))], 'mm_lw': [float(mp.re(mm_lw)), float(mp.im(mm_lw))],
                       'pp_real': [float(mp.re(pp_real)), float(mp.im(pp_real))], 'pp_lw': [float(mp.re(pp_lw)), float(mp.im(pp_lw))], 'pass': bool(p4)}

# ---------------------------------------------------------------- P1-P3: the mixed term, cutoff sweep, both frames
print("\n## P1-P3: the (+-) term, B_Lambda vs ln Lambda, rest and boosted frame (predicted slope, re-scored after the smoke: 1/(4 pi^2 s) in both frames)")
res = {}
ref = {}
for frame in ('rest', 'boost'):
    Bf = B_rest if frame == 'rest' else B_boost
    p0ref_sq = s_ref + (p * p if frame == 'boost' else 0)
    ref[frame] = []
    for lam in lams:
        t1 = time.time()
        vr = Bf(mp.sqrt(p0ref_sq), a_plus, a_minus, lam)
        ref[frame].append(vr)
        print(f"  s_ref={float(s_ref)} {frame:5s} Lambda={float(lam):<7} B_ref={fmt(vr)}  [{time.time() - t1:.1f} s]", flush=True)
for s in svals:
    for frame in ('rest', 'boost'):
        p0sq = s + (p * p if frame == 'boost' else 0)
        p0 = mp.sqrt(p0sq)
        pred = 1 / (4 * PI ** 2 * s)                       # re-scored after the smoke (header record); the first wording had p0sq
        pred_sub = (1 / (4 * PI ** 2)) * (1 / s - 1 / s_ref)
        Bf = B_rest if frame == 'rest' else B_boost
        vals, vals_ref = [], ref[frame]
        for lam in lams:
            t1 = time.time()
            v = Bf(p0, a_plus, a_minus, lam)
            vals.append(v)
            print(f"  s={float(s):<4} {frame:5s} Lambda={float(lam):<7} B={fmt(v)}  [{time.time() - t1:.1f} s]", flush=True)
        sl = slope(vals); sl_sub = slope([v - vr for v, vr in zip(vals, vals_ref)])
        R = [v - mp.log(lam) * pred for v, lam in zip(vals, lams)]
        err1 = abs(sl[-1] / pred - 1); err2 = abs(sl_sub[-1] / pred_sub - 1)
        print(f"    slopes      {[fmt(x) for x in sl]}   predicted {fmt(pred)}   rel err last {mp.nstr(err1, 3)}")
        print(f"    slopes(sub) {[fmt(x) for x in sl_sub]}   predicted {fmt(pred_sub)}   rel err last {mp.nstr(err2, 3)}")
        print(f"    R(Lambda) = B - ln(Lambda) pred: {[fmt(x) for x in R]}   |R_last - R_prev| = {mp.nstr(abs(R[-1] - R[-2]), 3)}")
        res[(float(s), frame)] = dict(p0sq=p0sq, pred=pred, pred_sub=pred_sub, slopes=sl, slopes_sub=sl_sub, R=R, err1=err1, err2=err2,
                                      B=vals, Bref=vals_ref)
tol1 = mp.mpf(10) ** (-3 if not args.smoke else -1)
p1 = all(r['err1'] < tol1 for r in res.values())
p2 = all(r['err2'] < tol1 for r in res.values())
print(f"  P1 {'PASS' if p1 else 'FAIL'} (slope -> 1/(4 pi^2 s) in both frames [re-scored], rel err < {mp.nstr(tol1, 2)})"); ok &= p1
print(f"  P2 {'PASS' if p2 else 'FAIL'} (B(s) - B(s_ref) still diverges, slope (1/4 pi^2)(1/s - 1/s_ref) [re-scored])"); ok &= p2
print("  P3 (reported): C = lim R(Lambda), rest vs boost at the same s:")
for s in svals:
    Cr, Cb = res[(float(s), 'rest')]['R'][-1], res[(float(s), 'boost')]['R'][-1]
    conv = max(abs(res[(float(s), 'rest')]['R'][-1] - res[(float(s), 'rest')]['R'][-2]), abs(res[(float(s), 'boost')]['R'][-1] - res[(float(s), 'boost')]['R'][-2]))
    print(f"    s={float(s)}: C_rest={fmt(Cr)}  C_boost={fmt(Cb)}  |diff|={mp.nstr(abs(Cr - Cb), 4)}  (convergence of R: {mp.nstr(conv, 2)})")
out['P123'] = {f"{k[0]}_{k[1]}": {kk: ([[float(mp.re(x)), float(mp.im(x))] for x in vv] if isinstance(vv, list) else [float(mp.re(vv)), float(mp.im(vv))])
                                  for kk, vv in r.items()} for k, r in res.items()}
out['checks']['P1'] = bool(p1); out['checks']['P2'] = bool(p2)

# ---------------------------------------------------------------- D4: the four-pair sum, rest frame
print("\n## D4 (reported): the symmetric four-pair sum, rest frame, s = 2.5: slopes of each term")
if not args.smoke:
    p0 = mp.sqrt(mp.mpf('2.5')); terms = {}
    for lab, (a1, a2) in {'++': (a_plus, a_plus), '--': (a_minus, a_minus), '+-': (a_plus, a_minus), '-+': (a_minus, a_plus)}.items():
        terms[lab] = [B_rest(p0, a1, a2, lam) for lam in lams]
        print(f"  ({lab}) slopes {[fmt(x) for x in slope(terms[lab])]}")
    tot = [sum(terms[l][j] for l in terms) for j in range(len(lams))]
    print(f"  sum  slopes {[fmt(x) for x in slope(tot)]}   predicted 2/(4 pi^2 p0^2) = {fmt(2 / (4 * PI ** 2 * p0 ** 2))}")
    out['D4'] = {l: [[float(mp.re(x)), float(mp.im(x))] for x in v] for l, v in terms.items()}
else:
    print("  (skipped in smoke)")

# ---------------------------------------------------------------- P6: smoothness in s
print("\n## P6 (reported): B_Lambda(s) across s_min(0.8) = 1.707 and 2 M^2 = 2, boosted frame, Lambda = 30, (+-)")
if not args.smoke:
    lam6 = mp.mpf(30); sgrid = [mp.mpf(x) for x in ('1.5', '1.7', '1.9', '2.1')]
    b6 = [B_boost(mp.sqrt(s + p * p), a_plus, a_minus, lam6) for s in sgrid]
    d2 = [b6[j + 1] - 2 * b6[j] + b6[j - 1] for j in range(1, len(b6) - 1)]
    print(f"  B: {[fmt(x) for x in b6]}\n  second differences: {[fmt(x) for x in d2]}")
    out['P6'] = [[float(mp.re(x)), float(mp.im(x))] for x in b6]
else:
    print("  (skipped in smoke)")

# ---------------------------------------------------------------- P7: the partial-fraction cancellation (1612.06510 c_00 = 0), rest frame
print("\n## P7: S x S with S = 1/x - (1/2)/(x - i) - (1/2)/(x + i) (1612.06510 eq. 2.2 at m = 0): weighted six-pair sum, log-slope should vanish")
p0 = mp.sqrt(mp.mpf('2.5')); l7 = [mp.mpf(30), mp.mpf(300)]
eps7 = mp.mpc(0, -1e-10)          # massless Feynman propagator, m^2 - i eps
S7 = {}
for lab, (a1, a2) in {'(0,0)': (eps7, eps7), '(0,+i)': (eps7, a_plus), '(0,-i)': (eps7, a_minus), '(+i,+i)': (a_plus, a_plus),
                      '(-i,-i)': (a_minus, a_minus), '(+i,-i)': (a_plus, a_minus)}.items():
    v = [B_rest(p0, a1, a2, l) for l in l7]
    S7[lab] = (v[1] - v[0]) / mp.log(l7[1] / l7[0])
    print(f"  {lab:8s} log-slope {fmt(S7[lab])}")
tot7 = S7['(0,0)'] - (S7['(0,+i)'] + S7['(0,-i)']) + (S7['(+i,+i)'] + S7['(-i,-i)']) / 4 + S7['(+i,-i)'] / 2
nl7 = -S7['(0,+i)'] + S7['(+i,-i)'] / 2
print(f"  weighted sum {fmt(tot7)}   (nonlocal part alone {fmt(nl7)}; bare pair {fmt(S7['(+i,-i)'])}, 1/(4 pi^2 s) = {fmt(1 / (4 * PI ** 2 * p0 ** 2))})")
p7 = abs(tot7) < mp.mpf(10) ** -5
print(f"  P7 {'PASS' if p7 else 'FAIL'} (|weighted log-slope| < 1e-5: c_00 = 0)"); ok &= p7
out['checks']['P7'] = {'slopes': {k: [float(mp.re(v)), float(mp.im(v))] for k, v in S7.items()}, 'total': [float(mp.re(tot7)), float(mp.im(tot7))], 'pass': bool(p7)}

print(f"\nstatus={'PASS' if ok else 'FAIL'}  (P1 {p1}, P2 {p2}, P4 {p4}, P5 {p5}, P7 {p7})")
if args.json:
    with open(args.json, 'w') as f:
        json.dump(out, f, indent=1)
print(f"END bubble_minkowski {time.time() - t0:.0f} s")
