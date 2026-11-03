#!/usr/bin/env python3
"""
Rung 2, the next technique (1969): the CLOP prescription (Cutkosky, Landshoff, Olive & Polkinghorne, Nucl. Phys. B 12, 281,
unseen; as described in Anselmi & Piva 1703.04584 section 4, read 2 Nov 2026) applied to the complex-ghost-pair bubble of
bubble_frame.py -- shift the Lee-Wick scale in ONE propagator, M -> M + delta, compute, let delta -> 0 from both sides.  Here the
pair is a1 = i M^2, a2 = -i M^2 (M = 1), and the shift acts on a1 -> i (1 + delta)^2 only (--shift 2 shifts a2 instead).

Why: Lee & Wick's reply to Nakanishi (Phys. Rev. D 3, 1046 (1971), unseen; hep-th/0609206 section 1) attributed the frame
dependence to real spatial momenta and reached for the CLOP rule for the coalescing singularities.  This run asks what the mass-shift
limit does to the SAME integral, in each method and each frame, before anything later is used.

Predictions (2 Nov 2026, written after bubble_frame.py's first run and before this script's first run; a = +-i, D = 4, UV subtraction at s = -1 in the same frame as in bubble_frame.py;
delta in {+-0.1, +-0.03, +-0.01}; s in {2.5, 4.0} above the rest-frame threshold 2 M^2, s = 1.0 below it):
  P1. Rest frame, method M.  With delta != 0 the pair is no longer complex conjugate: w1 + w2 is complex for all real k with
      Im(w1 + w2) = delta Re(1/w1) + O(delta^2) of DEFINITE sign, so the pinch is lifted to one side and the delta -> 0+ and 0- limits
      are the two i-epsilon values, B(+i eps) and B(-i eps) of bubble_lib (the sign assignment is printed, not predicted) -- NOT the
      principal value.  On a single (M, M*) pair CLOP's ambiguity is the full would-be discontinuity 2 f_d.  Scored: |B(delta) - B(sign
      i eps)| decreasing with |delta| for the matching sign, and the two limits differing by |B(+i eps) - B(-i eps)| within 10%.
  P2. Rest frame, method P (Feynman parameters).  With delta != 0 the real-axis branch point of the unshifted pair moves off the
      real t axis, the upper and lower continuations coincide (the path's Im t = +-10^-(dps-6) is far inside |Im t*| ~ delta), and the
      delta -> 0+- limits are f_+ and f_- of the unshifted pair.  Methods P = M for every delta (rule 3): the shifted pair has no
      pinch, so the two methods compute the same unambiguous function.
  P3. The full Lee-Wick propagator (both +-i M^2 poles on each line; the symmetric four-pair sum J = (1/4) sum_{a,b} B(a, b) of
      bubble_checks --lw, with the shift applied to line 1's pair): the (+-) and (-+) cross terms acquire OPPOSITE signs of
      Im(w1 + w2), their i pi residue terms cancel, and the two CLOP limits coincide with the principal-value (average) J: no ambiguity
      for the equal-mass pair, as 1703.04584 section 4 states on the 2 Nov reading of it (their +-pi survives only for unequal physical masses -- a
      2017-model statement, left for the back sweep).  Scored: |J(delta) - J_PV| decreasing with |delta| for both signs.
  P4. Boosted frame (|p| = 0.8), method M: there is no pinch at p_vec != 0 (bubble_frame.py), the shift is a regular perturbation, and
      delta -> 0 from either side reproduces the frame-dependent real-momentum integral of bubble_frame.py: CLOP does not restore
      Lorentz invariance of the Hamiltonian (real-k) integral.  Scored: |B_boost(delta) - B_boost(0)| -> 0 from both sides, with
      B_boost(0) != B_rest at s = 2.5.
  P5 (reading, not scored). Method P is Lorentz invariant by construction, so "CLOP + method P" is frame independent while "CLOP +
      method M" is not: the Lee-Wick reply restores Lorentz invariance only by replacing the real-momentum integral with the covariant
      parametric one plus a rule -- Nakanishi's objection, in one sentence.
Status PASS if P1-P4 hold.

RECORD (rule 11; 2 Nov evening).  P4's premise "there is no pinch at p_vec != 0" is wrong in the same way as bubble_frame.py's (see
its header): the pinch is a circle in the plane k_par = p/2 above s_min(|p|); the integral is conditionally convergent with the
i epsilon supplied by the geometry.  The prediction's content -- the shift is a regular perturbation and delta -> 0 from either side
returns the frame-dependent real-momentum integral -- is unchanged; the first smoke (pointwise quadrature, PASS, kept in
results/superseded/clop_sandbox_pointwise_quadrature.log) had P4's boosted numbers unconverged (the delta = 0 point by 0.5%, the
delta != 0 points by 1-3% between breakpoint sets, since the shifted pair's near-pole was not windowed either); re-run with
bubble_lib.B_momentum_boost, which windows the near-pole for non-conjugate pairs too (second review finding, 2 Nov).
"""
import time, json, argparse
import mpmath as mp
import bubble_lib as L

ap = argparse.ArgumentParser()
ap.add_argument('--dps', type=int, default=15)
ap.add_argument('--deltas', default='0.1,0.03,0.01')
ap.add_argument('--shift', type=int, default=1, choices=[1, 2])
ap.add_argument('--pboost', default='0.8')
ap.add_argument('--json', default=None)
args = ap.parse_args()
mp.mp.dps = args.dps
t0 = time.time()
I = L.I; PI = L.PI
s_ref = mp.mpf(-1)
tol = mp.mpf(10) ** (-(args.dps - 7))
deltas = [mp.mpf(d) for d in args.deltas.split(',')]
print(f"# bubble_clop: CLOP mass-shift limit on the complex-ghost-pair bubble; a = +-i, shift on line {args.shift}, deltas +-{args.deltas}, D = 4, s_ref = {s_ref}, dps = {args.dps}")


def pair(delta):
    """the (M, M*) pair with the shift M -> M + delta on one line: a = i (1+delta)^2 for the shifted line"""
    a1 = I * (1 + delta) ** 2 if args.shift == 1 else I
    a2 = -I * (1 + delta) ** 2 if args.shift == 2 else -I
    return a1, a2


def B_rest(s, a1, a2, prescription='PV'):
    s = mp.mpf(s)
    p0 = mp.sqrt(s) if s > 0 else I * mp.sqrt(-s)
    return L.B_momentum(4, p0, a1, a2, prescription=prescription, p0_ref=I)


def B_P(s, a1, a2, upper=True):
    f, info = L.continuation(4, mp.mpf(s), a1, a2, upper=upper)
    fref = L.I_straight(4, s_ref, a1, a2)
    return -I / (16 * PI ** 2) * (f - fref), info


def B_boost(s, p, a1, a2):
    return L.B_momentum_boost(s, p, a1, a2, s_ref)


ok = True
out = {}
a1_0, a2_0 = pair(mp.mpf(0))

print("# --- P1: rest frame, method M, the mass-shift limit against PV and the two i-epsilon values ---")
for s in (mp.mpf('2.5'), mp.mpf('4.0')):
    pv = B_rest(s, a1_0, a2_0, 'PV'); fey = B_rest(s, a1_0, a2_0, 'feynman'); mnu = B_rest(s, a1_0, a2_0, 'minus')
    print(f"#   s = {mp.nstr(s, 3)}: PV = {mp.nstr(pv, 10)}  +i eps = {mp.nstr(fey, 10)}  -i eps = {mp.nstr(mnu, 10)}  |gap| = {mp.nstr(abs(fey - mnu), 5)}")
    lim = {}
    for sign in (+1, -1):
        prev = None; mono = True
        for d in deltas:
            t1 = time.time()
            a1, a2 = pair(sign * d)
            b = B_rest(s, a1, a2, 'PV')          # no pinch for the shifted pair: prescription irrelevant
            dF, dM, dP = abs(b - fey), abs(b - mnu), abs(b - pv)
            print(f"#   delta = {'+' if sign > 0 else '-'}{mp.nstr(d, 3)}: B = {mp.nstr(b, 10)}  |B - (+i eps)| = {mp.nstr(dF, 4)}  |B - (-i eps)| = {mp.nstr(dM, 4)}  |B - PV| = {mp.nstr(dP, 4)} ({time.time()-t1:.0f}s)")
            near = min(dF, dM)
            if prev is not None:
                mono &= near < prev
            prev = near
            lim[sign] = b
        which = '+i eps' if abs(lim[sign] - fey) < abs(lim[sign] - mnu) else '-i eps'
        print(f"#   delta -> 0{'+' if sign > 0 else '-'} approaches {which}; approach monotone: {mono}")
        ok &= mono
    gap_ratio = abs(lim[+1] - lim[-1]) / abs(fey - mnu)
    print(f"#   |lim+ - lim-| / |(+i eps) - (-i eps)| = {mp.nstr(gap_ratio, 5)}   [P1: within 10% of 1 at the smallest delta]")
    ok &= abs(gap_ratio - 1) < mp.mpf('0.1')
    out[f'P1_s{float(s)}'] = dict(limp=str(lim[+1]), limm=str(lim[-1]), pv=str(pv), fey=str(fey), mnu=str(mnu))

print("# --- P2: rest frame, method P with the shifted pair: upper = lower continuation; delta -> 0+- against f_+, f_- ---")
for s in (mp.mpf('2.5'),):
    fav, fd, fp, fm, _ = L.average_continuation(4, s, a1_0, a2_0)
    fref0 = L.I_straight(4, s_ref, a1_0, a2_0)
    Bp = -I / (16 * PI ** 2) * (fp - fref0); Bm = -I / (16 * PI ** 2) * (fm - fref0)
    print(f"#   s = {mp.nstr(s, 3)} unshifted: B(f_+) = {mp.nstr(Bp, 10)}  B(f_-) = {mp.nstr(Bm, 10)}")
    for sign in (+1, -1):
        prev = None; mono = True
        for d in deltas:
            t1 = time.time()
            a1, a2 = pair(sign * d)
            bu, iu = B_P(s, a1, a2, True); bl, il = B_P(s, a1, a2, False)
            bm = B_rest(s, a1, a2, 'PV')
            same = abs(bu - bl)
            pm = abs(bu - bm)
            dp_, dm_ = abs(bu - Bp), abs(bu - Bm)
            print(f"#   delta = {'+' if sign > 0 else '-'}{mp.nstr(d, 3)}: P_upper = {mp.nstr(bu, 10)}  |upper - lower| = {mp.nstr(same, 3)} (ints {iu['ints']}/{il['ints']})  |P - M| = {mp.nstr(pm, 3)}  |P - B(f_+)| = {mp.nstr(dp_, 4)}  |P - B(f_-)| = {mp.nstr(dm_, 4)} ({time.time()-t1:.0f}s)")
            ok &= same < tol and pm < tol
            near = min(dp_, dm_)
            if prev is not None:
                mono &= near < prev
            prev = near
        print(f"#   delta -> 0{'+' if sign > 0 else '-'}: approach to one of f_+, f_- monotone: {mono}")
        ok &= mono

print("# --- P3: the full Lee-Wick propagator (four-pair sum), shift on line 1: the two limits against the PV sum ---")
for s in (mp.mpf('2.5'),):
    def Jsum(delta):
        a1s = [I * (1 + delta) ** 2, -I * (1 + delta) ** 2]
        return sum(L.B_momentum(4, mp.sqrt(s), a, b, 'PV', p0_ref=I) for a in a1s for b in (I, -I)) / 4
    jpv = Jsum(mp.mpf(0))
    print(f"#   s = {mp.nstr(s, 3)}: J_PV (unshifted, principal value on the cross terms) = {mp.nstr(jpv, 10)}")
    for sign in (+1, -1):
        prev = None; mono = True
        for d in deltas:
            t1 = time.time()
            j = Jsum(sign * d); dj = abs(j - jpv)
            print(f"#   delta = {'+' if sign > 0 else '-'}{mp.nstr(d, 3)}: J = {mp.nstr(j, 10)}  |J - J_PV| = {mp.nstr(dj, 4)} ({time.time()-t1:.0f}s)")
            if prev is not None:
                mono &= dj < prev
            prev = dj
        print(f"#   delta -> 0{'+' if sign > 0 else '-'}: approach to J_PV monotone: {mono}   [P3: both signs]")
        ok &= mono

print(f"# --- P4: boosted frame |p| = {args.pboost}, method M: the shift is a regular perturbation of the real-momentum integral ---")
p = mp.mpf(args.pboost); s = mp.mpf('2.5')
t1 = time.time()
b0 = B_boost(s, p, a1_0, a2_0)
brest = B_rest(s, a1_0, a2_0, 'PV')
print(f"#   s = {mp.nstr(s, 3)}: B_boost(delta = 0) = {mp.nstr(b0, 10)}  B_rest(PV) = {mp.nstr(brest, 10)}  |diff| = {mp.nstr(abs(b0 - brest), 4)} (frame dependence, bubble_frame.py) ({time.time()-t1:.0f}s)")
ok &= abs(b0 - brest) > tol * 100
for sign in (+1, -1):
    prev = None; mono = True
    for d in deltas:
        t1 = time.time()
        a1, a2 = pair(sign * d)
        b = B_boost(s, p, a1, a2); db = abs(b - b0)
        print(f"#   delta = {'+' if sign > 0 else '-'}{mp.nstr(d, 3)}: B_boost = {mp.nstr(b, 10)}  |B - B_boost(0)| = {mp.nstr(db, 4)} ({time.time()-t1:.0f}s)")
        if prev is not None:
            mono &= db < prev
        prev = db
    print(f"#   delta -> 0{'+' if sign > 0 else '-'}: approach to the unshifted boosted integral monotone: {mono}   [P4]")
    ok &= mono

if args.json:
    json.dump(out, open(args.json, 'w'), indent=1)
print(f"# END bubble_clop elapsed={time.time()-t0:.1f}s status={'PASS' if ok else 'FAIL'}")
