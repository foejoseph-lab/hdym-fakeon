#!/usr/bin/env python3
"""
Rung 2, technique of the day (1971): the complex-ghost-pair bubble at p_vec != 0 with the Lee-Wick energy contour and REAL spatial
loop momenta -- Nakanishi's computation (Phys. Rev. D 3, 811 (1971), unseen; restated by Nakanishi in hep-th/0609206 eqs. 3.3-3.6,
read 2 Nov 2026), the Lorentz-invariance crack as a number.  Setup and what counts as agreement:
notes/2026-11-02_rung2_1971_crack_summaries.md sections 1 and 4.

Model: the complex-conjugate pair a1 = M^2, a2 = conj(M^2) of bubble_lib (default a = +-i, the Lee-Wick pair of 1801.00915 eq. 5.5,
so the p_vec = 0 numbers are the committed results/lw.json; --gamma g switches to Nakanishi's M^2 = 1 + i g).  D = 4, metric
(+,-,-,-), B(p) = int d^4k/(2 pi)^4 1/((k^2 - a1)((k-p)^2 - a2)), UV-subtracted at the Euclidean point s_ref = p^2 = -1 IN THE
SAME FRAME (p0^2 -> s_ref + p_vec^2), which is Lorentz invariant below all thresholds.

Method M, boosted frame: the k0 integral along the Lee-Wick contour by residues (closing on the two right-hand poles, exactly as
bubble_lib.B_momentum does in the rest frame, with omega_2 = sqrt((k_vec - p_vec)^2 + a2) now):
    B/i = int d^3k/(2 pi)^3  (w1 + w2) / (2 w1 w2 [(w1 + w2)^2 - p0^2]),   real k_vec,
done as a 2-dimensional integral over (k_par, k_perp) by cylindrical symmetry.  In the rest frame w1 + w2 = 2 Re w1 is real for
all real k and the denominator vanishes on a whole sphere above threshold -- the pinch, where the Lee-Wick/fakeon prescription is
the principal value (bubble_lib, 30 Oct).  At p_vec != 0 the sum w1 + w2 is complex for real k except at the single point
k = p_vec/2, where it equals 2 Re E_{p/2}: no pinch, no prescription needed, the integral is unambiguous -- and that is the point:
its singular point sits at
    s_min(|p|) = 4 (Re E_{p/2})^2 - p_vec^2,   E_{p/2} = sqrt(p_vec^2/4 + a1),         (Nakanishi 2006 eq. 3.6)
which depends on |p_vec| at fixed s = p0^2 - p_vec^2, and lies below the rest-frame threshold 2 M^2 (= (M + M*)^2) for |p| > 0.

Method P (Feynman parameters, the average continuation, bubble_lib.average_continuation) is Lorentz invariant by construction and
was shown on 30 Oct to equal the rest-frame method M with the principal value (results/lw.log).  Rule 3: every number has P and M;
here the P-vs-M comparison in the boosted frame is not a check but the finding.

Predictions (2 Nov 2026, written before the first run; a = +-i, M = 1, s_ref = -1, |p| in {0.3, 0.8}, s in {1.0, 1.8, 2.5, 4.0}):
  P1. For s < s_min(|p|) the boosted integral equals the rest-frame value at the same s to the integration tolerance (both are the
      unique continuation of the Euclidean function along the real axis before any singularity): Lorentz invariance holds below the
      moving threshold.  Scored at s = 1.0 for both |p| (s_min(0.3) = 1.956, s_min(0.8) = 1.707, hand-evaluated from the formula).
  P2. For s > s_min(|p|) the boosted integral differs from the rest-frame value at the same s by more than the tolerance -- the crack;
      at s = 1.8 only |p| = 0.8 has crossed its threshold and differs, |p| = 0.3 has not and agrees; at s = 2.5 and 4.0 both differ.
  P3. The boosted integral is finite and continuous across s_min(|p|) (the denominator vanishes at one point of R^3, an integrable
      1/r singularity), so the crack is a non-analytic point, not a cut with a jump: checked by |B(s_min + h) - B(s_min - h)| -> 0 as
      h -> 0 at |p| = 0.8 (h = 0.1, 0.03, 0.01), not scored beyond "decreasing".
  P4. As |p| -> 0 at fixed s above threshold (s = 2.5; |p| = 0.3, 0.1, 0.03) the boosted integral tends to the rest-frame PRINCIPAL
      VALUE (Lee-Wick) value: Im(w1 + w2) is odd across the pinch sphere for small p_vec, so the two i-epsilon halves average out.
      The frame dependence is continuous at p_vec = 0 and grows with |p|.  Scored as monotone decrease of the difference.
  P5. Method P (average continuation) equals the rest-frame method M (PV) at every s in the grid (the 30 Oct result, re-run here as
      the anchor) and therefore differs from the boosted method M above s_min: P != M located at the undefined operation (the pinch
      prescription vs. the real-momentum integral), which is where rule 1 says the finding lives.
Status PASS if P1, P2 (as a pattern), P4 (monotone) and P5 hold; FAIL otherwise.  The crack is the EXPECTED outcome of this run;
a Lorentz-invariant boosted integral would be the surprise and would say the 1971 claim is wrong as restated.
"""
import sys, time, json, argparse
import mpmath as mp
import bubble_lib as L

ap = argparse.ArgumentParser()
ap.add_argument('--dps', type=int, default=15)
ap.add_argument('--gamma', type=float, default=None, help="Nakanishi's M^2 = 1 + i gamma instead of the a = +-i pair")
ap.add_argument('--pvals', default='0.3,0.8')
ap.add_argument('--svals', default='1.0,1.8,2.5,4.0')
ap.add_argument('--json', default=None)
args = ap.parse_args()
mp.mp.dps = args.dps
t0 = time.time()
I = L.I; PI = L.PI
a1 = mp.mpc(1, args.gamma) if args.gamma is not None else mp.mpc(0, 1)
a2 = mp.conj(a1)
s_ref = mp.mpf(-1)
tol = mp.mpf(10) ** (-(args.dps - 7))
print(f"# bubble_frame: complex-ghost-pair bubble at p_vec != 0, Lee-Wick contour, real k_vec; a1 = {a1}, a2 = {a2}, D = 4, s_ref = {s_ref}, dps = {args.dps}")


def E_half(p):
    return mp.sqrt(mp.mpf(p) ** 2 / 4 + a1)


def s_min(p):
    return 4 * mp.re(E_half(p)) ** 2 - mp.mpf(p) ** 2


def h_boost(kpar, kperp, p0sq, p):
    """integrand of B/i after the k0 residues, boosted frame; p0sq = p0^2 may be negative (Euclidean)"""
    k2 = kpar * kpar + kperp * kperp
    w1 = mp.sqrt(k2 + a1)
    w2 = mp.sqrt((kpar - p) ** 2 + kperp * kperp + a2)
    s = w1 + w2
    return s / (2 * w1 * w2 * (s * s - p0sq))


def B_boost(s, p):
    """UV-subtracted B(s, |p|) = B(p0^2 = s + p^2, p) - B(p0^2 = s_ref + p^2, p), method M in the boosted frame"""
    p = mp.mpf(p); s = mp.mpf(s)
    p0sq = s + p * p; p0sq_ref = s_ref + p * p

    def inner(kpar):
        f = lambda kperp: kperp * (h_boost(kpar, kperp, p0sq, p) - h_boost(kpar, kperp, p0sq_ref, p))
        return mp.quad(f, [0, mp.mpf(1) / 2, 2, mp.inf])

    pts = [-mp.inf, -2, 0, p / 2, p, p + 2, mp.inf]
    val = mp.quad(inner, pts)
    # d^3k/(2 pi)^3 = 2 pi kperp dkperp dkpar / (2 pi)^3
    return I * val / (4 * PI ** 2)


def B_rest(s, prescription='PV'):
    """rest frame, method M of bubble_lib with the same Euclidean subtraction (p0_ref = i, p0_ref^2 = -1)"""
    s = mp.mpf(s)
    p0 = mp.sqrt(s) if s > 0 else I * mp.sqrt(-s)
    return L.B_momentum(4, p0, a1, a2, prescription=prescription, p0_ref=I)


def B_methodP(s):
    """method P: the average continuation of the D = 4 Feynman-parameter integral, same subtraction (I_4(t) - I_4(-1))"""
    s = mp.mpf(s)
    fav = L.average_continuation(4, s, a1, a2)[0]
    fav_ref = L.I_straight(4, s_ref, a1, a2)
    return -I / (16 * PI ** 2) * (fav - fav_ref)


ok = True
out = {'a1': str(a1), 'rows': []}
pvals = [mp.mpf(x) for x in args.pvals.split(',')]
svals = [mp.mpf(x) for x in args.svals.split(',')]
print("# thresholds: rest frame 2 Re(sqrt(a1))^2 = (M + M*)^2 =", mp.nstr(4 * mp.re(mp.sqrt(a1)) ** 2, 8),
      "; s_min(|p|) =", ", ".join(f"{mp.nstr(p, 3)}: {mp.nstr(s_min(p), 8)}" for p in pvals), "  [P2 expects each below the rest-frame value]")

print("# --- P5 anchor: method P (average continuation) vs rest-frame method M (PV) ---")
restM = {}; restP = {}
for s in svals:
    t1 = time.time()
    restM[s] = B_rest(s); restP[s] = B_methodP(s)
    d = abs(restM[s] - restP[s])
    print(f"#   s = {mp.nstr(s, 4)}: M_rest(PV) = {mp.nstr(restM[s], 10)}  P = {mp.nstr(restP[s], 10)}  |diff| = {mp.nstr(d, 3)}  [P5: < {mp.nstr(tol, 2)}] ({time.time()-t1:.0f}s)")
    ok &= d < tol * 100

print("# --- boosted frame, method M with real k_vec: B(s, |p|) against the rest frame at the same s ---")
boost = {}
for p in pvals:
    for s in svals:
        t1 = time.time()
        b = B_boost(s, p); boost[(p, s)] = b
        d = abs(b - restM[s])
        below = s < s_min(p)
        tag = 'P1: agree' if below else 'P2: differ'
        print(f"#   |p| = {mp.nstr(p, 3)}, s = {mp.nstr(s, 4)} ({'below' if below else 'above'} s_min): B = {mp.nstr(b, 10)}  |B - B_rest| = {mp.nstr(d, 3)}  [{tag}] ({time.time()-t1:.0f}s)")
        out['rows'].append(dict(p=float(p), s=float(s), B=[float(mp.re(b)), float(mp.im(b))], diff=float(d), below=bool(below)))
        if below:
            ok &= d < tol * 100
        else:
            ok &= d > tol * 100

print("# --- P3: continuity across s_min at |p| = 0.8 ---")
p = mp.mpf('0.8'); sm = s_min(p); prev = None
for h in (mp.mpf('0.1'), mp.mpf('0.03'), mp.mpf('0.01')):
    t1 = time.time()
    jump = abs(B_boost(sm + h, p) - B_boost(sm - h, p))
    print(f"#   h = {mp.nstr(h, 2)}: |B(s_min + h) - B(s_min - h)| = {mp.nstr(jump, 4)} ({time.time()-t1:.0f}s)")
    if prev is not None and not (jump < prev):
        print("#   P3: not decreasing")
    prev = jump

print("# --- P4: |p| -> 0 at s = 2.5 against the rest-frame principal value and the two i-epsilon prescriptions ---")
s = mp.mpf('2.5')
pv = restM[s] if s in restM else B_rest(s)
fey = B_rest(s, 'feynman'); mnu = B_rest(s, 'minus')
print(f"#   rest frame: PV = {mp.nstr(pv, 10)}; +i eps = {mp.nstr(fey, 10)}; -i eps = {mp.nstr(mnu, 10)}")
prev = None; p4 = True
for p in (mp.mpf('0.3'), mp.mpf('0.1'), mp.mpf('0.03')):
    t1 = time.time()
    b = boost[(p, s)] if (p, s) in boost else B_boost(s, p)
    d = abs(b - pv)
    print(f"#   |p| = {mp.nstr(p, 3)}: B = {mp.nstr(b, 10)}  |B - PV| = {mp.nstr(d, 4)}  |B - (+i eps)| = {mp.nstr(abs(b - fey), 4)} ({time.time()-t1:.0f}s)")
    if prev is not None:
        p4 &= d < prev
    prev = d
print(f"#   P4 (difference to the PV decreases as |p| -> 0): {p4}")
ok &= p4

if args.json:
    json.dump(out, open(args.json, 'w'), indent=1)
print(f"# END bubble_frame elapsed={time.time()-t0:.1f}s status={'PASS' if ok else 'FAIL'}")
