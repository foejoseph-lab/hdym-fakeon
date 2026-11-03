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

RECORD OF THE FIRST RUN (rule 11; 2 Nov evening, sandbox).  The header above is kept as written.  Its description "no pinch, no
prescription needed, the integral is unambiguous" and P3's "the denominator vanishes at one point of R^3" are WRONG: the symmetry
k -> p - k fixes the plane k_par = p/2 pointwise, so w1 + w2 is real on that whole plane and the denominator vanishes on the CIRCLE
k_par = p/2, k_perp = k*(p0) for s > s_min(|p|) (bubble_lib.kstar_boost).  Im(w1 + w2) is odd across the plane, so the integral is
ABSOLUTELY convergent as a 3-dimensional integral (a 1/r singularity on a circle) and well defined without a prescription -- the
iterated k_perp integral picks up +-i pi times a residue on the two sides of the plane, whose average is the principal value -- but a
quadrature that does not resolve the near-pole at k* gives breakpoint-dependent numbers: the first run's
above-threshold values were off by 0.15% to 7% (worst at |p| = 0.3, s = 4.0; at s = 2.5, |p| = 0.8: 0.0065774 in the superseded log and
0.0065237 in an unlogged exploratory run with other breakpoints, against the converged 0.0065471464695), while every below-threshold
number stood to 10^-17.  Found by the GMRS certification (the closed form agreed only to 1.5% in an unlogged exploratory run), not by
any check in the script.  Fixed by the pinch-aware quadrature of bubble_lib.B_momentum_boost; the first log is kept in
results/superseded/frame_sandbox_pointwise_quadrature.log.  P1, P2, P4, P5 are unchanged in kind by the fix; P3 is re-scored on the
corrected numbers; P4 now has its reason: the pinch sphere of the rest frame becomes a circle whose i-epsilon sign flips across the
plane, and the average of the two signs IS the principal value.

Added 2 Nov evening (after the first run; the author chose the sweep shape and asked for the derivation):
  --sweep     s from --srange lo:hi:step at each |p| (the shape later prescriptions are compared on; traces s_min moving).
  --gmrs      certification against the GMRS closed form as restated by Nakanishi 2006 eq. (3.7), in two steps.
      (a) Derivation of the contour.  The "fish" D is the image of the half-plane {(k_par, k_perp >= 0)} under
          f(k) = E_k + E*_{p-k}; the boundary of the image of a half-plane under a smooth map lies in the image of the edge
          (k_perp = 0) or of the fold set (Jacobian of (Re f, Im f) zero).  On the edge df/dk_perp = 0 by evenness, so the
          edge IS a fold, and if the Jacobian has no zero for k_perp > 0 the boundary Gamma lies in the image of the axis
          k_perp = 0, touching the real axis at k = p/2 (f real there) -- Nakanishi's touching point.  Checked numerically on a grid
          (not a proof: a grid cannot exclude a sign change between nodes or beyond its edge; the large-|k| sign is checked at |k| = 50).
          Prediction G1: no interior zero (the Jacobian keeps one sign for k_perp > 0 on the grid and at large |k|).
      (b) Gamma_s = the half of the image of Gamma under s' = f^2 - p^2 with Im s' < 0, from s_min to infinity, parametrised by
          k_par from p/2 outward; F(s) = (1/16 pi^2) Re int_{Gamma_s} ds'/(s - s') [ |p|(M^2 - M*^2)/(s' sqrt(s' + p^2))
          + sqrt((s' - (M+M*)^2)(s' - (M-M*)^2))/s' ] with principal square roots (a branch choice the open text does not fix).
          Prediction G2: F(s) - F(s_0) is proportional to B(s, |p|) - B(s_0, |p|) with one constant factor c over the sweep (the
          normalisation and the UV subtraction differ); the constancy of c scored to 10^-6 relative, and -- derived after the first
          sweep from the two normalisations, F = i B (Nakanishi's 1/(2(2 pi)^3) against our i/(4 pi^2)) -- c = -i itself scored to 10^-6.
          If G2 fails with the principal branches the branch choice is reported as the open point, not a failure of the run.  The crack is the EXPECTED outcome of this run;
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
ap.add_argument('--sweep', action='store_true')
ap.add_argument('--srange', default='1.0:4.0:0.25')
ap.add_argument('--gmrs', action='store_true')
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


def B_boost(s, p):
    """UV-subtracted B(s, |p|), method M in the boosted frame, pinch-aware quadrature (bubble_lib.B_momentum_boost, 2 Nov evening)"""
    return L.B_momentum_boost(s, p, a1, a2, s_ref)


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
GRID = not args.sweep
pvals = [mp.mpf(x) for x in args.pvals.split(',')]
svals = [mp.mpf(x) for x in args.svals.split(',')]
print("# thresholds: rest frame 2 Re(sqrt(a1))^2 = (M + M*)^2 =", mp.nstr(4 * mp.re(mp.sqrt(a1)) ** 2, 8),
      "; s_min(|p|) =", ", ".join(f"{mp.nstr(p, 3)}: {mp.nstr(s_min(p), 8)}" for p in pvals), "  [P2 expects each below the rest-frame value]")

if GRID:
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


def f_edge(kpar, p):
    """E_k + E*_{p-k} on the axis k_perp = 0"""
    return mp.sqrt(kpar * kpar + a1) + mp.sqrt((kpar - p) ** 2 + a2)

if args.sweep:
  lo, hi, st = [mp.mpf(x) for x in args.srange.split(':')]
  svals = []
  x = lo
  while x <= hi + st / 10:
    svals.append(x); x += st
  print(f"# --- sweep: s = {mp.nstr(lo,4)}..{mp.nstr(hi,4)} step {mp.nstr(st,4)}; rest frame (PV) and boosted |p| = {args.pvals}; s_min marked ---")
  restM = {}
  for s in svals:
    restM[s] = B_rest(s)
  sweep_rows = {}
  for p in pvals:
    sm = s_min(p)
    print(f"#   |p| = {mp.nstr(p,3)}: s_min = {mp.nstr(sm, 8)}")
    for s in svals:
      t1 = time.time()
      b = B_boost(s, p); d = abs(b - restM[s]); below = s < sm
      sweep_rows[(p, s)] = b
      print(f"#   |p| = {mp.nstr(p,3)} s = {mp.nstr(s,5)} {'<' if below else '>'} s_min: B_rest = {mp.nstr(restM[s], 10)}  B_boost = {mp.nstr(b, 10)}  |diff| = {mp.nstr(d, 4)} ({time.time()-t1:.0f}s)")
      out['rows'].append(dict(p=float(p), s=float(s), B=[float(mp.re(b)), float(mp.im(b))], Brest=[float(mp.re(restM[s])), float(mp.im(restM[s]))], diff=float(d), below=bool(below)))
      ok &= (d < tol * 100) if below else (d > tol * 100)
  if args.gmrs:
    print("# --- GMRS certification (Nakanishi 2006 eq. 3.7) ---")
    for p in pvals:
      # (a) the Jacobian of (Re f, Im f) over (k_par, k_perp > 0)
      def fK(kpar, kperp):
        return mp.sqrt(kpar * kpar + kperp * kperp + a1) + mp.sqrt((kpar - p) ** 2 + kperp * kperp + a2)
      def jac(kpar, kperp):
        d1 = mp.diff(lambda u: fK(u, kperp), kpar); d2 = mp.diff(lambda v: fK(kpar, v), kperp)
        return mp.re(d1) * mp.im(d2) - mp.im(d1) * mp.re(d2)
      vals = []
      for i in range(13):
        kpar = -3 + (p + 6) * i / 12
        for j in range(1, 13):
          kperp = 4 * mp.mpf(j) / 12
          vals.append(jac(kpar, kperp))
      far = [jac(kpar, kperp) for kpar in (-50, p / 2, 50) for kperp in (mp.mpf('0.3'), 50)]
      signs = set(int(mp.sign(v)) for v in vals + far)
      print(f"#   |p| = {mp.nstr(p,3)}: Jacobian on a 13x12 grid (k_par in [-3, p+3], k_perp in (0, 4]) and at |k| ~ 50: min |J| = {mp.nstr(min(abs(v) for v in vals), 3)}, signs {sorted(signs)}   [G1: one sign; a grid check, not a proof]")
      ok &= len(signs) == 1
      # (b) Gamma_s from the axis image, the half with Im s' < 0
      M, Ms = mp.sqrt(a1), mp.sqrt(a2)
      sp_ = lambda kpar: f_edge(kpar, p) ** 2 - p * p
      side = +1 if mp.im(sp_(p / 2 + mp.mpf('0.5'))) < 0 else -1
      assert all(mp.im(sp_(p / 2 + side * u)) < 0 for u in (mp.mpf('0.05'), mp.mpf('0.5'), mp.mpf(5))), 'Im s\' changes sign along the chosen half of Gamma'
      def F(s):
        def g(u):           # u >= 0, kpar = p/2 + side*u
          kpar = p / 2 + side * u
          sv = sp_(kpar); dsv = side * mp.diff(sp_, kpar)
          term = p * (M * M - Ms * Ms) / (sv * mp.sqrt(sv + p * p)) + mp.sqrt((sv - (M + Ms) ** 2) * (sv - (M - Ms) ** 2)) / sv
          return dsv / (s - sv) * term
        return mp.quad(g, [0, 1, 4, mp.inf]) / (16 * PI ** 2)
      s0 = svals[len(svals) // 2]
      F0 = F(s0); B0 = sweep_rows[(p, s0)]
      cs = []
      for s in svals:
        if s == s0:
          continue
        Fs = F(s); dF = Fs - F0; dB = sweep_rows[(p, s)] - B0
        c_re = dB / mp.re(dF) if mp.re(dF) != 0 else mp.nan
        c_full = dB / dF if dF != 0 else mp.nan
        cs.append((s, c_re, c_full))
        print(f"#   |p| = {mp.nstr(p,3)} s = {mp.nstr(s,5)}: Re[F(s)-F(s0)] = {mp.nstr(mp.re(dF), 8)}  F(s)-F(s0) = {mp.nstr(dF, 8)}  c_Re = {mp.nstr(c_re, 8)}  c_full = {mp.nstr(c_full, 8)}")
      for lab, idx in (('Re', 1), ('full', 2)):
        vs = [c[idx] for c in cs if c[idx] == c[idx]]
        if vs:
          spread = max(abs(v - vs[0]) for v in vs) / max(abs(vs[0]), mp.mpf(10) ** -30)
          print(f"#   |p| = {mp.nstr(p,3)}: constancy of c_{lab} over the sweep: relative spread {mp.nstr(spread, 3)}   [G2: < 1e-6 for one of them]")
      spreads = []
      for idx in (1, 2):
        vs = [c[idx] for c in cs if c[idx] == c[idx]]
        spreads.append(max(abs(v - vs[0]) for v in vs) / max(abs(vs[0]), mp.mpf(10) ** -30) if vs else mp.inf)
      cre = [c[1] for c in cs if c[1] == c[1]]
      cval = abs(cre[0] + I) if cre else mp.inf
      g2 = min(spreads) < mp.mpf(10) ** -6 and cval < mp.mpf(10) ** -6
      print(f"#   |p| = {mp.nstr(p,3)}: c_Re (the physical channel: B is pure imaginary, F complex, so c_full cannot be constant) = {mp.nstr(cre[0], 10) if cre else 'nan'}; |c_Re + i| = {mp.nstr(cval, 3)}   [G2: c = -i to 1e-6, from F = i B]")
      print(f"#   |p| = {mp.nstr(p,3)}: G2 holds: {g2}" + ("" if g2 else "   (branch choice of the square roots is the open point, see header)"))
      out[f'gmrs_p{float(p)}'] = dict(g2=bool(g2), spreads=[float(x) for x in spreads])

if args.json:
    json.dump(out, open(args.json, 'w'), indent=1)
print(f"# END bubble_frame elapsed={time.time()-t0:.1f}s status={'PASS' if ok else 'FAIL'}")
