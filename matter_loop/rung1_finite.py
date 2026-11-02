#!/usr/bin/env python3
"""
Rung 1, the last owed item of Utiyama-DeWitt 1962: the FINITE p^4 log p^2 part of the matter-loop polarisation,
printed and fitted, and the 1/p^4 resummation of the abstract's "semiclassical graviton Green's function" shown.

Massless scalar (m = 0; the p^4 part at large p is mass-independent at leading order and m = 0 makes the x-integrals
elementary), one loop, second order in h, Euclidean.  The p^4 part of Q = Gamma_2/(2 Vol) is bubble-only because the
tadpole carries no p; it is nevertheless the full sum at this order, not a single-graph number.

Three regulators on the same graph:  (a) hard cutoff |l| < Lambda on the Feynman-shifted momentum;  (b1) Schwinger
proper time t >= 1/Lambda^2;  (c) dimensional regularization, n = 4 + eps_d, MS-bar, with the n-DIMENSIONAL angular
averages (the earlier dim-reg driver took them in d = 4, which is exact for the pole and wrong for the finite part).
For each: the p^4 part is written  c_L * log(Lambda^2) [or epsinv]  +  c_p * log(p^2)  +  K,  each a quadratic form
in h fitted to {Ric^2, R^2} and to the Barnes-Rivers forms {p^4 h.P2.h, p^4 h.P0s.h}.

Predictions (1 Nov 2026, written before the run):
  P1. c_p = -c_L in all three regulators (dictionary (1/16pi^2) log Lambda^2 <-> -1/eps for (c)): the p^4 term is
      c_L log(Lambda^2/p^2) + K.  This is a check of the bookkeeping, not physics; it must hold.  (Likewise the dim-reg
      "coefficient of log mu^2 equals c_L" line cannot fail: the same residue feeds both.  Dim reg's K below is the
      MS-bar constant with mu^2 identified with Lambda^2.)
  P2. K differs between the three regulators, and each difference is a combination of Ric^2 and R^2 only
      (local; absorbable in the quadratic-curvature counterterms) -- the regulator-dependent constant is shown and
      is unobservable.  The (a)-(b1) difference is expected to involve Euler's gamma; (c) is MS-bar.
      [Outcome, written after the runs: P2 as written FAILED for the hard cutoff (either measure) and for proper time
      with the plain measure, whose K are not transverse; it holds between dim reg and proper time with DeWitt's
      measure.  The END status label LOG-OK_FINITE-NONTRANSVERSE[...] and the rule that P2 pairs containing a
      non-transverse K are reported but not scored were introduced after the first cutoff-only run, which had ended
      status=OK with the same non-transverse K; they describe the result, not a prediction.]
  P3. Ric^2 = (1/4) p^4 h.P2.h + p^4 h.P0s.h and R^2 = 3 p^4 h.P0s.h (checked at the top of the run), so
      Ric^2 - R^2/3 is pure spin-2; with the log coefficients -(1/120) Ric^2 - [(1/4)(xi-1/6)^2 - 1/360] R^2
      the spin-0 coefficient of the p^4 log is -(1/120) - 3[(1/4)(xi-1/6)^2 - 1/360] = -(3/4)(xi-1/6)^2, which
      VANISHES AT xi = 1/6 (conformal coupling: no spin-0 running), and the spin-2 coefficient is -1/480 for all xi.
  P4. Tree level, S_E = -(M_P^2/2) int sqrt(g) R, quadratic part from (sqrt g R)^(2) derived in code: the spin-2
      coefficient of p^2 is positive and the spin-0 coefficient is negative (the Euclidean conformal-factor sign).
  P5. Resummed spin-2 inverse propagator, Euclidean:  G2^-1 = a M_P^2 p^2 + b p^4 log(p^2/Lambda^2) + k p^4, with
      a > 0 (P4) and b = +(1/480)(2/16pi^2) > 0 (P3): the induced spin-2 p^4 coefficient GROWS with p^2 -- the
      asymptotic-freedom sign of the Weyl-squared coefficient (positive c) -- and is therefore NEGATIVE for
      p^2 < Lambda^2 e^{-k/b}.  Consequences, to be printed and not interpreted further here:
      (i) for M_P^2 << p^2 << Lambda^2,  G2 ~ 1/(b p^4 log(p^2/Lambda^2)) -- the abstract's 1/p^4 with a running
          coefficient; (ii) with NO bare C^2 term at Lambda, G2^-1 has a Euclidean zero at p^2 = a M_P^2/(b log(Lambda^2/p^2) - k)
          for Lambda large enough (a Euclidean pole is a Lorentzian tachyon, not Stelle's ghost); Stelle's ghost
          (massive spin-2 with negative residue, mass^2 = a M_P^2/B) is what one gets instead when a bare
          Euclidean C^2 coefficient B > 0 at the scale p dominates the negative log.  Which of the two the
          1962 abstract meant by "semiclassical Green's function ~ 1/p^4" cannot be decided from the abstract.
      RECORD OF THE CORRECTION (rule 11).  The first draft of this header, written before any run of this script, read:
        "P5. Resummed spin-2 inverse propagator G2^-1 = a M_P^2 p^2 + b p^4 log(Lambda^2/p^2) + k p^4 with a, b > 0
         (Euclidean, p < Lambda): no Euclidean pole; for p^2 >> M_P^2 (and p << Lambda), G2 ~ 1/(b p^4 log(Lambda^2/p^2)),
         the abstract's 1/p^4 with a running coefficient; continued to Lorentzian timelike momentum the real part of G2^-1
         vanishes at |p^2| = a M_P^2 / (b log(Lambda^2/|p^2|) + k): the massive spin-2 pole of Stelle 1977 with a
         cutoff-dependent, logarithmically running mass.  Printed, not interpreted further here."
      It was replaced by the paragraph above in the same session, before the script was first executed (the first
      execution failed on the tree-level fit, an unrelated error), because it contradicts P3's sign: with c_L = -1/480 on
      spin-2, b log(Lambda^2/p^2) with b > 0 is the wrong sign convention and the induced term is negative below Lambda.
      The sandbox has no run of this script under the first wording.  What the log tests of P5 is only a > 0 and b > 0,
      which are P4 and P3 restated; the location of the zero is printed as a formula and not tested.
"""
import time, argparse, itertools
import sympy as sp
import ud1962_lib as U

ap = argparse.ArgumentParser()
ap.add_argument('--measure', choices=['plain', 'covariant'], default='plain')
ap.add_argument('--regulators', default='cutoff,propertime,dimreg')
ap.add_argument('--log', default=None)
args = ap.parse_args()
t0 = time.time()
out = open(args.log, 'w') if args.log else None
def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    if out: out.write(s + '\n'); out.flush()

xi = U.xi
pi2 = 16*sp.pi**2
name = f'rung1 finite {args.measure}'
say(f'# {name}: finite p^4 log p^2 part, massless scalar, regulators {args.regulators}')
vertex = U.V1cov if args.measure == 'covariant' else U.V1
ok = True

say('# --- spin projectors ---')
f2, f0 = U.projector_forms()
sR, rR = U.fit(U.Q_Ric, [f2, f0], ['P2', 'P0s'])
sS, rS = U.fit(U.Q_R2, [f2, f0], ['P2', 'P0s'])
say(f'#   Ric^2 = {sR} (resid {rR});  R^2 = {sS} (resid {rS})   [P3 expects 1/4, 1 and 0, 3]')
ok &= sR == {'P2': sp.Rational(1, 4), 'P0s': 1} and sS == {'P2': 0, 'P0s': 3}

say('# --- tree level: -(M_P^2/2) int sqrt(g) R at second order, derived from the metric ---')
MP = sp.Symbol('M_P', positive=True)
ME = U.metric_expansion()
sqrtg2, A2, sqrtgR2 = U.second_order_densities(ME)
Qtree = sp.expand(-MP**2/2*sqrtgR2)
# fit p^2 Q_tree onto the p^4 forms (p^2 h.P2.h alone is not a polynomial; the Einstein combination P2 - 2 P0s is)
sT, rT = U.fit(sp.expand(U.p2*Qtree), [f2, f0], ['P2', 'P0s'])
say(f'#   Q_tree = {sT}, residual {rT}   [P4: P2 coefficient > 0, P0s coefficient < 0]')
ok &= sT is not None and rT == 0 and (sT['P2']/MP**2) > 0 and (sT['P0s']/MP**2) < 0

results = {}
for reg in args.regulators.split(','):
    say(f'# --- regulator: {reg} ---')
    t1 = time.time()
    tab = U.bubble_finite(reg, covariant=(args.measure == 'covariant'))
    Q = sp.Integer(0)
    for (i, j, r, s) in itertools.product(range(4), repeat=4):
        Q += sp.Rational(1, 2)*U.H[i, j]*U.sym_get(tab, i, j, r, s)*U.H[r, s]
    Q = sp.expand(Q)
    if reg == 'dimreg':
        # coefficient of epsinv is already "per 1/eps_tHV" (dictionary (1/16pi^2) log Lambda^2 <-> -1/eps); the finite
        # part carries its own 1/16pi^2 and is put in the same x16pi^2 units as the other regulators
        Qp = sp.expand(Q.coeff(U.epsinv, 1)); Qf = sp.expand(Q - Qp*U.epsinv)
        Q = sp.expand(-Qp*2*U.lgL + pi2*Qf)       # pole written as its log-Lambda^2 equivalent: epsinv -> -log Lambda^2 (x16pi^2 units)
    else:
        Q = sp.expand(pi2*Q)
    Q = sp.expand(Q.subs(U.Ssym, U.p2))          # S stood for p^2 inside Delta; log p^2 stays as the symbol LS
    Q4 = U.deg_part(Q, 4)
    low = sp.expand(Q - Q4)
    low_noL = sp.expand(sum(t_ for t_ in sp.Add.make_args(low) if not t_.has(U.Lam)))
    low_log = sp.expand(sum(t_ for t_ in sp.Add.make_args(low) if t_.has(U.lgL) or t_.has(U.LS) or t_.has(U.epsinv)))
    say(f'#   m = 0: below p^4 only power divergences remain: Lambda-free part zero: {low_noL == 0}; log part zero: {low_log == 0}')
    ok &= low_noL == 0 and low_log == 0
    cL = sp.expand(Q4.coeff(U.lgL, 1)/2)
    rest = sp.expand(Q4 - 2*cL*U.lgL)
    if reg == 'dimreg':
        cmu = sp.expand(rest.coeff(sp.log(U.mu), 1)/2)
        rest = sp.expand(rest - 2*cmu*sp.log(U.mu))
        say(f'#   dim reg: pole reported through the dictionary 1/eps <-> -(1/16pi^2) log Lambda^2; coefficient of log mu^2 equals c_L: {sp.expand(cmu - cL) == 0}')
        ok &= sp.expand(cmu - cL) == 0
    cp = sp.expand(rest.coeff(U.LS, 1))
    K = sp.expand(rest - cp*U.LS)
    assert not K.has(U.LS) and not K.has(U.lgL) and not K.has(U.epsinv), 'leftover logs in the constant'
    fL, rL = U.fit(cL, [U.Q_Ric, U.Q_R2], ['Ric2', 'R2'])
    fp, rp = U.fit(cp, [U.Q_Ric, U.Q_R2], ['Ric2', 'R2'])
    fK, rK = U.fit(K, [U.Q_Ric, U.Q_R2], ['Ric2', 'R2'])
    say(f'#   c_L (log Lambda^2): {fL}, resid {rL}   [expected -1/120, -(60 xi^2 - 20 xi + 1)/240]')
    say(f'#   c_p (log p^2):      {fp}, resid {rp}   [P1 expects c_p = -c_L]')
    say(f'#   K  (constant):      {fK}, resid {rK if rK == 0 else "NONZERO (" + str(len(sp.Add.make_args(rK))) + " terms)"}')
    p1 = sp.expand(cp + cL) == 0
    say(f'#   P1 c_p = -c_L: {p1}')
    ok &= p1 and rL == 0 and rp == 0
    ok &= fL is not None and sp.simplify(fL['Ric2'] + sp.Rational(1, 120)) == 0 and sp.simplify(fL['R2'] + (60*xi**2 - 20*xi + 1)/240) == 0
    # spin decomposition of the log coefficient
    sL, _ = U.fit(cL, [f2, f0], ['P2', 'P0s'])
    say(f'#   c_L on projectors: P2 = {sL["P2"]}, P0s = {sp.factor(sL["P0s"])}   [P3: -1/480 and -(3/4)(xi-1/6)^2]')
    ok &= sp.simplify(sL['P2'] + sp.Rational(1, 480)) == 0 and sp.simplify(sL['P0s'] + sp.Rational(3, 4)*(xi - sp.Rational(1, 6))**2) == 0
    # linearized Ward identity, bilinear form: dK/dh_ij . (p_i e_j + p_j e_i) must vanish identically in h and e
    # (the weaker test K(pure gauge) = 0 is NOT equivalent -- the 1 Nov review exhibited (p h p)^2 - p^4 h^2 passing it)
    E_ = sp.symbols('e1:5', real=True)
    ward = sp.expand(sum(sp.diff(K, U.H[i, j])*(U.P[i]*E_[j] + U.P[j]*E_[i]) for i in range(4) for j in range(i, 4)))
    say(f'#   K transverse (first-order variation under h -> h + p e + e p vanishes identically): {ward == 0}')
    sK, _ = U.fit(K, [f2, f0], ['P2', 'P0s'])
    if sK is None:
        # project out the transverse part anyway for the P5 printout: fit K on {P2, P0s} plus the two non-transverse p^4 scalars
        b5 = sum(U.P[i]*U.P[j]*U.H[i, j] for i in range(4) for j in range(4))**2
        b6 = U.p2*sum(U.P[i]*U.P[j]*U.H[i, k]*U.H[j, k] for i in range(4) for j in range(4) for k in range(4))
        b7 = U.p2**2*U.h_tr**2
        sK6, rK6 = U.fit(K, [f2, f0, b5, b6, b7], ['P2', 'P0s', 'pph_sq', 'p2_phhp', 'p4_h2'])
        say(f'#   K on the complete basis {{P2, P0s, (p h p)^2, p^2 (p h h p), p^4 h^2}}: {sK6}, resid {rK6}')
        sK = {'P2': sK6['P2'], 'P0s': sK6['P0s']} if sK6 else {'P2': sp.nan, 'P0s': sp.nan}
    results[reg] = dict(cL=cL, cp=cp, K=K, fK=fK, sL=sL, sK=sK, ward=(ward == 0))
    say(f'#   ({time.time()-t1:.0f}s)')

say('# --- P2: regulator differences in the constant K ---')
regs = list(results)
for a_, b_ in itertools.combinations(regs, 2):
    d = sp.expand(results[a_]['K'] - results[b_]['K'])
    fd, rd = U.fit(d, [U.Q_Ric, U.Q_R2], ['Ric2', 'R2'])
    both_t = results[a_]['ward'] and results[b_]['ward']
    say(f'#   K[{a_}] - K[{b_}] = {fd}, resid {rd if rd == 0 else "NONZERO (" + str(len(sp.Add.make_args(rd))) + " terms)"}   '
        f'[P2: local and covariant means resid 0{"" if both_t else "; not scored: a non-transverse K is in the pair"}]')
    if both_t:
        ok &= fd is not None and rd == 0

say('# --- P5: resummed inverse propagator per spin sector (units: coefficient of h.P.h in 2Q; 1/16pi^2 restored) ---')
LamL = sp.Symbol('Lambda', positive=True); S = sp.Symbol('p2', positive=True)
for reg in regs:
    sL, sK = results[reg]['sL'], results[reg]['sK']
    for sec in ('P2', 'P0s'):
        tree = 2*sT[sec]*S
        loop = 2*(sL[sec]*sp.log(LamL**2/S) + sK[sec])*S**2/pi2
        say(f'#   [{reg}] G^-1_{sec}(p^2) = {sp.simplify(tree)} + {sp.factor(sp.simplify(loop))}')
    a_ = 2*sT['P2']; bb = -2*sL['P2']/pi2; kk = 2*sK['P2']/pi2   # G2^-1 = a p^2 + bb p^4 log(p^2/Lambda^2) + kk p^4
    say(f'#   [{reg}] spin-2: a = {sp.simplify(a_)}, b = {sp.simplify(bb)} (coefficient of p^4 log(p^2/Lambda^2)), k = {sp.simplify(kk)}'
        f'  [P5: a > 0, b > 0: {bool(a_/MP**2 > 0 and bb > 0)}]')
    say(f'#   [{reg}] spin-2, M_P^2 << p^2 << Lambda^2:  G_2 ~ 1/({sp.simplify(bb)} p^4 log(p^2/Lambda^2))')
    say(f'#   [{reg}] spin-2, no bare C^2 term: Euclidean zero of G_2^-1 at  p^2 = ({sp.simplify(a_)}) / ({sp.simplify(bb)} log(Lambda^2/p^2) - ({sp.simplify(kk)}))')
    say(f'#   [{reg}] spin-0 p^4 log coefficient {sp.factor(sL["P0s"])}: zero at xi = 1/6')

all_ward = all(v['ward'] for v in results.values())
status = 'OK' if (ok and all_ward) else ('LOG-OK_FINITE-NONTRANSVERSE[' + ','.join(k for k, v in results.items() if not v['ward']) + ']' if ok else 'MISMATCH')
say(f'# END {name} elapsed={time.time()-t0:.1f}s status={status}')
if out: out.close()
