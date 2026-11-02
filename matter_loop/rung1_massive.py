#!/usr/bin/env python3
"""
Rung 1, back pass to Sakharov 1967 (induced gravity): the MASSIVE scalar loop at the orders below p^4.

The same graphs as rung1_run.py / rung1_finite.py (bubble + tadpole, one real scalar, second order in h, Euclidean), with
m != 0 and the external momentum expanded about p = 0 (Delta = m^2 + x(1-x)p^2 Taylor-expanded; exact at each order in
p^2 since the integrals are elementary).  Reported per regulator and measure, at the orders p^0 (cosmological) and p^2
(Newton), each written as

    Q(p^d) = c_4 Lambda^4 + c_2 Lambda^2 + c_L log Lambda^2 [or the dim-reg pole, through the dictionary] + c_m log m^2 + K,

each coefficient a quadratic form in h fitted to the covariant density of that order (sqrt(g)^(2) at p^0, (sqrt g R)^(2) at
p^2), with the bilinear linearized Ward test on every piece at p^2, and the Lambda^2 term printed BESIDE the m^2 log term.
Regulators: (a) Feynman-shifted hard cutoff (1962); (b1) Schwinger proper time (1951/65); (c) dimensional regularization
with n-dimensional angular averages and dot products, MS-bar (1972-74).  Measures: plain (phi) and covariant (chi = g^{1/4} phi).

Why this run (handoff of 1 Nov, agreed with the author): Sakharov 1967 read the m^2 R and Lambda^2 R terms of this very loop
as the ORIGIN of the Einstein term (1/G induced by matter).  The line was taken up (Zel'dovich; Adler RMP 1982 [m]) and
died of a diagnosed cause (Visser 2002 [m]: the induced G of massless fields is regulator-dependent).  Treatment: no
literature hunt; one run; the "-1/12 discipline" applied -- a regulated finite number is meaningful only with the subtracted
piece shown beside it and the unobservability of the subtraction argued -- and a verdict paragraph plus one ledger row.

Predictions (2 Nov 2026, written before the first run of this script.  P0-P1 are bookkeeping and must hold; P2-P4 are
derived from the heat-kernel form of the same quantity, W = -1/2 int_{1/Lambda^2} dt/t Tr e^{-tD} with Tr e^{-tD} =
(16 pi^2 t^2)^-1 int sqrt g e^{-t m^2} (a_0 + t a_1 + t^2 a_2 + ...), a_1 = (1/6 - xi) R, and the MS-bar Gamma-function
constants; they test whether the GRAPHS, with their vertices, measures and angular averages, reproduce that form at the
massive orders, which rung 1 has so far checked only for the divergent parts):
  P0. The Lambda^4, Lambda^2 and log Lambda^2 parts at p^0 and p^2, and the log Lambda^2 part at p^4, reproduce the committed
      rung1_<regulator>_<measure>.log lines of 31 Oct: log p^0 = -1/4 m^4 sqrt(g)^(2); log p^2 = (1/6 - xi)/2 m^2 (sqrt g R)^(2);
      log p^4 = -(1/120) Ric^2 - ((60 xi^2 - 20 xi + 1)/240) R^2 (units 1/16pi^2), in every regulator and measure; Lambda^4
      and Lambda^2 covariant ONLY for proper time with chi (cutoff: neither; proper time with phi: neither; the 31 Oct logs).
      In dim reg the pole at p^0 and p^2 equals the cutoff's log Lambda^2 coefficient (same residue) -- the massive dim-reg
      run is new (rung1_dimreg.py of 31 Oct is massless).
  P1. c_m = -c_L at p^0 and at p^2 in all cases: the mass log and the cutoff log come as log(Lambda^2/m^2) [dim reg:
      log(mu^2/m^2)].  Bookkeeping; must hold.  (At p^4 the mass enters only as the low-energy expansion and the massless
      constant of rung1_finite.py is not re-derived here.)
  P2. The constants.  K_0 (p^0, on m^4 sqrt(g)^(2)) and K_2 (p^2, on m^2 (sqrt g R)^(2)) are covariant/transverse for
      dim reg (either measure) and for proper time with chi; NOT for the hard cutoff (either measure) and NOT for proper
      time with phi -- the pattern of the finite p^4 constant (results note 8.2) repeated at the massive orders.
  P3. Where transverse, proper time (chi) differs from dim reg by Euler's constant times the log coefficient and nothing
      else, at both orders:  K_d[propertime, chi] - K_d[dimreg] = -gamma_E * c_L^(d),  d = 0, 2  (as at p^4).
  P4. The dim-reg (MS-bar, mu^2 == Lambda^2 in the dictionary) constants are rational multiples of the log coefficients:
      K_0[dimreg] = (3/2) c_L^(0) = -3/8  (on m^4 sqrt(g)^(2)),   K_2[dimreg] = 1 * c_L^(2) = (1/6 - xi)/2  (on m^2 (sqrt g R)^(2)).
      [From Gamma(-2 + e) = (1/2)(1/e - gamma + 3/2) + ... and Gamma(-1 + e) = -(1/e - gamma + 1) + ... in MS-bar; written
      with the convention that the dictionary is 1/eps_tHV <-> -(1/16pi^2) log Lambda^2 as in rung1_finite.py.]
  P5. The Sakharov reading, printed and not interpreted further (a reading of a_1, not a test): with S_E = -(M_P^2/2) int sqrt g R,
      the induced Planck mass is M_P^2(ind) = -(2/16pi^2) [c_2 Lambda^2 + c_L m^2 log(Lambda^2/m^2) + K_2] per (sqrt g R)^(2),
      i.e. (1/16pi^2)(1/6 - xi)[Lambda^2 - m^2 log(Lambda^2/m^2) - m^2 (+ gamma_E m^2 for proper time)] for the covariant cases
      and (1/16pi^2)(1/6 - xi) m^2 [log(m^2/mu^2) - 1] in dim reg.  The Lambda^2 piece FAILS the first slot of the -1/12
      discipline (absent in dim reg; non-covariant under the hard cutoff with either placement -- 31 Oct and 1 Nov logs).
      The m^2 log piece PASSES it: its coefficient (1/6 - xi)/2 is the same in every regulator and measure, its sign is
      that of (1/6 - xi), zero at conformal coupling; it fixes only the running d M_P^2(ind)/d log mu^2 = -(1/16pi^2)(1/6 - xi) m^2,
      never the value, which is K_2-dependent (scheme).  Massless fields give nothing in dim reg (Visser's point, in the code).

RECORD OF THE FIRST RUN (rule 11; 2 Nov 2026, sandbox).  With the header as above, the first full plain-measure run gave
K_0[dimreg] = -1/8 (P4 FAILED at p^0; everything else as predicted, P4 at p^2 included), while the covariant-measure run, started
at the same time, gave -3/8 and passed P2-P4 entirely.  Cause, found by redoing the p^0 order exactly in n with the same
n-dimensional averages (the two scripts agreed on the bubble and differed on the tadpole): the tadpole symbol is built from the
4-dimensional x-space expansion and carries only k_1..k_4, so the flat trace term sqrt(g) delta^{mu nu} k_mu k_nu lacked the
k_perp^2 of the n - 4 extra loop-momentum components; k_perp^2 -> (n-4)/n k^2 has no pole and feeds the finite part through the
pole (c_1 A).  With DeWitt's measure the coefficient of k_perp^2 is (g^{-1/4} sqrt g g^{-1/4})^(2) = 0, which is why that run was
right.  Fixed in ud1962_lib.perp_symbol / tadpole_dimreg_full; the first plain log is kept as
results/superseded/rung1_massive_plain_sandbox_4d_tadpole_symbol.log.  Nothing in the 31 Oct and 1 Nov logs is touched: the
massless p^4 part is bubble-only and rung1_dimreg.py reports the pole only, which the k_perp^2 term does not change.  The
prediction P4 was not altered; the code was.  The primary, prediction-independent reason the fix is right: in dim reg the two
measures must agree exactly (the g^{-1/4} factors contribute int d^n k * 1 = 0), the first run broke that (-1/8 with phi against
-3/8 with chi), and after the fix the dim-reg blocks of the two logs are identical line for line.  Also for the record: a
cutoff-only smoke preceded the first full run, after which the status scoring was adjusted so that P3/P4 are "not scored" rather
than MISMATCH when dim reg is not in the run (bookkeeping, no prediction touched).
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

xi, m = U.xi, U.m
pi2 = 16*sp.pi**2
cov = args.measure == 'covariant'
name = f'rung1 massive {args.measure}'
say(f'# {name}: massive scalar loop at p^0 and p^2 (Sakharov 1967 back pass), regulators {args.regulators}')
ok = True

say('# metric expansion, densities, tadpole symbol' + (' (covariant measure)' if cov else ''))
ME = U.metric_expansion()
sqrtg2, A2, sqrtgR2 = U.second_order_densities(ME)
sigma2 = U.tadpole_symbol(cov, ME)

# bases for the non-covariant fallbacks
htr = U.h_tr
b_h2 = htr**2
b_hh = sum(U.H[i, j]**2 for i in range(4) for j in range(4))
b1 = U.p2*htr**2
b2 = U.p2*b_hh
b3 = htr*sum(U.P[i]*U.P[j]*U.H[i, j] for i in range(4) for j in range(4))
b4 = sum(U.P[i]*U.P[j]*U.H[i, k]*U.H[j, k] for i in range(4) for j in range(4) for k in range(4))
fR, rR = U.fit(sqrtgR2, [b1, b2, b3, b4], ['p2_h2', 'p2_hh', 'h_pph', 'phhp'])
say(f'#   (sqrt g R)^(2) on {{p^2 h^2, p^2 h_ij h_ij, h p_i p_j h_ij, p_i h_ik h_jk p_j}} = {fR}, resid {rR}')
E_ = sp.symbols('e1:5', real=True)
def ward(form):
    """first-order variation under h -> h + p e + e p, identically in h and e"""
    return sp.expand(sum(sp.diff(form, U.H[i, j])*(U.P[i]*E_[j] + U.P[j]*E_[i]) for i in range(4) for j in range(i, 4))) == 0
say(f'#   Ward: (sqrt g R)^(2) transverse {ward(sqrtgR2)}; sqrt(g)^(2) transverse {ward(sqrtg2)} (it is not: the p^0 density is '
    f'invariant only together with the second-order gauge transformation of sqrt(g)^(1); at p^0 the test is the fit to sqrt(g)^(2) alone)')
ok &= ward(sqrtgR2)

def pieces(Qd, dimreg):
    """split a degree-d part into Lambda^4, Lambda^2, log Lambda^2, log m^2, (log mu), constant"""
    Qd = sp.expand(Qd)
    c4 = sp.expand(Qd.coeff(U.Lam, 4)); c2 = sp.expand(Qd.coeff(U.Lam, 2))
    rest = sp.expand(Qd - c4*U.Lam**4 - c2*U.Lam**2)
    cL = sp.expand(rest.coeff(U.lgL, 1)/2)
    rest = sp.expand(rest - 2*cL*U.lgL)
    cmu = sp.Integer(0)
    if dimreg:
        cmu = sp.expand(rest.coeff(sp.log(U.mu), 1)/2)
        rest = sp.expand(rest - 2*cmu*sp.log(U.mu))
    cm = sp.expand(rest.coeff(U.LM, 1))
    K = sp.expand(rest - cm*U.LM)
    assert not any(K.has(s) for s in (U.Lam, U.lgL, U.LM, U.epsinv, U.mu)) and not K.has(sp.log), 'leftover regulator symbols or logs in the constant'
    return dict(c4=c4, c2=c2, cL=cL, cmu=cmu, cm=cm, K=K)

def rep(label, form, basis, names, fallback=None, fnames=None, expect=None):
    sol, res = U.fit(form, basis, names)
    if sol is not None and res == 0:
        txt = ' + '.join(f'({sol[n]}) {n}' for n in names)
        say(f'#   {label}: {txt}' + (f'   [expected {expect}]' if expect else ''))
        return sol
    if form == 0:
        say(f'#   {label}: 0' + (f'   [expected {expect}]' if expect else ''))
        return {n: sp.Integer(0) for n in names}
    if fallback is not None:
        sol2, res2 = U.fit(form, fallback, fnames)
        say(f'#   {label}: NOT a combination of {names}; on {fnames}: {sol2}, resid {res2}' + (f'   [expected {expect}]' if expect else ''))
    else:
        say(f'#   {label}: NOT a combination of {names}' + (f'   [expected {expect}]' if expect else ''))
    return None

results = {}
for reg in args.regulators.split(','):
    say(f'# --- regulator: {reg} ---')
    t1 = time.time()
    dimreg = reg == 'dimreg'
    tab = U.bubble_massive(reg, covariant=cov, order=2)
    Qb = sp.Integer(0)
    for (i, j, r, s) in itertools.product(range(4), repeat=4):
        Qb += sp.Rational(1, 2)*U.H[i, j]*U.sym_get(tab, i, j, r, s)*U.H[r, s]
    if dimreg:
        perp = U.perp_symbol(cov, ME)
        fperp, _ = U.fit(perp, [sqrtg2], ['sqrtg2'])
        say(f'#   k_perp^2 coefficient of the trace term (second order): {fperp} x sqrt(g)^(2)   (check of the 2 Nov correction: plain 1, covariant 0)')
        ok &= fperp is not None and fperp['sqrtg2'] == (0 if cov else 1)
        Qt = U.tadpole_dimreg_full(sigma2, cov, ME)
    else:
        one = U.radial_one_prop if reg == 'cutoff' else U.pt_one_prop
        Qt = sp.expand(U.tadpole_from_symbol(sigma2, one(0), one(2)).subs({sp.log(m): U.LM/2, sp.log(m**2): U.LM}))
        Qt = sp.expand(sp.expand_log(Qt, force=True).subs({sp.log(m): U.LM/2}))
    Q = sp.expand(Qb + Qt)
    if dimreg:
        Qp = sp.expand(Q.coeff(U.epsinv, 1)); Qf = sp.expand(Q - Qp*U.epsinv)
        Q = sp.expand(-Qp*2*U.lgL + pi2*Qf)      # 1/eps_tHV -> -(1/16pi^2) log Lambda^2; units 1/16pi^2
    else:
        Q = sp.expand(pi2*Q)
    Q = sp.expand(Q.subs(U.Ssym, U.p2))
    say(f'#   ({time.time()-t1:.0f}s for the graphs)')
    res = {}
    # ---- p^0: cosmological order
    d0 = pieces(U.deg_part(Q, 0), dimreg)
    say('#   -- p^0 (cosmological order; units 1/16pi^2) --')
    rep('Lambda^4', d0['c4'], [sqrtg2], ['sqrtg2'], [b_h2, b_hh], ['h2', 'hh'], '-1/4 for proper time chi; non-covariant otherwise (31 Oct)')
    rep('Lambda^2', d0['c2'], [m**2*sqrtg2], ['m2sqrtg2'], [m**2*b_h2, m**2*b_hh], ['m2h2', 'm2hh'], '1/2 for proper time chi')
    sL0 = rep('log Lambda^2', d0['cL'], [m**4*sqrtg2], ['m4sqrtg2'], expect='-1/4 (P0)')
    if dimreg:
        say(f'#   log mu^2 coefficient equals c_L: {sp.expand(d0["cmu"] - d0["cL"]) == 0}')
        ok &= sp.expand(d0['cmu'] - d0['cL']) == 0
    say(f'#   P1 at p^0, c_m = -c_L: {sp.expand(d0["cm"] + d0["cL"]) == 0}')
    ok &= sp.expand(d0['cm'] + d0['cL']) == 0
    sK0 = rep('K_0 (constant)', d0['K'], [m**4*sqrtg2], ['m4sqrtg2'], [m**4*b_h2, m**4*b_hh], ['m4h2', 'm4hh'],
              'P2: covariant for dim reg and proper time chi only; P4: -3/8 in dim reg')
    ok &= sL0 is not None and sL0['m4sqrtg2'] == -sp.Rational(1, 4)
    res['K0'] = d0['K']; res['K0fit'] = sK0; res['cL0'] = d0['cL']
    # ---- p^2: Newton order
    d2 = pieces(U.deg_part(Q, 2), dimreg)
    say('#   -- p^2 (Newton order; units 1/16pi^2) --')
    say(f'#   Lambda^4 at p^2 is zero: {d2["c4"] == 0}')
    ok &= d2['c4'] == 0
    s2 = rep('Lambda^2', d2['c2'], [sqrtgR2], ['sqrtgR2'], [b1, b2, b3, b4], ['p2_h2', 'p2_hh', 'h_pph', 'phhp'],
             '-(1/6 - xi)/2 for proper time chi; non-covariant otherwise (31 Oct)')
    say(f'#   Lambda^2 piece transverse: {ward(d2["c2"])}')
    sL2 = rep('log Lambda^2', d2['cL'], [m**2*sqrtgR2], ['m2sqrtgR2'], expect='(1/6 - xi)/2 = -(6 xi - 1)/12 (P0)')
    say(f'#   log piece transverse: {ward(d2["cL"])}')
    if dimreg:
        say(f'#   log mu^2 coefficient equals c_L: {sp.expand(d2["cmu"] - d2["cL"]) == 0}')
        ok &= sp.expand(d2['cmu'] - d2['cL']) == 0
    say(f'#   P1 at p^2, c_m = -c_L: {sp.expand(d2["cm"] + d2["cL"]) == 0}')
    ok &= sp.expand(d2['cm'] + d2['cL']) == 0
    wK = ward(d2['K'])
    sK2 = rep('K_2 (constant)', d2['K'], [m**2*sqrtgR2], ['m2sqrtgR2'], [m**2*b1, m**2*b2, m**2*b3, m**2*b4],
              ['m2p2_h2', 'm2p2_hh', 'm2h_pph', 'm2phhp'], 'P2: transverse for dim reg and proper time chi only; P4: (1/6 - xi)/2 in dim reg')
    say(f'#   K_2 transverse (first-order variation under h -> h + p e + e p vanishes identically): {wK}')
    ok &= sL2 is not None and sp.simplify(sL2['m2sqrtgR2'] - (sp.Rational(1, 6) - xi)/2) == 0
    res['K2'] = d2['K']; res['K2fit'] = sK2; res['ward2'] = wK; res['cL2'] = d2['cL']; res['c2'] = d2['c2']; res['c2fit'] = s2
    # ---- p^4: only the log, as a check of the expansion against the 31 Oct logs
    d4 = pieces(U.deg_part(Q, 4), dimreg)
    say('#   -- p^4 (check only) --')
    s4 = rep('log Lambda^2', d4['cL'], [U.Q_Ric, U.Q_R2], ['Ric2', 'R2'], expect='-1/120, -(60 xi^2 - 20 xi + 1)/240 (P0)')
    ok &= s4 is not None and sp.simplify(s4['Ric2'] + sp.Rational(1, 120)) == 0 and sp.simplify(s4['R2'] + (60*xi**2 - 20*xi + 1)/240) == 0
    say(f'#   P1 at p^4, c_m = -c_L: {sp.expand(d4["cm"] + d4["cL"]) == 0}   (the massive p^4 constant is the low-energy expansion, not reported)')
    ok &= sp.expand(d4['cm'] + d4['cL']) == 0
    results[reg] = res
    say(f'#   ({time.time()-t1:.0f}s)')

say('# --- P2: covariance/transversality pattern of the constants ---')
expect_t = {'dimreg': True, 'propertime': cov, 'cutoff': False}
p2ok = True
for reg, r in results.items():
    t0_ = r['K0fit'] is not None; t2_ = r['ward2']
    say(f'#   [{reg}, {args.measure}] K_0 covariant: {t0_}; K_2 transverse: {t2_}   [P2 expects {expect_t[reg]} for both]')
    p2ok &= (t0_ == expect_t[reg]) and (t2_ == expect_t[reg])
say(f'#   P2 pattern as predicted: {p2ok}')

say('# --- P3, P4: the transverse constants against the heat-kernel/MS-bar values ---')
p3ok = True; p4ok = True
if 'dimreg' in results and results['dimreg']['K0fit'] is not None and results['dimreg']['K2fit'] is not None:
    k0 = results['dimreg']['K0fit']['m4sqrtg2']; k2 = results['dimreg']['K2fit']['m2sqrtgR2']
    say(f'#   K_0[dimreg] = {k0}   [P4: -3/8];   K_2[dimreg] = {sp.factor(k2)}   [P4: (1/6 - xi)/2]')
    p4ok &= sp.simplify(k0 + sp.Rational(3, 8)) == 0 and sp.simplify(k2 - (sp.Rational(1, 6) - xi)/2) == 0
elif 'dimreg' in results:
    p4ok = False
    say('#   P4 FAILS: the dim-reg constants are not covariant')
else:
    p4ok = None
    say('#   P4 not scorable: dim reg not in this run')
say(f'#   P4 holds: {p4ok}')
if 'dimreg' in results and 'propertime' in results and results['propertime']['ward2'] and results['propertime']['K0fit'] is not None:
    for d_, key, cl in ((0, 'K0', 'cL0'), (2, 'K2', 'cL2')):
        diff = sp.expand(results['propertime'][key] - results['dimreg'][key] + U.EulerGamma*results['dimreg'][cl])
        say(f'#   K_{d_}[propertime] - K_{d_}[dimreg] + gamma_E c_L^({d_}) = {diff}   [P3: 0]')
        p3ok &= diff == 0
    say(f'#   P3 holds: {p3ok}')
else:
    say('#   P3 not scorable (needs dim reg and a transverse proper-time constant: expected only with the covariant measure)')
    p3ok = None

say('# --- P5: the Sakharov reading (printed, not interpreted further) ---')
MP = sp.Symbol('M_P', positive=True); LamS = sp.Symbol('Lambda', positive=True); muS = sp.Symbol('mu', positive=True)
for reg, r in results.items():
    if r['K2fit'] is None:
        say(f'#   [{reg}] K_2 not transverse: no covariant induced M_P^2 can be read off (the -1/12 discipline fails at the first slot)')
        continue
    c2v = r['c2fit']['sqrtgR2'] if r['c2fit'] is not None else sp.Symbol('c2_noncov')
    cLv = (sp.Rational(1, 6) - xi)/2
    K2v = r['K2fit']['m2sqrtgR2']
    if reg == 'dimreg':
        ind = -2*(cLv*m**2*sp.log(muS**2/m**2) + K2v*m**2)/pi2
    else:
        ind = -2*(c2v*LamS**2 + cLv*m**2*sp.log(LamS**2/m**2) + K2v*m**2)/pi2
    say(f'#   [{reg}] M_P^2(induced) = {sp.collect(sp.expand(sp.expand_log(sp.simplify(ind), force=True)), [sp.log(LamS), sp.log(muS), sp.log(m)])}')
    say(f'#   [{reg}] d M_P^2(induced) / d log(scale^2), from the m^2 log term alone = {sp.simplify(-2*cLv*m**2/pi2)}   (universal: sign of (1/6 - xi), zero at xi = 1/6)')
say('#   Lambda^2 term: absent in dim reg, non-covariant under the hard cutoff -> fails the first slot; m^2 log term: same coefficient in all'
    ' regulators and measures -> passes; its constant K_2 is scheme-dependent -> the induced value is not fixed, only its running.')

all_t = all(r['ward2'] and r['K0fit'] is not None for r in results.values())
nont = [k for k, r in results.items() if not (r['ward2'] and r['K0fit'] is not None)]
scored = ok and p2ok and (p4ok in (True, None)) and (p3ok in (True, None))
status = 'OK' if (scored and all_t) else (f'LOG-OK_FINITE-NONTRANSVERSE[{",".join(nont)}]' if scored else 'MISMATCH')
say(f'# END {name} elapsed={time.time()-t0:.1f}s status={status}')
if out: out.close()
