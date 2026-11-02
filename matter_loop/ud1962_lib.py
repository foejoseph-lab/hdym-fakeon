"""
Shared machinery for rung 1 (Utiyama-DeWitt 1962, reconstructed): one real scalar on g = delta + h, Euclidean,
flat-space perturbation theory to second order in h, one loop, hard momentum cutoff Lambda.

Conventions (identical to ud1962_bubble.py, which this supersedes):
  S = 1/2 int sqrt(g) [ g^{mu nu} d_mu phi d_nu phi + m^2 phi^2 + xi R phi^2 ],  R = +2 on the unit S^2.
  Fourier f(x) = int_p f(p) e^{ipx}.  Quadratic forms are written  Gamma_2 = int_p (1/2) h(p) . Pi(p) . h(-p)
  and, in the real-symbol bookkeeping used throughout, h(-p) is represented by the same symbols h_ij as h(p);
  every form we match is even in p, so this is consistent.
"""
import itertools
import sympy as sp

m, xi, Lam, x = sp.symbols('m xi Lambda x', positive=True)
P = sp.symbols('p1:5', real=True)
L = sp.symbols('l1:5', real=True)
X = sp.symbols('x1:5', real=True)
H = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'h{min(i,j)+1}{max(i,j)+1}', real=True))
HSYMS = [H[i, j] for i in range(4) for j in range(i, 4)]
p2 = sum(pi**2 for pi in P)
delta = sp.eye(4)
lsym = sp.Symbol('ell', positive=True)
Dl = sp.Symbol('Delta', positive=True)
lgL = sp.log(Lam)


def dot(a, b):
    return sum(ai*bi for ai, bi in zip(a, b))


# ---------------------------------------------------------------- vertices
def V1(k, q, p):
    """first-order vertex: S_1 = 1/2 int phi(k) phi(q) h_{mu nu}(p) V_{mu nu}(k,q;p), k+q+p=0"""
    kq = dot(k, q)
    pp = dot(p, p)
    return sp.Matrix(4, 4, lambda i, j:
        sp.Rational(1, 2)*(k[i]*q[j] + k[j]*q[i])
        - sp.Rational(1, 2)*delta[i, j]*(kq - m**2)
        + xi*(pp*delta[i, j] - p[i]*p[j]))


def V1cov(k, q, p):
    """first-order vertex for the covariant measure, phi = g^(-1/4) chi: the operator becomes
    g^(-1/4) K g^(-1/4) = K_0 + [K_1 - (h/4) K_0 - K_0 (h/4)] + O(h^2), h = trace;  K_0(k) = k^2 + m^2."""
    V = V1(k, q, p)
    corr = sp.Rational(1, 4)*((dot(k, k) + m**2) + (dot(q, q) + m**2))
    return V - corr*delta


# ---------------------------------------------------------------- angular averages in d=4
def _dfact(n):
    r = 1
    for i in range(n - 1, 0, -2):
        r *= i
    return r


def ang_avg_monomial(exps):
    if any(e % 2 for e in exps):
        return sp.Integer(0)
    A = sum(exps) // 2
    num = 1
    for e in exps:
        num *= _dfact(e)
    den = 1
    for j in range(A):
        den *= (4 + 2*j)
    return sp.Rational(num, den)


def ang_avg(expr, vars_=L):
    poly = sp.Poly(sp.expand(expr), *vars_)
    res = sp.Integer(0)
    for monom, coeff in poly.terms():
        c = ang_avg_monomial(monom)
        if c != 0:
            res += coeff*c*lsym**sum(monom)
    return sp.expand(res)


# ---------------------------------------------------------------- radial integrals, hard cutoff
def radial_two_prop(power):
    """(1/8pi^2) int_0^Lam dl l^(3+power)/(l^2+Delta)^2, large-Lambda, keeping Lam^4, Lam^2, log, const.
    Elementary: t=l^2, T=Lam^2, (1/2) int_0^T t^(1+k)/(t+Delta)^2 dt."""
    T = Lam**2
    lg = 2*lgL - sp.log(Dl)
    if power == 0:
        ser = sp.Rational(1, 2)*(lg - 1)
    elif power == 2:
        ser = sp.Rational(1, 2)*(T - 2*Dl*lg + Dl)
    elif power == 4:
        ser = sp.Rational(1, 2)*(T**2/2 - 2*Dl*T + 3*Dl**2*lg - Dl**2)
    else:
        raise ValueError(power)
    return sp.expand(ser/(8*sp.pi**2))


def radial_one_prop(power):
    """(1/8pi^2) int_0^Lam dl l^(3+power)/(l^2+m^2), large-Lambda.  (1/2) int_0^T t^(1+k)/(t+m^2) dt."""
    T = Lam**2
    M2 = m**2
    lg = 2*lgL - sp.log(M2)
    if power == 0:
        ser = sp.Rational(1, 2)*(T - M2*lg)
    elif power == 2:
        ser = sp.Rational(1, 2)*(T**2/2 - M2*T + M2**2*lg)
    else:
        raise ValueError(power)
    return sp.expand(ser/(8*sp.pi**2))


def check_radials(say):
    """numerical check of the large-Lambda forms against the exact integrals"""
    for pw in (0, 2, 4):
        exact = sp.integrate(lsym**(3 + pw)/(lsym**2 + Dl)**2, (lsym, 0, Lam))
        ser = radial_two_prop(pw)*8*sp.pi**2
        sub = {Lam: 1000, Dl: sp.Rational(7, 10)}
        rel = abs((sp.N(exact.subs(sub), 30) - sp.N(ser.subs(sub), 30))/sp.N(exact.subs(sub), 30))
        say(f'#   radial two-prop k={pw//2}: rel.diff vs exact at Lambda=1e3, Delta=0.7: {sp.N(rel, 3)}')
    for pw in (0, 2):
        exact = sp.integrate(lsym**(3 + pw)/(lsym**2 + m**2), (lsym, 0, Lam))
        ser = radial_one_prop(pw)*8*sp.pi**2
        sub = {Lam: 1000, m: sp.Rational(4, 5)}
        rel = abs((sp.N(exact.subs(sub), 30) - sp.N(ser.subs(sub), 30))/sp.N(exact.subs(sub), 30))
        say(f'#   radial one-prop k={pw//2}: rel.diff vs exact at Lambda=1e3, m=0.8: {sp.N(rel, 3)}')


# ---------------------------------------------------------------- splitting by divergence order
def split(expr):
    """dict L4, L2, log (coefficient of log Lambda^2), fin"""
    e = sp.expand(expr)
    d = {}
    d['L4'] = e.coeff(Lam, 4)
    d['L2'] = e.coeff(Lam, 2)
    rest = sp.expand(e - d['L4']*Lam**4 - d['L2']*Lam**2)
    d['log'] = sp.expand(rest.coeff(lgL, 1)/2)
    d['fin'] = sp.expand(rest - 2*d['log']*lgL)
    assert not d['fin'].has(Lam), 'leftover Lambda dependence'
    return d


# ---------------------------------------------------------------- the bubble
def bubble_raw(radial_fn, vertex=None):
    """Pi_{ij rs}(p) for i<=j, r<=s as dict of raw expressions in x (Feynman parameter), with
    1/((k^2+m^2)((k-p)^2+m^2)) -> int_0^1 dx [radial_fn(power)] where radial_fn(2n) stands for
    int d^4l/(2pi)^4 l^(2n) * (regularized two-propagator weight) with Delta -> m^2 + x(1-x)p^2.
    The -1/2 of Pi is included; the x-integral is NOT done."""
    Delta = m**2 + x*(1 - x)*p2
    k = [li + x*pi for li, pi in zip(L, P)]
    kmp = [ki - pi for ki, pi in zip(k, P)]
    vertex = vertex or V1
    Va = vertex([-ki for ki in k], kmp, P)
    Vb = vertex([-ki for ki in kmp], k, [-pi for pi in P])
    R = {pw: radial_fn(pw) for pw in (0, 2, 4)}
    out = {}
    for (i, j, r, s) in itertools.product(range(4), repeat=4):
        if i > j or r > s:
            continue
        num = ang_avg(Va[i, j]*Vb[r, s])
        poly = sp.Poly(num, lsym)
        val = sp.Integer(0)
        for (pw,), c in poly.terms():
            val += c*R[pw]
        out[(i, j, r, s)] = sp.expand(-sp.Rational(1, 2)*val.subs(Dl, Delta))
    return out


def bubble_Pi(radial_fn=radial_two_prop, vertex=None):
    """as bubble_raw, then split by divergence order and x-integrated (L4, L2, log); fin kept with x."""
    raw = bubble_raw(radial_fn, vertex)
    out = {}
    for key, val in raw.items():
        d = split(val)
        out[key] = {
            'L4': sp.expand(sp.integrate(d['L4'], (x, 0, 1))),
            'L2': sp.expand(sp.integrate(d['L2'], (x, 0, 1))),
            'log': sp.expand(sp.integrate(d['log'], (x, 0, 1))),
            'fin': d['fin'],
        }
    return out


def sym_get(table, i, j, r, s):
    a, b = (i, j) if i <= j else (j, i)
    c, d_ = (r, s) if r <= s else (s, r)
    return table[(a, b, c, d_)]


def quad_form(table, part):
    """Gamma_2 contribution = 1/2 sum h_ij Pi_ijrs h_rs for the given part"""
    G = sp.Integer(0)
    for (i, j, r, s) in itertools.product(range(4), repeat=4):
        G += sp.Rational(1, 2)*H[i, j]*sym_get(table, i, j, r, s)[part]*H[r, s]
    return sp.expand(G)


# ---------------------------------------------------------------- the tadpole (seagull)
_Hp = sp.Matrix(4, 4, lambda i, j: sp.Symbol(f'hp{min(i,j)+1}{max(i,j)+1}'))
_eps = sp.Symbol('epsilon')
KS = sp.symbols('k1:5', real=True)


def _trunc(expr, order=2):
    """keep terms up to _eps^order"""
    e = sp.expand(expr)
    return sum(e.coeff(_eps, n)*_eps**n for n in range(order + 1))


def metric_expansion():
    """x-dependent expansion of the metric objects to second order in eps, h(x) = Hp e^{ipx} + H e^{-ipx}.
    Returns dict with: hx, trh, sqrtg (1+eps..+eps^2..), G = sqrt(g) g^{mu nu}, sqrtgR = sqrt(g) R, gm14 = g^(-1/4)."""
    px = dot(P, X)
    ep, em = sp.exp(sp.I*px), sp.exp(-sp.I*px)
    hx = (_Hp*ep + H*em)
    D = lambda f, a: sp.diff(f, X[a])
    Gam1 = [[[sp.Rational(1, 2)*(D(hx[a, c], b) + D(hx[a, b], c) - D(hx[b, c], a))
              for c in range(4)] for b in range(4)] for a in range(4)]
    Gam2 = [[[sp.expand(-sum(hx[a, d]*Gam1[d][b][c] for d in range(4)))
              for c in range(4)] for b in range(4)] for a in range(4)]
    Ric1 = sp.Matrix(4, 4, lambda b, c: sp.expand(sum(D(Gam1[a][b][c], a) - D(Gam1[a][b][a], c) for a in range(4))))
    Ric2 = sp.Matrix(4, 4, lambda b, c: sp.expand(
        sum(D(Gam2[a][b][c], a) - D(Gam2[a][b][a], c) for a in range(4))
        + sum(Gam1[a][a][d]*Gam1[d][b][c] - Gam1[a][c][d]*Gam1[d][b][a] for a in range(4) for d in range(4))))
    R1x = Ric1.trace()
    R2x = sp.expand(Ric2.trace() - sum(hx[b, c]*Ric1[b, c] for b in range(4) for c in range(4)))
    trh = hx.trace()
    trhh = (hx*hx).trace()
    sqrtg = 1 + _eps*trh/2 + _eps**2*(trh**2/8 - trhh/4)
    sqrtgR = _eps*R1x + _eps**2*sp.expand(R2x + trh*R1x/2)
    G = (sp.eye(4) + _eps*(-hx + trh/2*sp.eye(4)) + _eps**2*((trh**2/8 - trhh/4)*sp.eye(4) - (trh/2)*hx + hx*hx))
    # g^(-1/4) = exp(-1/4 ln det g),  ln det g = eps trh - eps^2 trhh/2
    gm14 = 1 - _eps*trh/4 + _eps**2*(trh**2/32 + trhh/8)
    return dict(hx=hx, trh=trh, sqrtg=sqrtg, sqrtgR=sqrtgR, G=G, gm14=gm14, ep=ep, em=em)


def _cross(expr, ME, std=True):
    """x-independent cross term (coefficient of Hp*H) of expr, then Hp -> H; std: divide by 2 (see tadpole_form)."""
    ep, em = ME['ep'], ME['em']
    e = sp.expand(expr)
    e = e.subs({ep**2: 0, em**2: 0}).subs({ep*em: 1})
    e = sp.expand(e)
    e = sum(t for t in sp.Add.make_args(e) if not t.has(*X))
    e = e.subs({_Hp[i, j]: H[i, j] for i in range(4) for j in range(4)})
    e = sp.expand(e)
    return sp.expand(e/2) if std else e


def second_order_densities(ME=None):
    """standard second-order densities sqrt(g)^(2), (sqrt(g) g^{mu nu})^(2), (sqrt(g) R)^(2) as quadratic forms in h."""
    ME = ME or metric_expansion()
    c2 = lambda e: _cross(sp.expand(e).coeff(_eps, 2), ME)
    return c2(ME['sqrtg']), ME['G'].applyfunc(c2), c2(ME['sqrtgR'])


def tadpole_symbol(covariant_measure=False, ME=None):
    """Symbol sigma_2(k) of the second-order part of the operator (K, or g^(-1/4) K g^(-1/4) with the covariant
    measure), x-averaged (standard normalization), as a polynomial in k1..k4:
        Tr(w(K_0) O_2) = int d^4x int_k w(k^2+m^2) sigma_2(x,k),  sigma(x,k) = e^{-ikx} O (e^{ikx}).
    K f = -d_mu (G_{mu nu} d_nu f) + sqrt(g)(m^2 + xi R) f."""
    ME = ME or metric_expansion()
    kx = dot(KS, X)
    chi = sp.exp(sp.I*kx)
    f = _trunc(ME['gm14']*chi) if covariant_measure else chi
    G = ME['G']
    Kf = 0
    for mu in range(4):
        inner = sum(G[mu, nu]*sp.diff(f, X[nu]) for nu in range(4))
        Kf += -sp.diff(_trunc(inner), X[mu])
    Kf += (ME['sqrtg']*(m**2) + xi*ME['sqrtgR'])*f
    Kf = _trunc(Kf)
    if covariant_measure:
        Kf = _trunc(ME['gm14']*Kf)
    sigma = sp.expand(Kf/chi)
    sigma = sp.expand(sp.powsimp(sigma))
    s2 = sigma.coeff(_eps, 2)
    return _cross(s2, ME)


def tadpole_from_symbol(sigma2, J0, J2):
    """Q_tad = 1/2 int_k w(k) sigma2(k) with sigma2 polynomial in k of degree <= 2: angular average,
    k^0 -> J0, k^2 -> J2."""
    poly = sp.Poly(sp.expand(sigma2), *KS)
    val = 0
    for monom, coeff in poly.terms():
        c = ang_avg_monomial(monom)
        if c == 0:
            continue
        deg = sum(monom)
        val += coeff*c*(J0 if deg == 0 else J2 if deg == 2 else None)
    return sp.expand(sp.Rational(1, 2)*val)


def tadpole_raw(J0, J2, dens=None):
    """Q_tad = 1/2 [ A2_{mu mu} J2/4 + B2 J0 ] with given one-propagator integrals J0 = int_k w(k), J2 = int_k k^2 w(k)."""
    sqrtg2, A2, sqrtgR2 = dens if dens is not None else second_order_densities()
    B2 = m**2*sqrtg2 + xi*sqrtgR2
    return sp.expand(sp.Rational(1, 2)*(A2.trace()*J2/4 + B2*J0)), (sqrtg2, A2, sqrtgR2)


def tadpole_form():
    """Gamma_tad as a quadratic form 1/2 h Pi_tad h, split by divergence order.
    Gamma_tad = 1/2 Tr(G K_2) = 1/2 int d^4x int_k [A2_{mu nu}(x) k_mu k_nu + B2(x)]/(k^2+m^2);
    in the convention Gamma_2 = int_p 1/2 h Pi h, the quadratic form Q = 1/2 h Pi h for a local term
    c int d^4x F^(2)[h] equals c * F^(2)_std[h]; and Q_tad = 1/2 [ A2_{mu mu} J2/4 + B2 J0 ] with the
    standard (halved) densities, J2 = int_k k^2/(k^2+m^2), J0 = int_k 1/(k^2+m^2), k_mu k_nu -> delta k^2/4."""
    G, dens = tadpole_raw(radial_one_prop(0), radial_one_prop(2))
    return split(G), dens


# ---------------------------------------------------------------- proper-time (1965) regularization
EulerGamma = sp.EulerGamma
tsym = sp.Symbol('t', positive=True)


def pt_two_prop(power):
    """int d^4l/(2pi)^4 l^(2n) int_{1/Lam^2}^inf dt t e^{-t Delta}/... : the Schwinger form of the two-propagator
    product after Feynman parametrization (u = x t), Gaussian in l:
        int d^4l/(2pi)^4 l^(2n) e^{-t l^2} = Gamma(n+2)/(16 pi^2 t^(n+2)),
    so the weight is Gamma(n+2)/(16 pi^2) int_{1/Lam^2}^inf dt e^{-t Delta} t^(-n-1) = Gamma(n+2) Delta^n Gamma(-n, Delta/Lam^2)/(16 pi^2).
    Large-Lambda expansions of the incomplete gamma, O(1/Lam^2) dropped."""
    n = power // 2
    lg = 2*lgL - sp.log(Dl)
    g = EulerGamma
    if n == 0:
        G = lg - g
    elif n == 1:
        G = Lam**2 - Dl*lg + Dl*(g - 1)
    elif n == 2:
        G = Lam**4/2 - Dl*Lam**2 + Dl**2*lg/2 + Dl**2*(sp.Rational(3, 4) - g/2)
    else:
        raise ValueError(power)
    return sp.expand(sp.factorial(n + 1)*G/(16*sp.pi**2))


def pt_two_prop_integrand(power):
    """the same before the t-integration: Gamma(n+2) e^{-t Delta} / (16 pi^2 t^(n+1)) -- for the small-t expansion
    of Tr e^{-tD} (the heat-kernel coefficients)."""
    n = power // 2
    return sp.factorial(n + 1)*sp.exp(-tsym*Dl)/(16*sp.pi**2*tsym**(n + 1))


def pt_one_prop(power):
    """int d^4k/(2pi)^4 k^(2n) int_{1/Lam^2}^inf dt e^{-t(k^2+m^2)} = Gamma(n+2)/(16pi^2) int dt e^{-t m^2} t^(-n-2)"""
    n = power // 2
    lg = 2*lgL - sp.log(m**2)
    g = EulerGamma
    M = m**2
    if n == 0:
        G = Lam**2 - M*lg + M*(g - 1)
    elif n == 1:
        G = Lam**4/2 - M*Lam**2 + M**2*lg/2 + M**2*(sp.Rational(3, 4) - g/2)
    else:
        raise ValueError(power)
    return sp.expand(sp.factorial(n + 1)*G/(16*sp.pi**2))


def pt_one_prop_integrand(power):
    n = power // 2
    return sp.factorial(n + 1)*sp.exp(-tsym*m**2)/(16*sp.pi**2*tsym**(n + 2))


def check_pt_radials(say):
    for pw in (0, 2, 4):
        n = pw // 2
        exact = sp.factorial(n + 1)*Dl**n*sp.uppergamma(-n, Dl/Lam**2)/(16*sp.pi**2)
        sub = {Lam: 1000, Dl: sp.Rational(7, 10)}
        ex = sp.N(exact.subs(sub), 30); se = sp.N(pt_two_prop(pw).subs(sub), 30)
        say(f'#   proper-time two-prop n={n}: rel.diff vs exact at Lambda=1e3, Delta=0.7: {sp.N(abs((ex-se)/ex), 3)}')
    for pw in (0, 2):
        n = pw // 2
        exact = sp.factorial(n + 1)*(m**2)**(n + 1)*sp.uppergamma(-n - 1, m**2/Lam**2)/(16*sp.pi**2)
        sub = {Lam: 1000, m: sp.Rational(4, 5)}
        ex = sp.N(exact.subs(sub), 30); se = sp.N(pt_one_prop(pw).subs(sub), 30)
        say(f'#   proper-time one-prop n={n}: rel.diff vs exact at Lambda=1e3, m=0.8: {sp.N(abs((ex-se)/ex), 3)}')


# ---------------------------------------------------------------- first-order curvature forms
h_tr = sum(H[i, i] for i in range(4))
R1 = sp.expand(-sum(P[i]*P[j]*H[i, j] for i in range(4) for j in range(4)) + p2*h_tr)
R1mn = sp.Matrix(4, 4, lambda i, j: sp.Rational(1, 2)*(
    -sum(P[a]*P[j]*H[a, i] for a in range(4)) - sum(P[a]*P[i]*H[a, j] for a in range(4))
    + P[i]*P[j]*h_tr + p2*H[i, j]))
Q_R2 = sp.expand(R1*R1)
Q_Ric = sp.expand(sum(R1mn[i, j]**2 for i in range(4) for j in range(4)))


def deg_part(expr, deg):
    poly = sp.Poly(sp.expand(expr), *P)
    return sp.expand(sum(c*sp.prod([Pi_**e for Pi_, e in zip(P, mon)]) for mon, c in poly.terms() if sum(mon) == deg))


def fit(expr, basis, names):
    """fit expr = sum c_i basis_i identically in (p, h); return dict or None, plus residual"""
    cs = sp.symbols(' '.join(names)) if len(names) > 1 else (sp.Symbol(names[0]),)
    resid = sp.expand(expr - sum(c*b for c, b in zip(cs, basis)))
    eqs = sp.Poly(resid, *P, *HSYMS).coeffs()
    sol = sp.solve(eqs, list(cs), dict=True)
    if not sol:
        return None, resid
    return {n: sp.factor(sol[0][c]) for n, c in zip(names, cs)}, sp.expand(resid.subs(sol[0]))


# ---------------------------------------------------------------- dimensional regularization (1972-74)
epsinv = sp.Symbol('epsinv')   # 1/eps in 't Hooft-Veltman's convention, 1/eps = 1/(8 pi^2 (n-4))


def _pole(expr_n, nsym):
    """pole part at n -> 4 of expr_n(n), expressed as coefficient * epsinv with epsinv = 1/(8 pi^2 (n-4))"""
    d = sp.Symbol('d')
    ser = sp.series(expr_n.subs(nsym, 4 + d), d, 0, 1).removeO()
    c = sp.expand(ser).coeff(d, -1)
    return sp.expand(c*8*sp.pi**2*epsinv)


def dimreg_two_prop(power):
    """pole part of int d^n l/(2pi)^n (l^2)^a/(l^2+Delta)^2, a = power/2, in units of epsinv.
    Angular averages are taken in d=4 (their O(n-4) parts only affect finite terms)."""
    a = power // 2
    n = sp.Symbol('n')
    expr = sp.gamma(a + n/2)*sp.gamma(2 - a - n/2)/(sp.gamma(n/2)*sp.gamma(2))*Dl**(a + n/2 - 2)/(4*sp.pi)**(n/2)
    return _pole(expr, n)


def dimreg_one_prop(power):
    """pole part of int d^n k/(2pi)^n (k^2)^a/(k^2+m^2), in units of epsinv."""
    a = power // 2
    n = sp.Symbol('n')
    expr = sp.gamma(a + n/2)*sp.gamma(1 - a - n/2)/(sp.gamma(n/2)*sp.gamma(1))*(m**2)**(a + n/2 - 1)/(4*sp.pi)**(n/2)
    return _pole(expr, n)


# ================================================================ additions of 1 Nov 2026 (rung 1 finish)
# ---------------------------------------------------------------- hard cutoff on the UNSHIFTED loop momentum
Tsym = sp.Symbol('T', positive=True)        # T = Lambda^2
_radial_cache = {}


def radial_hard(a, n, M2):
    """(1/8pi^2) int_0^Lam dk k^(3+2a)/(k^2+M2)^n  =  (1/16pi^2) int_0^T dt t^(1+a)/(t+M2)^n,  T = Lam^2,
    large-T asymptotics keeping T^2, T, log T, T^0 and dropping O(1/T).  Exact elementary integral by t = s - M2:
    int_{M2}^{T+M2} (s-M2)^(1+a) s^(-n) ds = sum_k C(1+a,k) (-M2)^(1+a-k) int s^(k-n) ds."""
    key = (a, n, M2)
    if key in _radial_cache:
        return _radial_cache[key]
    s = sp.Symbol('s', positive=True)
    res = sp.Integer(0)
    for k in range(a + 2):
        c = sp.binomial(a + 1, k)*(-M2)**(a + 1 - k)
        q = k - n
        if q == -1:
            # log(T+M2) - log(M2) = log T + log(1 + M2/T) - log M2 -> log T - log M2 + O(1/T)
            term = 2*lgL - sp.log(M2) if M2 != 0 else 2*lgL
        else:
            # [(T+M2)^(q+1) - M2^(q+1)]/(q+1); expand (T+M2)^(q+1) for large T to the order needed
            if M2 == 0:
                term = Tsym**(q + 1)/(q + 1)
            else:
                ser = sp.series((Tsym + M2)**(q + 1), Tsym, sp.oo, 4).removeO() if q + 1 < 0 else sp.expand((Tsym + M2)**(q + 1))
                term = (ser - M2**(q + 1))/(q + 1)
        res += c*term
    res = sp.expand(res)

    def _nonneg_in(expr, sym):
        """drop the terms of an expanded sum that carry a negative power of sym"""
        kept_ = sp.Integer(0)
        for t_ in sp.Add.make_args(expr):
            powers = [e_.exp for e_ in t_.atoms(sp.Pow) if e_.base == sym]
            if any(pw_ < 0 for pw_ in powers):
                continue
            kept_ += t_
        return kept_

    kept = _nonneg_in(res, Tsym)
    kept = sp.expand(kept.subs(Tsym, Lam**2)/(16*sp.pi**2))
    kept = sp.expand(_nonneg_in(kept, Lam))
    _radial_cache[key] = kept
    return kept


def check_radial_hard(say):
    """numerical check of radial_hard against the exact integral for the (a, n) pairs the unshifted bubble uses"""
    worst = 0
    for a in range(0, 5):
        for n in range(2, 7):
            exact = sp.integrate(lsym**(3 + 2*a)/(lsym**2 + m**2)**n, (lsym, 0, Lam))/(8*sp.pi**2)
            sub = {Lam: 1000, m: sp.Rational(4, 5)}
            ex = sp.N(exact.subs(sub), 40); se = sp.N(radial_hard(a, n, m**2).subs(sub), 40)
            if ex != 0:
                rel = abs((ex - se)/ex)
                worst = max(worst, rel)
    say(f'#   radial_hard: worst rel.diff vs exact over a=0..4, n=2..6 at Lambda=1e3, m=0.8: {sp.N(worst, 3)}')


_eps_p = sp.Symbol('eps_p')


def bubble_raw_unshifted(vertex=None, pmax=4):
    """Pi_{ij rs}(p) with the hard cutoff |k| < Lambda on the propagator momentum k itself (no Feynman parameter).
    1/((k^2+m^2)((k-p)^2+m^2)) = (k^2+m^2)^-2 sum_j (-u)^j,  u = (p^2 - 2 k.p)/(k^2+m^2); every power of p lowers
    the UV degree by one, the p^0 term has UV degree 4, so terms of p-degree > 4 carry no Lambda dependence and are
    dropped (pmax).  Exact in Lambda^4, Lambda^2, log Lambda^2; the finite part from this expansion is NOT the
    finite part of the integral (the large-k expansion is not uniform in the infrared) and is not reported."""
    vertex = vertex or V1
    k = list(L)
    kmp = [ki - pi for ki, pi in zip(k, P)]
    Va = vertex([-ki for ki in k], kmp, P)
    Vb = vertex([-ki for ki in kmp], k, [-pi for pi in P])
    k2 = dot(k, k)
    kp = dot(k, P)
    M2 = m**2
    Dinv = sp.Symbol('Dinv')      # stands for 1/(k^2+m^2)
    u = (_eps_p**2*p2 - 2*_eps_p*kp)*Dinv
    geo = sp.expand(sum((-u)**j for j in range(pmax + 1)))
    geo = sum(geo.coeff(_eps_p, d)*_eps_p**d for d in range(pmax + 1))
    out = {}
    for (i, j, r, s) in itertools.product(range(4), repeat=4):
        if i > j or r > s:
            continue
        num = sp.expand(Va[i, j]*Vb[r, s]).subs({pi: _eps_p*pi for pi in P})
        num = sp.expand(num)
        num = sum(num.coeff(_eps_p, d)*_eps_p**d for d in range(pmax + 1))
        tot = sp.expand(num*geo*Dinv**2)
        tot = sum(tot.coeff(_eps_p, d)*_eps_p**d for d in range(pmax + 1))
        tot = sp.expand(tot.subs(_eps_p, 1))
        # angular average over k-hat, then radial by (a, n)
        poly = sp.Poly(tot, *L, Dinv)
        val = sp.Integer(0)
        for monom, coeff in poly.terms():
            kexps, n = monom[:4], monom[4]
            c = ang_avg_monomial(kexps)
            if c == 0:
                continue
            a = sum(kexps)//2
            val += coeff*c*radial_hard(a, n, M2)
        out[(i, j, r, s)] = sp.expand(-sp.Rational(1, 2)*val)
    return out


# ---------------------------------------------------------------- spin projectors (Barnes-Rivers), quadratic forms
def _theta(i, j):
    return delta[i, j] - P[i]*P[j]/p2


def projector_forms():
    """p^4 * (h . P2 . h) and p^4 * (h . P0s . h) as polynomials in (p, h);  theta = delta - p p/p^2,
    P2 = 1/2 (theta theta + theta theta) - 1/3 theta theta,  P0s = 1/3 theta_{mu nu} theta_{rho sigma}."""
    th = sp.Matrix(4, 4, lambda i, j: _theta(i, j))
    hth = sum(H[i, j]*th[i, j] for i in range(4) for j in range(4))             # theta . h
    hthth = sum(H[i, j]*th[i, r]*th[j, s]*H[r, s] for i in range(4) for j in range(4) for r in range(4) for s in range(4))
    P0s = sp.Rational(1, 3)*hth**2
    P2 = hthth - P0s
    f2 = sp.expand(sp.cancel(sp.together(p2**2*P2)))
    f0 = sp.expand(sp.cancel(sp.together(p2**2*P0s)))
    return f2, f0


# ---------------------------------------------------------------- massless finite part, three regulators
Ssym = sp.Symbol('S', positive=True)     # stands for p^2 inside logs
LS = sp.Symbol('logS')                   # log p^2, kept as a symbol for fitting
epsd = sp.Symbol('epsd')                 # n = 4 + epsd (dim reg; finite parts need the n-dimensional angular averages)
mu = sp.Symbol('mu', positive=True)      # dim reg scale


def ang_avg_monomial_n(exps, n):
    """n-dimensional angular average of a monomial in the first 4 components of l (only those appear);
    <l_{i1}...l_{i2A}> = (l^2)^A  prod (e_i - 1)!! / prod_{j<A} (n + 2j)."""
    if any(e % 2 for e in exps):
        return sp.Integer(0)
    A = sum(exps)//2
    num = 1
    for e in exps:
        num *= _dfact(e)
    den = sp.Integer(1)
    for j in range(A):
        den *= (n + 2*j)
    return num/den


def dimreg_two_prop_full(power):
    """int d^n l/(2pi)^n (l^2)^a/(l^2+Delta)^2, n = 4 + epsd, MS-bar (mu^(4-n) e^{gamma epsd/2} (4 pi)^(-epsd/2) stripped),
    expanded to O(epsd^0).  Returned in NATURAL units as  A/epsd + B  (A, B carry their own 1/16pi^2).  The n-dimensional
    angular averages are applied by the caller BEFORE the epsd-expansion (bubble_finite)."""
    a = power//2
    n = 4 + epsd
    expr = (mu**(-epsd)*sp.gamma(a + n/2)*sp.gamma(2 - a - n/2)/sp.gamma(n/2)*Dl**(a + n/2 - 2)/(4*sp.pi)**(n/2)
            *(4*sp.pi)**(epsd/2)*sp.exp(-EulerGamma*epsd/2))
    return expr    # kept exact in epsd; expanded after the angular average


lsq = sp.Symbol('lsq', positive=True)     # the n-dimensional l^2 (the loop momentum has components beyond the 4 of p and h)


def _mom(a, b):
    """momentum a*l + b*p, with l the n-dimensional loop momentum and p the 4-dimensional external one"""
    return (sp.sympify(a), sp.sympify(b))


def _dot_n(u, v):
    (a, b), (c, d) = u, v
    return a*c*lsq + (a*d + b*c)*dot(L, P) + b*d*p2


def _comp(u, i):
    a, b = u
    return a*L[i] + b*P[i]


def V1_n(k, q, p_, covariant=False):
    """first-order vertex with n-dimensional dot products; k, q as (a, b) pairs, p_ = (0, +-1)"""
    kq = _dot_n(k, q)
    pp = _dot_n(p_, p_)
    V = sp.Matrix(4, 4, lambda i, j:
        sp.Rational(1, 2)*(_comp(k, i)*_comp(q, j) + _comp(k, j)*_comp(q, i))
        - sp.Rational(1, 2)*delta[i, j]*(kq - m**2)
        + xi*(pp*delta[i, j] - _comp(p_, i)*_comp(p_, j)))
    if covariant:
        corr = sp.Rational(1, 4)*((_dot_n(k, k) + m**2) + (_dot_n(q, q) + m**2))
        V = V - corr*delta
    return V


def ang_avg_n(expr, n):
    """angular average over the n-dimensional l of a polynomial in l1..l4 and lsq: l_{i1}..l_{i2A} -> lsq^A * (n-dim tensor
    average); returns a polynomial in lsq."""
    poly = sp.Poly(sp.expand(expr), *L, lsq)
    res = sp.Integer(0)
    for monom, coeff in poly.terms():
        c = ang_avg_monomial_n(monom[:4], n)
        if c == 0:
            continue
        res += coeff*c*lsq**(sum(monom[:4])//2 + monom[4])
    return sp.expand(res)


def bubble_finite(regulator, covariant=False):
    """Massless bubble, all 256 components, Feynman-parametrised, x-integrated; returns per component the p^4-relevant
    expression with log Lambda^2 (via lgL) or the dim-reg pole (via epsinv = 1/(8 pi^2 epsd)) and log p^2 as the symbol LS.
    regulator: 'cutoff' (Feynman-parametrised, shifted), 'propertime', 'dimreg'.
    n-dimensional dot products throughout (lsq); for 'cutoff' and 'propertime' n = 4 exactly."""
    Delta = x*(1 - x)*Ssym
    k = _mom(1, x); kmp = _mom(1, x - 1)
    mk = _mom(-1, -x); mkmp = _mom(-1, 1 - x)
    pp_ = _mom(0, 1); mp_ = _mom(0, -1)
    Va = V1_n(mk, kmp, pp_, covariant).subs(m, 0)
    Vb = V1_n(mkmp, k, mp_, covariant).subs(m, 0)
    n = 4 + epsd if regulator == 'dimreg' else 4
    if regulator == 'cutoff':
        R = {pw: radial_two_prop(pw) for pw in (0, 2, 4)}
    elif regulator == 'propertime':
        R = {pw: pt_two_prop(pw) for pw in (0, 2, 4)}
    else:
        # expand each radial integral ONCE: R = A/epsd + B + O(epsd)
        R = {}
        for pw in (0, 2, 4):
            ser = sp.series(dimreg_two_prop_full(pw), epsd, 0, 1).removeO()
            ser = sp.expand(sp.expand_log(sp.expand(ser), force=True))
            R[pw] = (sp.expand(ser.coeff(epsd, -1)), sp.expand(ser.coeff(epsd, 0)))
    out = {}
    for (i, j, r, s) in itertools.product(range(4), repeat=4):
        if i > j or r > s:
            continue
        avg = ang_avg_n(Va[i, j]*Vb[r, s], n)
        poly = sp.Poly(avg, lsq)
        val = sp.Integer(0)
        for (pw,), coeff in poly.terms():
            if regulator == 'dimreg':
                # coeff is rational in epsd through the n-dimensional averages: coeff = c0 + c1 epsd + ...
                c0 = coeff.subs(epsd, 0)
                c1 = sp.diff(coeff, epsd).subs(epsd, 0)
                A, B = R[2*pw]
                val += c0*A*8*sp.pi**2*epsinv + sp.expand(c0*B + c1*A)
            else:
                val += coeff*R[2*pw]
        val = sp.expand(-sp.Rational(1, 2)*val.subs(Dl, Delta))
        val = sp.expand(sp.expand_log(val, force=True))
        val = val.subs({sp.log(Ssym): LS, sp.log(x): _lx, sp.log(1 - x): _l1x})
        out[(i, j, r, s)] = _xint(val)
    return out


_lx, _l1x = sp.symbols('logx log1mx')
_xint_cache = {}


def _xint_monomial(k, i, j):
    """int_0^1 x^k log(x)^i log(1-x)^j dx, i, j in {0, 1}"""
    key = (k, i, j)
    if key not in _xint_cache:
        _xint_cache[key] = sp.simplify(sp.integrate(x**k*sp.log(x)**i*sp.log(1 - x)**j, (x, 0, 1)))
    return _xint_cache[key]


def _xint(expr):
    """x-integral over [0, 1] of a polynomial in x, logx, log1mx (and anything else), by the monomial table"""
    poly = sp.Poly(sp.expand(expr), x, _lx, _l1x)
    res = sp.Integer(0)
    for (k, i, j), c in poly.terms():
        res += c*_xint_monomial(k, i, j)
    return sp.expand(res)
