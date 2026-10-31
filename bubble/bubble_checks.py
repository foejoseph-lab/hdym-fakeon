"""bubble_checks.py -- the one-loop bubble (the instrument check; rung 2 at p=0 on the 31 Oct ladder) of a fakeon scalar against the
published numbers, computed two independent ways (bubble_lib.py: method P Feynman parameters with the average
continuation, method M momentum-space residues with the principal value at the pinch), then the Lee-Wick and
principal-value (FW) prescriptions on the same integral.  Untrusted generator; bubble_cert.py recomputes the published
closed forms in ball arithmetic and compares.

Sections (all run by default; each is a separate log when run alone):
    --fakeon   ordinary masses, average continuation: 1801.00915 (5.4), 2001.01942 (6.6), (4.1); D = 2, 3, 4
               method P = method M; (6.3) the fakeon as the E -> 0 limit of a Lee-Wick pair
    --lw       the Lee-Wick pair a = +-i M^2 (1801.00915 eq. 5.5): (5.8) D = 4, the D = 2 closed form, (5.12) D = 3,
               the half-difference f_d, f_+ = f_- below the LW threshold, method P = method M
    --fw       principal-value (FW) mixed term 2001.01942: (2.6) nonlocal divergence, (4.3) pseudothreshold imaginary
               part, (3.1) finite part with the cutoff |k| <= Lambda
    --fig13    1703.04584 propagator (2.1), m = 3, M = 1: Im J(p0) across the LW threshold, Re J = 0, the kink at
               p0 = sqrt(2), the half-difference against the rigid-contour residue term of sec. 3

Runs (Spyder):
    import sys
    %cd C:\\GitHub\\hdym-fakeon\\bubble
    !"{sys.executable}" -B -u bubble_checks.py --fakeon > results\\fakeon.log
    !"{sys.executable}" -B -u bubble_checks.py --lw     > results\\lw.log
    !"{sys.executable}" -B -u bubble_checks.py --fw     > results\\fw.log
    !"{sys.executable}" -B -u bubble_checks.py --fig13  > results\\fig13.log
Outputs: results/<section>.json (every number printed, for bubble_cert.py); the log ends
'# END <section> elapsed=<s>s status=PASS|FAIL' with the worst deviation of the section.
"""
import argparse
import json
import os
import sys
import time

import mpmath as mp

import bubble_lib as L

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

I = L.I
PI = L.PI
mp.mp.dps = 20

# tolerances: method P carries the end-point offset delta = 10^-(dps-6) = 1e-14 and O(sqrt delta) at a threshold;
# quadratures are good to ~1e-18; (2.6) and (3.1) are approached as 1/Lambda^2.
TOL_EXACT = mp.mpf("1e-12")
TOL_LAM = mp.mpf("2e-5")
args_schwarz = False


def say(*a):
    print(*a)
    sys.stdout.flush()


def ns(x, n=15):
    return mp.nstr(x, n)


def cplx(x):
    """A complex number for the json, as two decimal strings at the working precision (floats would lose digits)."""
    x = mp.mpc(x)
    return [mp.nstr(mp.re(x), mp.mp.dps + 2), mp.nstr(mp.im(x), mp.mp.dps + 2)]


def num(x):
    return mp.nstr(mp.mpf(x), mp.mp.dps + 2)


class Section:
    def __init__(self, name):
        self.name = name
        self.t0 = time.time()
        self.worst = mp.mpf(0)
        self.worst_what = ""
        self.fails = []
        self.data = {}
        say(f"# {name}: start, mp.dps = {mp.mp.dps}")

    def check(self, what, dev, tol):
        dev = abs(mp.mpf(dev)) if not isinstance(dev, mp.mpc) else abs(dev)
        ok = dev <= tol
        if dev > self.worst:
            self.worst, self.worst_what = dev, what
        if not ok:
            self.fails.append(what)
        say(f"    [{'ok' if ok else 'FAIL'}] {what}: deviation {ns(dev, 3)} (tol {ns(tol, 2)})")
        return ok

    def end(self):
        el = time.time() - self.t0
        status = "PASS" if not self.fails else "FAIL(" + ",".join(self.fails) + ")"
        self.data["_meta"] = dict(worst=float(self.worst), worst_what=self.worst_what, elapsed=el, status=status,
                                  dps=mp.mp.dps)
        with open(os.path.join(RES, f"{self.name}.json"), "w") as fh:
            json.dump(self.data, fh, indent=1)
        say(f"# worst deviation {ns(self.worst, 3)} at {self.worst_what}")
        say(f"# END {self.name} elapsed={el:.1f}s status={status}")


# ------------------------------------------------------------------------------------------------------------------
def run_fakeon():
    S = Section("fakeon")
    pref = mp.mpc(0, 1)
    tref = pref ** 2
    Ts = [mp.mpf(x) for x in ["0.5", "3", "4.5", "6", "10", "20"]]
    say("## D = 4, ordinary masses, method P f_+-(T) = straight integral at T +- i0, f_AV = (f_+ + f_-)/2;")
    say("##   method M PV / Feynman momentum integrals, UV-subtracted at t_ref = -1 in both methods")
    for (a1, a2, lab) in [(mp.mpf(1), mp.mpf(1), "m1=m2=1"), (mp.mpf(1), mp.mpf("2.25"), "m1=1,m2=1.5")]:
        say(f"### {lab}  (threshold (m1+m2)^2 = {ns((mp.sqrt(a1)+mp.sqrt(a2))**2, 6)})")
        rows = []
        for T in Ts:
            p0 = mp.sqrt(T)
            fp = L.I_pm_real_masses(4, T, a1, a2, +1)
            fm = L.I_pm_real_masses(4, T, a1, a2, -1)
            fav = (fp + fm) / 2
            Iref = L.I_straight(4, tref, a1, a2)
            Bpar_pv = -(I / (16 * PI ** 2)) * (fav - Iref)
            Bpar_F = -(I / (16 * PI ** 2)) * (fp - Iref)
            Bmom_pv = L.B_momentum(4, p0, a1, a2, "PV", p0_ref=pref)
            Bmom_F = L.B_momentum(4, p0, a1, a2, "feynman", p0_ref=pref)
            row = dict(T=float(T), f_plus=cplx(fp), f_minus=cplx(fm), f_AV=cplx(fav), I_ref=cplx(Iref),
                       B_PV_par=cplx(Bpar_pv), B_PV_mom=cplx(Bmom_pv), B_F_par=cplx(Bpar_F), B_F_mom=cplx(Bmom_F))
            say(f"  T={ns(T,4)}: f_AV = {ns(fav)}   f_+ - f_- = {ns(fp - fm, 10)}")
            S.check(f"{lab} T={T} D=4 Im f_AV = 0", mp.im(fav), TOL_EXACT)
            S.check(f"{lab} T={T} D=4 PV method P = method M", Bpar_pv - Bmom_pv, TOL_EXACT)
            S.check(f"{lab} T={T} D=4 Feynman method P = method M", Bpar_F - Bmom_F, TOL_EXACT)
            if a1 == a2 == 1:
                # 1801.00915 eq. 5.4: f_AV = (1/2) int ln[(1 - z x(1-x))^2], z = p^2/m^2 (subdivided at the roots)
                pts = [mp.mpf(0)]
                if T > 4:
                    pts += [(1 - mp.sqrt(1 - 4 / T)) / 2, mp.mpf(1) / 2, (1 + mp.sqrt(1 - 4 / T)) / 2]
                else:
                    pts += [mp.mpf(1) / 2]
                pts.append(mp.mpf(1))
                c54 = mp.quad(lambda x: mp.log((1 - T * x * (1 - x)) ** 2), pts) / 2
                row["eq54"] = num(c54)
                S.check(f"{lab} T={T} eq. 5.4", fav - c54, TOL_EXACT)
            # 2001.01942 eq. 6.6 and 4.1 (hold for the UV-subtracted B too: the subtraction is real-analytic)
            d66 = L.fakeon_D4_diff_closed(T, a1, a2)
            im41 = L.feynman_D4_im_closed(T, a1, a2)
            row["eq66"] = num(d66)
            row["eq41"] = num(im41)
            say(f"           B_PV - B_F = {ns(Bmom_pv - Bmom_F, 12)}  eq. 6.6 = {ns(d66, 12)};   2 Im(-i B_F) = "
                f"{ns(2 * mp.im(-I * Bmom_F), 12)}  eq. 4.1 = {ns(im41, 12)}")
            S.check(f"{lab} T={T} eq. 6.6", (Bmom_pv - Bmom_F) - d66, TOL_EXACT)
            S.check(f"{lab} T={T} eq. 4.1", 2 * mp.im(-I * Bmom_F) - im41, TOL_EXACT)
            # the discontinuity of ln Delta: f_+ - f_- = -2 pi i (x2 - x1) above threshold
            if T > (mp.sqrt(a1) + mp.sqrt(a2)) ** 2:
                x1, x2 = sorted([mp.re(r) for r in L.delta_roots(T, a1, a2)])
                S.check(f"{lab} T={T} f_+ - f_- = -2 pi i (x2-x1)", (fp - fm) + 2 * PI * I * (x2 - x1), TOL_EXACT)
            rows.append(row)
        S.data[f"D4_{lab}"] = rows
    say("## D = 3 and D = 2, ordinary masses: method P average = method M PV, and the Feynman prescriptions")
    for D in (3, 2):
        pref_D = {2: I / (4 * PI), 3: I / (8 * PI)}[D]
        rows = []
        for (a1, a2, lab) in [(mp.mpf(1), mp.mpf(1), "m1=m2=1"), (mp.mpf(1), mp.mpf("2.25"), "m1=1,m2=1.5")]:
            for T in [mp.mpf(3), mp.mpf(6), mp.mpf(20)]:
                p0 = mp.sqrt(T)
                fp = L.I_pm_real_masses(D, T, a1, a2, +1)
                fm = L.I_pm_real_masses(D, T, a1, a2, -1)
                Bpar_pv = pref_D * (fp + fm) / 2
                Bpar_F = pref_D * fp
                Bmom_pv = L.B_momentum(D, p0, a1, a2, "PV")
                Bmom_F = L.B_momentum(D, p0, a1, a2, "feynman")
                say(f"  D={D} {lab} T={ns(T,3)}: B_PV = {ns(Bmom_pv)}   B_F = {ns(Bmom_F)}")
                S.check(f"{lab} T={T} D={D} PV P=M", Bpar_pv - Bmom_pv, TOL_EXACT)
                S.check(f"{lab} T={T} D={D} Feynman P=M", Bpar_F - Bmom_F, TOL_EXACT)
                rows.append(dict(D=D, lab=lab, T=float(T), B_PV_par=cplx(Bpar_pv), B_PV_mom=cplx(Bmom_pv),
                                 B_F_par=cplx(Bpar_F), B_F_mom=cplx(Bmom_F)))
        S.data[f"D{D}"] = rows
    say("## 1801.00915 eq. 6.3: the fakeon propagator as the E -> 0 limit of (1/2)[1/(p^2-m^2+iE^2) + 1/(p^2-m^2-iE^2)]")
    say("##   J_E(T=6) - B_PV(T=6), m = 1, both by method M; expected -> 0 (as E^2 in D = 4, E^4 in D = 3)")
    T = mp.mpf(6)
    p0 = mp.sqrt(T)
    rows = []
    for D in (3, 4):
        ref = L.B_momentum(D, p0, 1, 1, "PV", p0_ref=(pref if D == 4 else None))
        prev = None
        for E in [mp.mpf(x) for x in ["0.2", "0.1", "0.05", "0.025"]]:
            ap, am = 1 + I * E * E, 1 - I * E * E
            J = sum(L.B_momentum(D, p0, a, b, "PV", p0_ref=(pref if D == 4 else None))
                    for a in (ap, am) for b in (ap, am)) / 4
            dv = J - ref
            ratio = (abs(prev) / abs(dv)) if prev is not None else None
            say(f"  D={D} E={ns(E,4)}: J_E - B_PV = {ns(dv, 5)}" + (f"   ratio to previous E: {ns(ratio, 4)}" if ratio else ""))
            rows.append(dict(D=D, E=float(E), J_E=cplx(J), B_PV=cplx(ref), diff=cplx(dv)))
            prev = dv
        # convergence check: the last step must shrink by at least 3 (E^2 gives 4, E^4 gives 16)
        S.check(f"D={D} eq. 6.3 E->0 convergence ratio >= 3", 0 if ratio >= 3 else ratio, mp.mpf(0))
    S.data["eq63"] = rows
    S.end()


# ------------------------------------------------------------------------------------------------------------------
def run_lw():
    S = Section("lw")
    Ts = [mp.mpf(x) for x in ["0.5", "1", "1.5", "1.9", "2.1", "3", "5", "8"]]
    pref = mp.mpc(0, 1)
    tref = pref ** 2
    say("## Lee-Wick pair a = +-i (M = 1), 1801.00915 eq. 5.5; f = (+-) + (-+) terms, g = (++) + (--) terms")
    say("## method P: continuation around t = 0 on |t| = 1, along the real axis, a semicircle (radius min(0.1, half the distance to the nearest requested T) = 0.05 here) around the LW threshold t = 2; monodromy tracked by continuity")
    for D in (4, 2):
        say(f"### D = {D}")
        res = {}
        for a1, a2, nm in [(I, -I, "pm"), (-I, I, "mp"), (I, I, "pp"), (-I, -I, "mm")]:
            res[nm] = (L.continuation_sweep(D, Ts, a1, a2, True, log=say), L.continuation_sweep(D, Ts, a1, a2, False, log=say))
        pref_D = {2: I / (4 * PI), 4: -I / (16 * PI ** 2)}[D]
        if D == 4:
            Iref = sum(L.I_straight(4, tref, a, b) for a in (I, -I) for b in (I, -I))
        rows = []
        for T in Ts:
            p0 = mp.sqrt(T)
            fp = res["pm"][0][T][0] + res["mp"][0][T][0]
            fm = res["pm"][1][T][0] + res["mp"][1][T][0]
            gp = res["pp"][0][T][0] + res["mm"][0][T][0]
            gm = res["pp"][1][T][0] + res["mm"][1][T][0]
            fav, fd, gav, gd = (fp + fm) / 2, (fp - fm) / 2, (gp + gm) / 2, (gp - gm) / 2
            worst_ratio = max(res[nm][s][T][1]["max_ratio"] for nm in res for s in (0, 1))
            if D == 4:
                f0, g, f, fav_pub = L.lw_closed_D4(T)
                fd_pub = (2 * PI / T) * mp.sqrt(T * T - 4) if T > 2 else mp.mpf(0)
                Jpar = pref_D * ((fav + gav) - Iref) / 4
                Jmom = sum(L.B_momentum(4, p0, a, b, "PV", p0_ref=pref) for a in (I, -I) for b in (I, -I)) / 4
                lab = "eq. 5.8"
            else:
                f0, g, f, fav_pub = L.lw_closed_D2(T)
                # the published 4 pi/sqrt(4 - t^2) continues to +-4 pi i/sqrt(t^2 - 4) above the threshold: f_AV = f0 there
                fd_pub = 4 * PI / mp.sqrt(T * T - 4) if T > 2 else mp.mpf(0)
                Jpar = pref_D * (fav + gav) / 4
                Jmom = sum(L.B_momentum(2, p0, a, b, "PV") for a in (I, -I) for b in (I, -I)) / 4
                lab = "sec. 5 D=2"
            say(f"  T={ns(T,4)}: f_AV = {ns(fav)}  published {ns(fav_pub)} | g_AV = {ns(gav)}  g = {ns(g)} | "
                f"f_d = {ns(fd, 8)} (|continued published term| {ns(fd_pub, 8)}) | g_d = {ns(gd, 3)} | tracking worst ratio {ns(worst_ratio, 2)}")
            S.check(f"D={D} T={T} f_AV vs {lab}", fav - fav_pub, TOL_EXACT * (10 if abs(T - 2) < mp.mpf("0.2") else 1))
            S.check(f"D={D} T={T} g_AV vs g", gav - g, TOL_EXACT)
            S.check(f"D={D} T={T} |f_d| vs continued published term", abs(fd) - fd_pub, TOL_EXACT * (10 if abs(T - 2) < mp.mpf("0.2") else 1))
            S.check(f"D={D} T={T} g_d = 0", gd, TOL_EXACT)
            S.check(f"D={D} T={T} Re f_d = 0, Im f_AV = 0", abs(mp.re(fd)) + abs(mp.im(fav)), TOL_EXACT)
            S.check(f"D={D} T={T} J_LW method P = method M", Jpar - Jmom, TOL_EXACT)
            S.check(f"D={D} T={T} tracking unambiguous (ratio < 0.5)", 0 if worst_ratio < mp.mpf("0.5") else worst_ratio, mp.mpf(0))
            rows.append(dict(T=float(T), f_plus=cplx(fp), f_minus=cplx(fm), f_AV=cplx(fav), f_d=cplx(fd), g_AV=cplx(gav),
                             g_d=cplx(gd), f0_pub=num(f0), g_pub=num(g), f_AV_pub=num(fav_pub), fd_pub=num(fd_pub),
                             J_par=cplx(Jpar), J_mom=cplx(Jmom), worst_ratio=float(worst_ratio)))
        if D == 4:
            S.data["D4_Iref"] = cplx(Iref)
        S.data[f"D{D}"] = rows
    say("### D = 3: method M PV against the continued closed forms of eq. 5.12 (Re[f0 + g] at T +- i0)")
    rows = []
    # Euclidean anchor first
    t = mp.mpf("-1.69")
    f0c, gc = L.lw_closed_D3_euclid(t)
    Jm = sum(L.B_momentum(3, mp.mpc(0, mp.sqrt(-t)), a, b) for a in (I, -I) for b in (I, -I)) / 4
    say(f"  Euclidean t={t}: closed f0+g = {ns(f0c + gc)}   method M = {ns(Jm * 32 * PI / I)}")
    S.check("D=3 Euclidean t=-1.69 eq. 5.12 closed = method M", Jm * 32 * PI / I - (f0c + gc), TOL_EXACT)
    S.data["D3_euclid"] = dict(t=float(t), f0=cplx(f0c), g=cplx(gc), J_mom=cplx(Jm))

    def closed_D3(tc):
        tc = mp.mpc(tc)
        r = mp.sqrt(-tc)
        f0 = (2 * I / r) * mp.log((mp.sqrt(2) - I * r) / (mp.sqrt(2) + I * r))
        si, smi = mp.sqrt(I), mp.sqrt(-I)
        g = (I / r) * (mp.log((2 * si + r) / (2 * si - r)) + mp.log((2 * smi - r) / (2 * smi + r)))
        return f0, g

    dl = mp.mpf(10) ** (-(mp.mp.dps - 6))
    for T in Ts:
        p0 = mp.sqrt(T)
        Jm = sum(L.B_momentum(3, p0, a, b, "PV") for a in (I, -I) for b in (I, -I)) / 4
        fg_mom = Jm * 32 * PI / I
        fp0, gp0 = closed_D3(mp.mpc(T, dl))
        fm0, gm0 = closed_D3(mp.mpc(T, -dl))
        re_p, re_m = mp.re(fp0 + gp0), mp.re(fm0 + gm0)
        say(f"  T={ns(T,4)}: method M [f0+g] = {ns(fg_mom)}   Re closed(T+i0) = {ns(re_p)}   Re closed(T-i0) = {ns(re_m)}"
            f"   Im f0(T+i0) = {ns(mp.im(fp0), 8)} (2 pi/sqrt T = {ns(2 * PI / mp.sqrt(T), 8)} above threshold)")
        S.check(f"D=3 T={T} eq. 5.12 Re closed = method M", fg_mom - re_p, TOL_EXACT * 10)
        if args_schwarz:
            S.check(f"D=3 T={T} Schwarz reflection of the closed form (consistency of the transcription only)", L.schwarz_gap(re_p, re_m), TOL_EXACT)
        rows.append(dict(T=float(T), fg_mom=cplx(fg_mom), f0_plus=cplx(fp0), g_plus=cplx(gp0), J_mom=cplx(Jm)))
    S.data["D3"] = rows
    S.end()


# ------------------------------------------------------------------------------------------------------------------
def run_fw():
    S = Section("fw")
    say("## principal-value (FW) mixed term Sigma' (one +i eps, one -i eps propagator), D = 4, cutoff |k| <= Lambda")
    a1, a2 = mp.mpf(1), mp.mpf("0.25")
    say(f"### m1 = 1, m2 = 0.5: pseudothreshold (m1-m2)^2 = 0.25, threshold 2.25")
    rows = []
    Lams = [mp.mpf(50), mp.mpf(500), mp.mpf(5000)]
    for T in [mp.mpf("0.1"), mp.mpf(1), mp.mpf(5)]:
        p0 = mp.sqrt(T)
        vals = {}
        for Lam in Lams:
            pv, dl = L.Sigma_FW_mixed(4, p0, a1, a2, Lam)
            vals[Lam] = pv + dl
        slope0 = (vals[Lams[1]] - vals[Lams[0]]) / mp.log(100)
        slope = (vals[Lams[2]] - vals[Lams[1]]) / mp.log(100)
        c26 = I * (a1 - a2) / ((4 * PI) ** 2 * T)
        im = 2 * mp.im(-I * vals[Lams[2]])
        c43 = L.fw_D4_im_closed(T, a1, a2)
        say(f"  T={ns(T,3)}: d Sigma'/d ln Lambda^2 (Lambda 50->500) = {ns(slope0, 10)}, (500->5000) = {ns(slope, 10)}   eq. 2.6 = {ns(c26, 10)}"
            f"   deviation ratio {ns(abs(slope0 - c26) / abs(slope - c26), 5)} (100 for 1/Lambda^2)"
            f"   2 Im(-i Sigma') = {ns(im, 12)}   eq. 4.3 = {ns(c43, 12)}")
        S.check(f"T={T} eq. 2.6 nonlocal divergence coefficient", slope - c26, TOL_LAM * abs(c26))
        S.check(f"T={T} eq. 2.6 approached as 1/Lambda^2 (deviation ratio in [90, 110])",
                0 if 90 < abs(slope0 - c26) / abs(slope - c26) < 110 else 1, mp.mpf(0))
        S.check(f"T={T} eq. 4.3 pseudothreshold imaginary part", im - c43, TOL_EXACT)
        rows.append(dict(T=float(T), Sigma=[[float(Lm), cplx(v)] for Lm, v in vals.items()], slope=cplx(slope), slope0=cplx(slope0), eq26=cplx(c26),
                         im=num(im), eq43=num(c43)))
    S.data["m1_m2"] = rows
    say("### m2 = 0, m1 = m = 1: eq. 3.1  i(4pi)^2 Sigma' = -(m^2/p^2) ln(4 Lambda^2/m^2) - ((p^2-m^2)/p^2) ln((p^2-m^2+i eps)/p^2)")
    a1, a2 = mp.mpf(1), mp.mpf(0)
    rows = []
    for T in [mp.mpf("0.5"), mp.mpf(3)]:
        p0 = mp.sqrt(T)
        for Lam in [mp.mpf(200), mp.mpf(2000)]:
            pv, dl = L.Sigma_FW_mixed(4, p0, a1, a2, Lam)
            lhs = I * (4 * PI) ** 2 * (pv + dl)
            rhs = -(1 / T) * mp.log(4 * Lam ** 2) - ((T - 1) / T) * mp.log(mp.mpc(T - 1, mp.mpf("1e-30")) / T)
            say(f"  T={ns(T,3)} Lambda={ns(Lam,4)}: i(4pi)^2 Sigma' = {ns(lhs, 12)}   eq. 3.1 = {ns(rhs, 12)}   diff = {ns(lhs - rhs, 4)}")
            S.check(f"T={T} Lambda={Lam} eq. 3.1 (O(1/Lambda^2))", lhs - rhs, mp.mpf(10) / Lam ** 2)
            rows.append(dict(T=float(T), Lam=float(Lam), lhs=cplx(lhs), eq31=cplx(rhs)))
    S.data["m2_0"] = rows
    say("### D = 3 and D = 2 (2001.01942 sec. 7), m1 = 1, m2 = 0.5: 2 Im(-i Sigma') = -theta(u-)/(4 sqrt p^2) (D=3), -theta(u-)/sqrt(u+ u-) (D=2)")
    a1, a2 = mp.mpf(1), mp.mpf("0.25")
    rows = []
    for D in (3, 2):
        for T in [mp.mpf("0.1"), mp.mpf("0.2"), mp.mpf(1)]:
            p0 = mp.sqrt(T)
            pv, dl = L.Sigma_FW_mixed(D, p0, a1, a2, mp.mpf(10) ** 6)
            im = 2 * mp.im(-I * (pv + dl))
            up, um = (mp.sqrt(a1) + mp.sqrt(a2)) ** 2 - T, (mp.sqrt(a1) - mp.sqrt(a2)) ** 2 - T
            pub = (-(1 / (4 * mp.sqrt(T))) if D == 3 else -(1 / mp.sqrt(up * um))) if um > 0 else mp.mpf(0)
            say(f"  D={D} T={ns(T,3)}: 2 Im(-i Sigma') = {ns(im, 12)}   published {ns(pub, 12)}")
            S.check(f"D={D} T={T} sec. 7 imaginary part", im - pub, TOL_EXACT)
            rows.append(dict(D=D, T=float(T), im=num(im), pub=num(pub)))
    S.data["D3_D2"] = rows
    S.end()


# ------------------------------------------------------------------------------------------------------------------
def run_fig13():
    S = Section("fig13")
    say("## 1703.04584 eq. 2.1: D(p) = M^4/((p^2 - m^2 + i eps)((p^2)^2 + M^4)), m = 3, M = 1, D = 4 (their fig. 13 setting)")
    say("##   J = sum_ab c_a c_b B_ab, poles a = m^2, +i, -i with partial-fraction weights c; LW pinch averaged (PV), physical pinch Feynman")
    m2 = mp.mpf(9)
    poles = [m2, I, -I]
    c = {m2: 1 / (m2 ** 2 + 1), I: 1 / ((I - m2) * (2 * I)), -I: 1 / ((-I - m2) * (-2 * I))}
    say(f"   weights: c_m = {ns(c[m2], 8)}, c_+ = {ns(c[I], 8)}, c_- = {ns(c[-I], 8)}  (sum = {ns(sum(c.values()), 4)}, "
        f"sum c_a a = {ns(sum(c[a] * a for a in poles), 4)}: the propagator falls as 1/p^6)")
    pref = mp.mpc(0, 1)

    def J(p0, presc_lw="PV"):
        tot = mp.mpc(0)
        for a in poles:
            for b in poles:
                presc = "feynman" if (a == m2 and b == m2) else presc_lw
                tot += c[a] * c[b] * L.B_momentum(4, p0, a, b, presc, p0_ref=pref)
        return tot

    grid = [mp.mpf(x) / 10 for x in range(5, 26)]
    rows = []
    say("   p0      Im J(p0) - Im J(i)         Re J   (Re J = 0 below p0 = 2m is a Schwarz-reflection identity of the averaged terms; a consistency check of the code, not of the prescription)")
    for p0 in grid:
        v = J(p0)
        say(f"   {ns(p0,3):6s}  {ns(mp.im(v), 14):24s}  {ns(mp.re(v), 3)}")
        S.check(f"p0={p0} Re J = 0", mp.re(v), TOL_EXACT)
        rows.append(dict(p0=float(p0), J=cplx(v)))
    S.data["grid"] = rows
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fine = [mp.mpf(5) / 10 + mp.mpf(k) * (mp.mpf(2)) / 160 for k in range(161)]
        ys = [float(mp.im(J(p0)) - mp.im(J(I))) for p0 in fine]
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot([float(p) for p in fine], ys, lw=1.5)
        ax.axvline(float(mp.sqrt(2)), ls="--", lw=0.8, color="gray")
        ax.set_xlabel("p0"); ax.set_ylabel("Im J(p0) - Im J(i)")
        ax.set_title("1703.04584 eq. 2.1, m=3, M=1 (their fig. 13): kink at p0 = sqrt 2")
        fig.tight_layout(); fig.savefig(os.path.join(RES, "fig13.png"), dpi=150); plt.close(fig)
        say("   plot saved: results/fig13.png (compare by eye with 1703.04584 fig. 13)")
    except Exception as e:  # matplotlib absent: the numbers above are the record
        say(f"   plot not saved ({e})")
    r2 = mp.sqrt(2)
    hs = [mp.mpf("0.02"), mp.mpf("0.01"), mp.mpf("0.005")]
    say("   one-sided derivatives of Im J at the LW threshold p0 = sqrt 2 (p^2 = 2 M^2):")
    ders = []
    for h in hs:
        dl = (mp.im(J(r2)) - mp.im(J(r2 - h))) / h
        dr = (mp.im(J(r2 + h)) - mp.im(J(r2))) / h
        say(f"     h={ns(h,3)}: left {ns(dl, 8)}   right {ns(dr, 8)}   jump {ns(dr - dl, 8)}")
        ders.append(dict(h=float(h), left=float(dl), right=float(dr)))
    S.data["kink"] = ders
    jump = ders[-1]["right"] - ders[-1]["left"]
    S.check("derivative of Im J jumps at p0 = sqrt 2 (|jump| > 1e-3)", 0 if abs(jump) > mp.mpf("1e-3") else 1, mp.mpf(0))
    # a smooth point for contrast
    p0s = mp.mpf("1.9")
    dl = (mp.im(J(p0s)) - mp.im(J(p0s - hs[-1]))) / hs[-1]
    dr = (mp.im(J(p0s + hs[-1])) - mp.im(J(p0s))) / hs[-1]
    say(f"     contrast at p0 = 1.9: left {ns(dl, 8)} right {ns(dr, 8)} jump {ns(dr - dl, 3)}")
    S.check("Im J smooth away from the threshold (|jump| at p0 = 1.9 < 1e-5)", 0 if abs(dr - dl) < mp.mpf("1e-5") else abs(dr - dl), mp.mpf(0))
    S.data["smooth"] = dict(p0=float(p0s), left=float(dl), right=float(dr))
    say("   half-difference of the two continuations above the LW threshold against the rigid-contour residue of")
    say("   1703.04584 sec. 3, (1/16 pi) M^4/(m^4 + M^4) sqrt(1 - 4 M^4/(p^2)^2)  [given the D = 4 LW f_d above, the ratio 1/2 is")
    say("   an algebraic identity between the two readings: J+ - J- equals that residue; the independent content is f_d itself]:")
    rows = []
    for p0 in [mp.mpf("1.6"), mp.mpf(2), mp.mpf("2.5")]:
        t = p0 ** 2
        jd = (J(p0, "feynman") - J(p0, "minus")) / 2
        pub = (1 / (16 * PI)) * (1 / (m2 ** 2 + 1)) * mp.sqrt(1 - 4 / t ** 2)
        say(f"     p0={ns(p0,3)}: J_d = {ns(jd, 10)}   published magnitude {ns(pub, 10)}   ratio |J_d|/published = {ns(abs(jd) / pub, 8)}")
        rows.append(dict(p0=float(p0), J_d=cplx(jd), pub=float(pub)))
    S.data["halfdiff"] = rows
    S.check("J_d is real (the LW discontinuity is on the real part of J)", mp.im(jd), TOL_EXACT)
    S.end()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    for s in ("fakeon", "lw", "fw", "fig13"):
        ap.add_argument(f"--{s}", action="store_true")
    ap.add_argument("--dps", type=int, default=20)
    ap.add_argument("--schwarz", action="store_true", help="also run the Schwarz-reflection consistency check (not standing: a tool to call when something looks odd)")
    args = ap.parse_args()
    mp.mp.dps = args.dps
    args_schwarz = args.schwarz
    sel = [s for s in ("fakeon", "lw", "fw", "fig13") if getattr(args, s)] or ["fakeon", "lw", "fw", "fig13"]
    for s in sel:
        globals()[f"run_{s}"]()
