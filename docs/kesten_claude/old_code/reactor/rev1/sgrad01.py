import numpy as np

def dp_func(x, y, z):
    return 14.7 * y / z * (x / 492.0) ** 1.823 * (1.0 - np.exp(-0.0672 * z * 492.0 / (14.7 * x)))

def kcf(a, b, c, d, e):
    return 0.61 * a / b * (c / b * d) ** -0.667 * (a / (e * c)) ** -0.41

def eval1(a, b):
    return b**3 / 3.0 - a**3 / 3.0

def eval2(a, b):
    return b**2 / 2.0 - a**2 / 2.0

def unbar(table, T):
    return np.interp(T, np.linspace(300, 2000, len(table)), table)

def dp3f(temp, diff, press):
    return 14.7 * diff / press * (temp / 492.0) ** 1.823 * (
        1.0 - np.exp(-0.0672 * press * 492.0 / (14.7 * temp))
    )

def update_tps(t, h4, kc4, ci4, h3, dp3, dcpdx, hc):
    tps = t - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / hc
    return max(tps, 1.0)

def iterate_tps_ci3_convergence(t, ci3, kc3, dp3, kc4, ci4, h4tbl, h3tbl, d03, p, hc, a, tol=0.05):
    cps = ci3 / 2
    dcpdx = kc3 / dp3 * (ci3 - cps)
    h4 = unbar(h4tbl, t)
    h3 = unbar(h3tbl, t)
    tps = update_tps(t, h4, kc4, ci4, h3, dp3, dcpdx, hc)
    tmtp0 = t - tps

    for _ in range(50):
        h3 = unbar(h3tbl, tps)
        dp3 = dp3f(tps, d03, p)
        dcpdx = kc3 / dp3 * (ci3 - cps)
        tps_new = update_tps(t, h4, kc4, ci4, h3, dp3, dcpdx, hc)
        tmtp1 = t - tps_new
        if abs(tmtp0 - tmtp1) / tmtp1 < tol:
            break
        tps = 0.5 * (tps + tps_new)
        tmtp0 = tmtp1

    return cps, dcpdx, dp3, h3, tps

def solve_cp_profile(ci3, kc3, dp3, a, k0, gamma, beta, en3, max_iter=50, nx=24):
    nx1 = nx + 1
    xoa = 0.0
    vnu = -kc3 / dp3
    delxoa = (1.0 - xoa) / nx

    cpx = np.zeros(nx1 + 1)
    cpox = np.zeros(nx1 + 1)
    pcpox = np.zeros(nx1 + 1)
    rhet = np.zeros(nx1 + 1)
    dx = np.zeros(nx1 + 1)

    lp2 = 1
    for _ in range(max_iter):
        xa = xoa
        for i in range(nx1):
            cpx[i] = (xa - xoa) / (1.0 - xoa) * ci3
            denom = 1.0 + beta * (1.0 - cpx[i] / ci3)
            rhet[i] = k0 * ci3 ** (1.0 - en3) * cpx[i] ** en3 * np.exp(gamma * beta * (1.0 - cpx[i] / ci3) / denom)
            dx[i] = xa
            xa += delxoa

        for i in range(nx):
            cpx[i] = (cpx[i] + cpx[i + 1]) / 2.0
            rhet[i] = (rhet[i] + rhet[i + 1]) / 2.0

        dxl = xoa
        dxu = dxl + delxoa
        rr1 = 0.0
        ctrm = (a * vnu + 1.0) / (a * vnu)
        for i in range(nx):
            rr1 += rhet[i] * (eval2(dxl, dxu) - ctrm * eval1(dxl, dxu))
            dxl = dxu
            dxu += delxoa

        cpox[0] = max(ci3 - a**2 / (dp3 * rr1), 0.0)

        int1 = 1
        r1, r2, ps1, ps2 = 0.0, 0.0, 0.0, 0.0
        xoa = 0.0
        xa = xoa + delxoa
        for k in range(1, nx):
            for i in range(int1):
                r1 += rhet[i] * eval1(xoa, xa)
                xoa = xa
                xa += delxoa
            r1 *= (1.0 / xoa - ctrm)

            xad = xa
            xa -= delxoa
            for i in range(int1, nx - 1):
                ps1 += rhet[i + 1] * eval2(xa, xad)
                ps2 += rhet[i + 1] * eval1(xa, xad)
                xa = xad
                xad += delxoa
            r2 = ps1 - ctrm * ps2

            cpox[k] = max(ci3 - a**2 / dp3 * (r1 + r2), 0.0)
            xoa = 0.0
            xa = xoa + delxoa
            int1 += 1
            r1 = r2 = ps1 = ps2 = 0.0

        dxl = xoa
        dxu = dxl + delxoa
        rr2 = 0.0
        for i in range(nx):
            rr2 += rhet[i] * eval1(dxl, dxu)
            dxl = dxu
            dxu += delxoa
        cpox[nx1] = max(ci3 - a**2 / dp3 * (1.0 - ctrm) * rr2, 0.0)

        return cpox[: nx1 + 1], cpx[: nx1 + 1], rhet[: nx1 + 1]

    return cpox, cpx, rhet

def sgrad(temp_vars, tables, constants, shared, results):
    T = temp_vars["TEMP"]
    P = constants["PRES"]
    CI1, CI2, CI3, CI4 = temp_vars["C1"], temp_vars["C2"], temp_vars["C3"], temp_vars["C4"]
    ALPHA2 = constants["ALPHA2"]
    D03, D04 = constants["DIF3"], constants["DIF4"]
    G = temp_vars["G"]
    AP = temp_vars["AP"]
    A = shared["A"]
    KP = constants["KP"]
    BGM = constants["BGM"]
    EN3 = constants["EN3"]

    RHO = CI1 + CI2 + CI3 + CI4

    DI3 = D03 * 14.7 / P * (T / 492.0) ** 1.823
    DI4 = D04 * 14.7 / P * (T / 492.0) ** 1.823

    MU = unbar(tables["VISVST"], T)
    CF4 = unbar(tables["CFTBL4"], T)
    CF3 = unbar(tables["CFTBL3"], T)
    CF2 = unbar(tables["CFTBL2"], T)
    CF1 = unbar(tables["CFTBL1"], T)
    H4 = unbar(tables["H4TBL"], T)
    H3 = unbar(tables["H3TBL"], T)

    KC3 = kcf(G, RHO, MU, DI3, AP)
    KC4 = kcf(G, RHO, MU, DI4, AP)

    CFBAR = (CI1 * CF1 + CI2 * CF2 + CI3 * CF3 + CI4 * CF4) / RHO
    HC = 0.74 * G * CFBAR * (G / (AP * MU)) ** -0.41

    DP3 = dp3f(T, D03, P)

    CPS, DCPDX, DP3, H3, TPS = iterate_tps_ci3_convergence(
        T, CI3, KC3, DP3, KC4, CI4, tables["H4TBL"], tables["H3TBL"], D03, P, HC, A
    )

    GAMMA = BGM / TPS
    BETA = -CPS * H3 * DP3 / (KP * TPS)
    K0 = ALPHA2 * np.exp(-GAMMA) * CI1 ** EN3

    CPOX, CPX, RHET = solve_cp_profile(CI3, KC3, DP3, A, K0, GAMMA, BETA, EN3)

    DCPDX = KC3 / DP3 * (CI3 - CPOX[-1])
    GRAD = DCPDX * DP3
    TGRAD = HC * (T - TPS)

    results["GRAD"] = GRAD
    results["TGRAD"] = TGRAD
    results["CPOX"] = CPOX
    results["CPX"] = CPX
    return GRAD, TGRAD
