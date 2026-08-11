import numpy as np


def _dp3f(temp, d0, pres):
    temp = max(float(temp), 1e-12)
    pres = max(float(pres), 1e-12)
    return 14.7 * d0 / pres * (temp / 492.0) ** 1.823 * (
        1.0 - np.exp(-0.0672 * pres * 492.0 / (14.7 * temp))
    )


def _kcf(g, rho, mu, diff, ap):
    rho = max(float(rho), 1e-30)
    mu = max(float(mu), 1e-30)
    diff = max(float(diff), 1e-30)
    ap = max(float(ap), 1e-30)
    return 0.61 * g / rho * (mu / (rho * diff)) ** (-0.667) * (g / (ap * mu)) ** (-0.41)


def _eval1(a, b):
    return b ** 3 / 3.0 - a ** 3 / 3.0


def _eval2(a, b):
    return b ** 2 / 2.0 - a ** 2 / 2.0


def _rhetf(ci3, cpx, gamma, beta, k0, order_n):
    if ci3 <= 0.0 or cpx <= 0.0 or k0 <= 0.0:
        return 0.0
    ratio = 1.0 - cpx / ci3
    denom = 1.0 + beta * ratio
    if abs(denom) < 1e-30:
        return 0.0
    return k0 * ci3 ** (1.0 - order_n) * cpx ** order_n * np.exp(gamma * beta * ratio / denom)


def _cpxf(x, x0a, cps):
    if x <= x0a:
        return 0.0
    denom = 1.0 - x0a
    if abs(denom) < 1e-30:
        return cps
    return (x - x0a) / denom * cps


def _trapp(u, v, npart, x0a, cps, ci3, gamma, beta, k0, order_n):
    n = max(int(npart) - 1, 1)
    part = float(npart)
    h = (v - u) / part
    uph = u + h
    cpx1 = _cpxf(u, x0a, cps)
    cpx2 = _cpxf(v, x0a, cps)
    rhet1 = _rhetf(ci3, cpx1, gamma, beta, k0, order_n)
    rhet2 = _rhetf(ci3, cpx2, gamma, beta, k0, order_n)
    trm1 = (u ** 2 * rhet1) / 2.0
    trm2 = (v ** 2 * rhet2) / 2.0
    total = 0.0
    for _ in range(n):
        cpx = _cpxf(uph, x0a, cps)
        rhet = _rhetf(ci3, cpx, gamma, beta, k0, order_n)
        total += uph ** 2 * rhet
        uph += h
    return h * (trm1 + total + trm2)


def SGRAD(TEMP, G, concentrations, A, AP, target_reaction, target_species, cfg, props, pres=None, verbose=False):
    """
    Closer port of Kesten's vapor SGRAD routine.
    Returns ammonia mass gradient (GRAD) and heat-transfer term (TGRAD).
    """
    if target_species is None or target_species.name != "NH3":
        return 0.0, 0.0, {"note": "SGRAD only ported for NH3 path"}

    temp = float(TEMP)
    p = float(cfg.PRES if pres is None else pres)
    if temp <= 0.0 or p <= 0.0 or A <= 0.0 or AP <= 0.0 or G <= 0.0:
        return 0.0, 0.0, {"note": "degenerate conditions"}

    ci1 = max(float(concentrations.get("H2", 0.0)), 0.0)
    ci2 = max(float(concentrations.get("N2", 0.0)), 0.0)
    ci3 = max(float(concentrations.get("NH3", 0.0)), 0.0)
    ci4 = max(float(concentrations.get("N2H4", 0.0)), 0.0)
    if ci3 <= 0.0:
        return 0.0, 0.0, {"note": "no NH3 present"}

    rho = max(ci1 + ci2 + ci3 + ci4, 1e-30)
    mu = max(float(props.VISVST(temp)), 1e-30)
    cf1 = props.CFTBL1(temp)
    cf2 = props.CFTBL2(temp)
    cf3 = props.CFTBL3(temp)
    cf4 = props.CFTBL4(temp)
    cfbar = (ci1 * cf1 + ci2 * cf2 + ci3 * cf3 + ci4 * cf4) / rho

    di3 = cfg.DIF3 * 14.7 / p * (temp / 492.0) ** 1.823
    di4 = cfg.DIF4 * 14.7 / p * (temp / 492.0) ** 1.823
    kc3 = _kcf(G, rho, mu, di3, AP)
    kc4 = _kcf(G, rho, mu, di4, AP)
    hc = 0.74 * G * cfbar * (G / (AP * mu)) ** (-0.41)

    dp3 = _dp3f(temp, cfg.DIF3, p)
    h4 = props.H4TBL(temp)
    h3 = props.H3TBL(temp)
    h3p = h3
    dp3p = dp3
    tps = temp
    tpsp = 0.0
    tpspp = 0.0
    x0 = 0.0
    x0p = 0.0
    cmcpn = ci3 / 2.0
    gamma = cfg.BGM / max(temp, 1e-12)
    beta = 0.0
    k0 = 0.0
    order_n = float(getattr(cfg, "EN2", 1.0))
    h2_floor = max(ci1, 1e-30)

    found_x0 = False
    for waf1 in np.arange(0.8, 1.0, 0.05):
        waf2 = 1.0 - waf1
        ltflg = 0
        lp1 = 1
        cps = ci3 / 2.0
        cmcpn = ci3 - cps
        dcpdx = kc3 / max(dp3, 1e-30) * (ci3 - cps)
        h3 = props.H3TBL(temp)
        h3p = h3
        dp3p = dp3
        tps = temp
        tpsp = 0.0
        tpspp = 0.0
        x0 = 0.0
        x0p = 0.0

        while lp1 <= 25:
            if lp1 != 1:
                tpspp = tpsp
                tpsp = tps

            tps = temp - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
            tps = max(tps, 1.0)
            h3 = props.H3TBL(tps)
            dp3 = _dp3f(tps, cfg.DIF3, p)
            dp3p = dp3
            h3p = h3
            tmtpn = temp - tps
            gamma = cfg.BGM / max(tps, 1e-12)
            beta = -cps * h3 * dp3 / max(cfg.KP * tps, 1e-30)
            k0 = cfg.ALPHA2 * np.exp(-gamma) * h2_floor ** cfg.EN3

            x0p = x0
            x0 = A - cps / max(dcpdx, 1e-30)
            if x0 < 0.0:
                x0 = 0.0
                x0a = 0.0
                denom = 1.0 + dp3 / max(A * kc3, 1e-30)
                cps = ci3 / denom
                dcpdx = ci3 / max(A, 1e-30)
                tps = temp - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
                tps = max(tps, 1.0)
            else:
                x0a = x0 / max(A, 1e-30)

            riesum = _trapp(x0a, 1.0, 50, x0a, cps, ci3, gamma, beta, k0, order_n)
            cpsp = cps
            cmcpo = cmcpn
            cps = ci3 - A * riesum / max(kc3, 1e-30)

            if ltflg != 1:
                if cps <= 0.25 * ci3:
                    ltflg = 1
                    x00 = waf1 * x0p + waf2 * x0
                    denom = kc3 * A - kc3 * x00
                    if abs(denom) < 1e-30:
                        cps = 0.0
                    else:
                        cps = ci3 / (1.0 + dp3 / denom)
                    dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cps)
                    cmcpn = ci3 - cps
                    h3 = h3p
                    lp1 += 1
                    continue
                cmcpn = ci3 - cps
            else:
                ltflg = 0
                if cps < 0.0:
                    cps = 0.0
                    cps = 0.2 * cps + 0.8 * cpsp
                    dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cps)
                    cmcpn = ci3 - cps
                    h3 = h3p
                    lp1 += 1
                    continue
                cmcpn = ci3 - cps

            dcpdx = kc3 / max(dp3, 1e-30) * (ci3 - cps)
            tgrad = hc * (temp - tps)
            tpspp = tpsp
            tpsp = tps
            tmtpo = tmtpn

            tps = temp - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
            tps = max(tps, 1.0)
            h3 = props.H3TBL(tps)
            dp3 = _dp3f(tps, cfg.DIF3, p)
            tmtpn = temp - tps
            gamma = cfg.BGM / max(tps, 1e-12)
            beta = -cps * h3 * dp3 / max(cfg.KP * tps, 1e-30)
            k0 = cfg.ALPHA2 * np.exp(-gamma) * h2_floor ** cfg.EN3

            temp_ok = abs(tmtpo - tmtpn) / max(abs(tmtpn), 1e-30) <= 0.05
            conc_ok = abs(cmcpo - cmcpn) / max(abs(cmcpn), 1e-30) <= 0.05
            if temp_ok and conc_ok:
                found_x0 = True
                break

            temp_loop = (
                min(tps, tpsp, tpspp) < tpsp or max(tps, tpsp, tpspp) > tpsp
            )
            if temp_loop:
                tpspp = tpsp
                tpsp = tps
                tps = 0.5 * (tpsp + tpspp)
                h3 = props.H3TBL(tps)
                dp3 = _dp3f(tps, cfg.DIF3, p)
                dp3p = dp3
                tmtpn = temp - tps
                dcpdx = (hc * (temp - tps) - h4 * kc4 * ci4) / max(h3 * dp3, 1e-30)
                cpsp = cps
                cmcpo = cmcpn
                cps = ci3 - dp3 / max(kc3, 1e-30) * dcpdx
                cps = max(cps, 0.0)
                cmcpn = ci3 - cps
                lp1 += 1
                continue

            cps = 0.2 * cps + 0.8 * cpsp
            dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cps)
            cmcpn = ci3 - cps
            h3 = h3p
            lp1 += 1

        if found_x0:
            break

    if not found_x0:
        grad = kc3 * max(ci3 - cps, 0.0)
        tgrad = hc * (temp - max(tps, 1.0))
        return grad, tgrad, {"note": "x0 iteration did not converge", "kc3": kc3, "kc4": kc4}

    nx = 24
    nx1 = nx + 1
    delx0a = (1.0 - x0 / max(A, 1e-30)) / float(nx)
    cpx = np.zeros(nx1, dtype=float)
    pcpox = np.zeros(nx1, dtype=float)
    cpox = np.zeros(nx1, dtype=float)
    dx = np.zeros(nx1, dtype=float)
    rhet = np.zeros(nx1, dtype=float)
    lp2 = 1
    cpx_surface = cps

    while lp2 <= 300:
        x0a = x0 / max(A, 1e-30)
        xa = x0a
        vnu = -kc3 / max(dp3, 1e-30)
        ctrm = (A * vnu + 1.0) / (A * vnu)

        for i in range(nx1):
            if lp2 == 1:
                cpx[i] = _cpxf(xa, x0a, cps)
            rhet[i] = _rhetf(ci3, cpx[i], gamma, beta, k0, order_n)
            dx[i] = xa
            xa += delx0a

        for i in range(nx):
            cpx[i] = 0.5 * (cpx[i] + cpx[i + 1])
            rhet[i] = 0.5 * (rhet[i] + rhet[i + 1])

        dxl = x0a
        dxu = dxl + delx0a
        rr1 = 0.0
        for i in range(nx):
            rr1 += rhet[i] * (_eval2(dxl, dxu) - ctrm * _eval1(dxl, dxu))
            dxl = dxu
            dxu += delx0a
        cpox[0] = max(ci3 - A * A / max(dp3, 1e-30) * rr1, 0.0)

        int1 = 1
        for k in range(1, nx):
            r1 = 0.0
            r2 = 0.0
            ps1 = 0.0
            ps2 = 0.0
            x0a_k = x0 / max(A, 1e-30)
            xa = x0a_k + delx0a
            for i in range(int1):
                r1 += rhet[i] * _eval1(x0a_k, xa)
                x0a_k = xa
                xa += delx0a
            r1 *= (1.0 / max(x0a_k, 1e-30) - ctrm)
            xad = xa
            xa = xa - delx0a
            for i in range(int1, nx - 1):
                ps1 += rhet[i + 1] * _eval2(xa, xad)
                ps2 += rhet[i + 1] * _eval1(xa, xad)
                xa = xad
                xad += delx0a
            r2 = ps1 - ctrm * ps2
            cpox[k] = max(ci3 - A * A / max(dp3, 1e-30) * (r1 + r2), 0.0)
            int1 += 1

        dxl = x0 / max(A, 1e-30)
        dxu = dxl + delx0a
        rr2 = 0.0
        for i in range(nx):
            rr2 += rhet[i] * _eval1(dxl, dxu)
            dxl = dxu
            dxu += delx0a
        cpox[nx] = max(ci3 - A * A / max(dp3, 1e-30) * (1.0 - ctrm) * rr2, 0.0)

        dcpdx = kc3 / max(dp3, 1e-30) * (ci3 - cpox[nx])
        h3p = h3
        dp3p = dp3
        tps = temp - (h4 * kc4 * ci4 + h3 * dp3 * dcpdx) / max(hc, 1e-30)
        tps = max(tps, 1.0)
        h3 = props.H3TBL(tps)
        dp3 = _dp3f(tps, cfg.DIF3, p)
        tmtpo = temp - tps if lp2 == 1 else tmtpn
        tmtpn = temp - tps

        if lp2 > 1:
            cmcpo = cmcpn
            cmcpn = ci3 - cpox[nx]
            temp_ok = abs(tmtpo - tmtpn) / max(abs(tmtpn), 1e-30) <= 0.05
            conc_ok = abs(cmcpo - cmcpn) / max(abs(cmcpn), 1e-30) <= 0.05
            if temp_ok and conc_ok:
                cpx_surface = cpox[nx]
                break

        for i in range(nx1):
            if lp2 % 5 != 0:
                cpx[i] = 0.8 * cpx[i] + 0.2 * cpox[i]
            else:
                cpx[i] = 0.5 * (cpox[i] + pcpox[i])
            pcpox[i] = cpox[i]

        cmcpn = ci3 - cpx[nx]
        dcpdx = kc3 / max(dp3p, 1e-30) * (ci3 - cpx[nx])
        tps = temp - (h4 * kc4 * ci4 + h3p * dp3p * dcpdx) / max(hc, 1e-30)
        tps = max(tps, 1.0)
        h3 = props.H3TBL(tps)
        dp3 = _dp3f(tps, cfg.DIF3, p)
        tmtpn = temp - tps
        gamma = cfg.BGM / max(tps, 1e-12)
        beta = -cpx[nx] * h3 * dp3 / max(cfg.KP * tps, 1e-30)
        k0 = cfg.ALPHA2 * np.exp(-gamma) * h2_floor ** cfg.EN3
        cpx_surface = cpx[nx]
        lp2 += 1

    grad = kc3 / max(dp3, 1e-30) * (ci3 - cpx_surface) * dp3
    tgrad = hc * (temp - tps)
    diag = {
        "tps": tps,
        "x0": x0,
        "x0a": x0 / max(A, 1e-30),
        "cps": cpx_surface,
        "kc3": kc3,
        "kc4": kc4,
        "hc": hc,
        "gamma": gamma,
        "beta": beta,
        "k0": k0,
        "lp2": lp2,
    }
    if verbose:
        print(f"[SGRAD] x0={x0:.6e} cps={cpx_surface:.6e} tps={tps:.2f} lp2={lp2}")
    return grad, tgrad, diag
