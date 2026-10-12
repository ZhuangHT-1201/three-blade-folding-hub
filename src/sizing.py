import json, os
import numpy as np
from scipy.optimize import minimize

import params as P
import loads as L

HERE = os.path.dirname(os.path.abspath(__file__))

MG, MP, MPIN = P.MAT[P.MAT_GRIP], P.MAT[P.MAT_PLATE], P.MAT[P.MAT_PIN]
QR_PIN_KD = 0.85
KT_LUG = 2.6
GAP = 0.5e-3
T_CUFF = 3.0e-3
L_CUFF = 0.06
T_BULKHEAD = 6.0e-3
ROOT_PER = 0.502
ROOT_AREA = 4.97e-3
THETA_ALLOW_DEG = 0.25
P_FRET = 100e6
D_MIN = 0.010
NAMES = ['d', 'k', 't_t', 'R_t', 'a_t', 't_p', 'R_p', 'b_a']


def layout():
    f = os.path.join(HERE, 'layout_opt.json')
    if os.path.exists(f):
        j = json.load(open(f))
        return j['r_h'], j['s'], j.get('y_h', 0.0)
    return 0.20, 0.08, 0.0


def aero_flap_moment(rpm, r_h):
    b = L.bemt(rpm)
    return np.trapezoid(b['dT'] * (L.r - r_h), L.r), b


def blade_Izz(r_h):
    return np.trapezoid(L.MPL * (L.r) ** 2, L.r)


def design_loads(r_h, s, m_grip, r_cg_grip, z_h):
    out = {}
    for tag, rpm in (('lim', P.RPM_MAX), ('nom', P.RPM_NOM)):
        H = L.hinge_loads(rpm, r_h, m_grip, r_cg_grip)
        F, M = H['F'].copy(), H['M'].copy()
        M[1] -= z_h * F[0]
        Ma, _ = aero_flap_moment(rpm, r_h)
        out[tag] = dict(F=F, M=M, Ma=Ma, Om=H['bemt']['Om'])
    lim = out['lim']; nom = out['nom']
    spin_acc = lim['Om'] / 3.0
    Mz_spin = (blade_Izz(r_h) + m_grip * r_cg_grip ** 2) * spin_acc
    LL = dict(
        Fx=lim['F'][0], Fy=abs(lim['F'][1]), Fz=abs(lim['F'][2]) * P.GUST_FACTOR,
        My=abs(lim['M'][1]) + (P.GUST_FACTOR - 1 + P.ALT_FLAP_RATIO) * abs(lim['Ma']),
        Mz=abs(lim['M'][2]) + Mz_spin, Mx=abs(lim['M'][0]),
        Fx_m=nom['F'][0], My_m=abs(nom['M'][1]), My_a=P.ALT_FLAP_RATIO * abs(nom['Ma']),
        Mz_m=abs(nom['M'][2]), Mz_a=0.3 * abs(nom['M'][2]) + 0.05 * abs(nom['Ma']),
    )
    return LL


def balance_offset(r_h, m_grip, r_cg_grip):
    H = L.hinge_loads(P.RPM_NOM, r_h, m_grip, r_cg_grip)
    return H['M'][1] / H['F'][0]


def geometry_mass(x, r_h, s):
    d, k, t_t, R_t, a_t, t_p, R_p, b_a = x
    W_t = s + 2 * R_t
    W_p = s + 2 * R_p
    L_t = (P.R_ROOT - r_h) + a_t
    v_pin = np.pi / 4 * d ** 2 * (1 - k ** 2) * (t_t + 2 * t_p + 1.2 * d)
    v_tongue = t_t * W_t * L_t - 2 * np.pi / 4 * d ** 2 * t_t
    v_cuff = ROOT_PER * L_CUFF * T_CUFF + ROOT_AREA * T_BULKHEAD
    R_core = 0.075
    v_plate = t_p * (np.pi * R_core ** 2 + 3 * (0.5 * (b_a + W_p) * (r_h - R_core) + W_p * R_p * 0.8)
                     - 3 * 2 * np.pi / 4 * d ** 2)
    m_pin = MPIN['rho'] * v_pin
    m_grip = MG['rho'] * (v_tongue + v_cuff)
    r_cg_grip = (r_h * v_tongue + (P.R_ROOT + L_CUFF / 2) * v_cuff) / (v_tongue + v_cuff)
    m_plates = 2 * MP['rho'] * v_plate
    total = 3 * (2 * m_pin + m_grip) + m_plates
    return dict(total=total, m_pin=m_pin, m_grip=m_grip, m_plates=m_plates, r_cg_grip=r_cg_grip,
                W_t=W_t, W_p=W_p)


def constraints(x, r_h, s, LL):
    d, k, t_t, R_t, a_t, t_p, R_p, b_a = x
    U = P.SF_ULT
    W_t, W_p = s + 2 * R_t, s + 2 * R_p
    Px = LL['Fx'] / 2 + LL['Mz'] / s
    Py = LL['Fy'] / 2
    Pp = np.hypot(Px, Py)
    A_pin = np.pi / 4 * d ** 2 * (1 - k ** 2)
    Z_pin = np.pi / 32 * d ** 3 * (1 - k ** 4)
    h = t_t + t_p + 2 * GAP
    Mf = LL['My'] / 2
    M_pin = Pp * h / 4 + Mf / 2
    V_pin = Pp / 2 + Mf / h
    c = {}
    c['pin_shear'] = U * V_pin / A_pin / (MPIN['tau_u'] * QR_PIN_KD)
    c['pin_bend'] = U * M_pin / Z_pin / (MPIN['su'] * QR_PIN_KD)
    c['pin_yield'] = np.sqrt((M_pin / Z_pin) ** 2 + 3 * (4 / 3 * V_pin / A_pin) ** 2) / MPIN['sy']
    c['tongue_shearout'] = U * Pp / (2 * t_t * (a_t - d / 2 * np.cos(np.radians(40))) * MG['tau_u'])
    c['tongue_bearing'] = U * (Pp / (d * t_t) + 6 * Mf / (d * t_t ** 2)) / (1.5 * MG['su'])
    c['tongue_net'] = U * LL['Fx'] / ((W_t - 2 * d) * t_t) / MG['su']
    sec_t = (U * (LL['Fx'] / (W_t * t_t) + 6 * LL['My'] / (W_t * t_t ** 2) + 6 * LL['Mz'] / (t_t * W_t ** 2)))
    c['tongue_root'] = sec_t / MG['su']
    Ppl = Pp / 2 + Mf / h
    c['plate_shearout'] = U * Ppl / (2 * t_p * (R_p - d / 2 * np.cos(np.radians(40))) * MP['tau_u'])
    c['plate_bearing'] = U * Ppl / (d * t_p) / (1.5 * MP['su'])
    c['plate_net'] = U * 2 * Ppl / ((W_p - 2 * d) * t_p) / MP['su']
    sig_arm = (LL['Fx'] / 2 + LL['My'] / h) / (b_a * t_p) + 6 * (LL['Mz'] / 2) / (t_p * b_a ** 2)
    c['plate_arm'] = U * sig_arm / MP['su']
    def good(sm, sa, mat):
        return P.SF_FAT * (sa / mat['se'] + sm / mat['su'])
    Pm = LL['Fx_m'] / 2 + LL['Mz_m'] / s
    Pa = LL['Mz_a'] / s
    c['fat_plate_hole'] = good(KT_LUG * Pm / 2 / ((W_p / 2 - d) * t_p), KT_LUG * Pa / 2 / ((W_p / 2 - d) * t_p), MP)
    c['fat_tongue_hole'] = good(KT_LUG * Pm / ((W_t / 2 - d) * t_t), KT_LUG * Pa / ((W_t / 2 - d) * t_t), MG)
    c['fat_tongue_root'] = good(1.3 * (LL['Fx_m'] / (W_t * t_t) + 6 * LL['My_m'] / (W_t * t_t ** 2)),
                                1.3 * 6 * LL['My_a'] / (W_t * t_t ** 2), MG)
    c['fat_plate_arm'] = good(1.3 * (LL['Fx_m'] / 2 + LL['My_m'] / h) / (b_a * t_p),
                              1.3 * (LL['My_a'] / h) / (b_a * t_p), MP)
    Mpm = Pm * h / 4 + LL['My_m'] / 4
    Mpa = Pa * h / 4 + LL['My_a'] / 4
    c['fat_pin'] = good(1.5 * Mpm / Z_pin, 1.5 * Mpa / Z_pin, MPIN)
    c['fat_plate_lug'] = good(KT_LUG * (Pm / 2 + LL['My_m'] / 2 / h) / ((W_p / 2 - d) * t_p),
                              KT_LUG * (Pa / 2 + LL['My_a'] / 2 / h) / ((W_p / 2 - d) * t_p), MP)
    th_allow = np.radians(THETA_ALLOW_DEG)
    I_t = W_t * t_t ** 3 / 12
    c['stiff_tongue'] = LL['My'] * (P.R_ROOT - r_h) / (MG['E'] * I_t) / th_allow
    I_a = 2 * b_a * t_p * (h / 2) ** 2 + 2 * b_a * t_p ** 3 / 12
    c['stiff_arm'] = LL['My'] * (r_h - 0.075) / (MP['E'] * I_a) / th_allow
    c['fret_tongue'] = (Pm / (d * t_t) + 6 * LL['My_m'] / 2 / (d * t_t ** 2)) / P_FRET
    c['fret_plate'] = (Pm / 2 + LL['My_m'] / 2 / h) / (d * t_p) / P_FRET
    c['edge_t'] = 1.25 * d / a_t
    c['edge_lobe_t'] = 1.0 * d / R_t
    c['edge_p'] = 1.25 * d / R_p
    c['pin_space'] = 2.5 * d / s
    c['cuff_clear'] = (R_p + P.CLEARANCE_MIN + 0.02) / (P.R_ROOT - r_h)
    c['arm_width'] = b_a / (W_p * 1.2)
    return c


def optimise(z_h=None, r_h=None, s=None):
    r0, s0, y_h = layout()
    r_h = r0 if r_h is None else r_h
    s = s0 if s is None else s
    x0 = np.array([0.018, 0.0, 0.026, 0.022, 0.030, 0.006, 0.022, 0.10])
    lb = np.array([D_MIN, 0.0, 0.010, 0.010, 0.012, 0.004, 0.010, 0.030])
    ub = np.array([0.030, 0.6, 0.050, 0.050, 0.060, 0.020, 0.050, 0.200])
    scale = ub
    m_grip, r_cg = 0.6, 0.25
    for it in range(4):
        zh = balance_offset(r_h, m_grip, r_cg) if z_h is None else z_h
        LL = design_loads(r_h, s, m_grip, r_cg, zh)
        obj = lambda xs: geometry_mass(xs * scale, r_h, s)['total']
        cons = [{'type': 'ineq', 'fun': (lambda xs, k=k: 1 - constraints(xs * scale, r_h, s, LL)[k])}
                for k in constraints(x0, r_h, s, LL)]
        res = minimize(obj, x0 / scale, method='SLSQP', bounds=list(zip(lb / scale, ub / scale)),
                       constraints=cons, options=dict(maxiter=500, ftol=1e-10))
        x0 = res.x * scale
        gm = geometry_mass(x0, r_h, s)
        m_grip, r_cg = gm['m_grip'], gm['r_cg_grip']
    return dict(x=x0, res=res, LL=LL, z_h=zh, mass=gm, r_h=r_h, s=s, y_h=y_h,
                util=constraints(x0, r_h, s, LL))


def round_design(x):
    d, k, t_t, R_t, a_t, t_p, R_p, b_a = x
    up = lambda v, q: np.ceil(v / q - 1e-9) * q
    d = up(d, 0.002)
    R_t, a_t, R_p = max(R_t, d), max(a_t, 1.25 * d), max(R_p, 1.25 * d)
    return np.array([d, np.floor(k * 10) / 10, up(t_t, 0.001), up(R_t, 0.001), up(a_t, 0.001),
                     up(t_p, 0.0005), up(R_p, 0.001), up(b_a, 0.005)])


def report(sol, title):
    print('=====', title)
    print('z_h = %.1f mm' % (sol['z_h'] * 1000))
    for k, v in sol['LL'].items():
        print('  %-6s %10.1f' % (k, v))
    for n, v in zip(NAMES, sol['x']):
        print('  %-4s %8.2f %s' % (n, v * (1 if n == 'k' else 1000), '' if n == 'k' else 'mm'))
    m = sol['mass']
    print('  mass total %.3f kg (pin %.3f x6, grip %.3f x3, plates %.3f)' % (m['total'], m['m_pin'], m['m_grip'], m['m_plates']))
    act = {k: round(v, 3) for k, v in sol['util'].items() if v > 0.9}
    print('  active constraints', act)


if __name__ == '__main__':
    sol = optimise()
    report(sol, 'optimised, CF-balanced hinge plane')
    xr = round_design(sol['x'])
    gm = geometry_mass(xr, sol['r_h'], sol['s'])
    LLr = design_loads(sol['r_h'], sol['s'], gm['m_grip'], gm['r_cg_grip'], sol['z_h'])
    ut = constraints(xr, sol['r_h'], sol['s'], LLr)
    rounded = dict(x=xr, LL=LLr, z_h=sol['z_h'], mass=gm, util=ut, r_h=sol['r_h'], s=sol['s'], y_h=sol['y_h'])
    report(rounded, 'rounded (manufacturable)')
    print('  max util %.3f (%s)' % (max(ut.values()), max(ut, key=ut.get)))

    base = optimise(z_h=0.0)
    report(base, 'reference: hinge plane on pitch axis (no CF balance)')

    def ser(s_):
        return dict(x={n: float(v) for n, v in zip(NAMES, s_['x'])}, z_h=float(s_['z_h']),
                    LL={k: float(v) for k, v in s_['LL'].items()},
                    mass={k: float(v) for k, v in s_['mass'].items()},
                    util={k: float(v) for k, v in s_['util'].items()},
                    r_h=s_['r_h'], s=s_['s'], y_h=s_['y_h'])
    json.dump(dict(opt=ser(sol), rounded=ser(rounded), ref=ser(base)),
              open(os.path.join(HERE, 'sizing.json'), 'w'), indent=1)
