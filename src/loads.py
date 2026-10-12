import json, os
import numpy as np
from scipy.optimize import brentq

import params as P
from blade_geom import T as T_RAW, section, section_props

HERE = os.path.dirname(os.path.abspath(__file__))
T_BLADE = P.flip(T_RAW)

def build_sections(n=60):
    xs = np.linspace(P.R_ROOT + 5e-4, P.R_TIP - 5e-4, n)
    rows = []
    for x in xs:
        s = section(T_BLADE, x)
        sp = section_props(s)
        pts = s.reshape(-1, 2)
        per = np.linalg.norm(s[:, 1] - s[:, 0], axis=1).sum()
        rows.append(dict(r=x, c=sp['chord'], A=sp['A'], cy=sp['cy'], cz=sp['cz'],
                         I1=sp['I1'], I2=sp['I2'], per=per,
                         ymin=pts[:, 0].min(), ymax=pts[:, 0].max(),
                         zmin=pts[:, 1].min(), zmax=pts[:, 1].max()))
    return rows


SEC = build_sections()
r = np.array([s['r'] for s in SEC])
c = np.array([s['c'] for s in SEC])
A = np.array([s['A'] for s in SEC])
per = np.array([s['per'] for s in SEC])
cy = np.array([s['cy'] for s in SEC])
cz = np.array([s['cz'] for s in SEC])
TH_ROOT, TH_TIP = np.radians(29.30), np.radians(6.80)
theta = TH_ROOT + (TH_TIP - TH_ROOT) * (r - r[0]) / (r[-1] - r[0])


def mass_per_length():
    skin = P.RHO_SKIN * P.SKIN_T * per
    spar = P.RHO_SPAR * P.SPAR_FRAC * A
    core = P.RHO_CORE * np.clip(A - P.SKIN_T * per - P.SPAR_FRAC * A, 0, None)
    return P.PAINT_FACTOR * (skin + spar + core)


MPL = mass_per_length()
BLADE_MASS = np.trapezoid(MPL, r)
BLADE_RCG = np.trapezoid(MPL * r, r) / BLADE_MASS


def polar(alpha):
    a0, a_s = 5.73, np.radians(13.0)
    al = np.abs(alpha)
    cl = np.where(al < a_s, a0 * al, a0 * a_s * (1 - 0.6 * np.clip((al - a_s) / np.radians(10), 0, 1)))
    cl = np.sign(alpha) * cl
    cd = 0.0087 - 0.0216 * al + 0.40 * al ** 2
    cd = np.where(al < a_s, cd, cd + 1.0 * (al - a_s))
    return cl, cd


def bemt(rpm):
    Om = rpm * 2 * np.pi / 60
    Vt = Om * P.R_TIP
    rr = r / P.R_TIP
    sig = P.NB * c / (np.pi * P.R_TIP)
    lam = np.zeros_like(rr)
    for k, (x, s_, th) in enumerate(zip(rr, sig, theta)):
        def f(l):
            phi = np.arctan2(l, x)
            cl, cd = polar(th - phi)
            fP = P.NB / 2 * (1 - x) / max(l / x, 1e-6) if x < 1 else 0
            F = 2 / np.pi * np.arccos(np.exp(-fP)) if x < 1 else 1e-3
            mom = 4 * F * l ** 2 * x
            be = 0.5 * s_ * (cl * np.cos(phi) - cd * np.sin(phi)) * (x ** 2 + l ** 2)
            return mom - be
        lam[k] = brentq(f, 1e-6, 0.6)
    phi = np.arctan2(lam, rr)
    alpha = theta - phi
    cl, cd = polar(alpha)
    W2 = Vt ** 2 * (rr ** 2 + lam ** 2)
    q = 0.5 * P.RHO_AIR * W2 * c
    dT = q * (cl * np.cos(phi) - cd * np.sin(phi))
    dD = q * (cl * np.sin(phi) + cd * np.cos(phi))
    T1 = np.trapezoid(dT, r)
    Q1 = np.trapezoid(dD * r, r)
    return dict(rpm=rpm, Om=Om, Vt=Vt, dT=dT, dD=dD, alpha=alpha, lam=lam,
                T=P.NB * T1, Q=P.NB * Q1, Pw=P.NB * Q1 * Om,
                CT=P.NB * T1 / (P.RHO_AIR * np.pi * P.R_TIP ** 2 * Vt ** 2))


def hinge_loads(rpm, r_h, m_grip, r_cg_grip, thrust_factor=1.0):
    b = bemt(rpm)
    Om2 = b['Om'] ** 2
    fx = MPL * Om2 * r
    fy = MPL * Om2 * cy - b['dD']
    fz = b['dT'] * thrust_factor
    pos = np.stack([r - r_h, cy, cz], 1)
    f = np.stack([fx, fy, fz], 1)
    F = np.trapezoid(f, r, axis=0)
    M = np.trapezoid(np.cross(pos, f), r, axis=0)
    k2 = (np.array([s['I1'] + s['I2'] for s in SEC]) / A)
    Ixx_flat = MPL * k2 * 0.97
    M[0] += -np.trapezoid(Om2 * Ixx_flat * np.sin(theta) * np.cos(theta), r)
    F[0] += m_grip * Om2 * r_cg_grip
    return dict(F=F, M=M, bemt=b)


if __name__ == '__main__':
    print(f'blade mass {BLADE_MASS:.3f} kg, r_cg {BLADE_RCG:.3f} m')
    print(f'chord root {c[0]*1e3:.1f} tip {c[-1]*1e3:.1f} mm')
    tab = []
    for rpm in (600, 700, 800, 900, 1000, 1100, 1200):
        b = bemt(rpm)
        tab.append(dict(rpm=rpm, Vt=b['Vt'], T=b['T'], Q=b['Q'], P=b['Pw'], CT=b['CT'],
                        amax=float(np.degrees(b['alpha'].max())), amin=float(np.degrees(b['alpha'].min()))))
        print(f"rpm {rpm}: Vt {b['Vt']:.1f} T {b['T']:.0f} N ({b['T']/9.81:.0f} kgf) Q {b['Q']:.0f} N.m P {b['Pw']/1e3:.1f} kW "
              f"CT {b['CT']:.5f} alpha {np.degrees(b['alpha'].min()):.1f}..{np.degrees(b['alpha'].max()):.1f} deg")
    out = dict(blade_mass=BLADE_MASS, blade_rcg=BLADE_RCG, table=tab,
               r=r.tolist(), c=c.tolist(), mpl=MPL.tolist(), theta=np.degrees(theta).tolist(),
               cy=cy.tolist(), cz=cz.tolist(), A=A.tolist())
    for rpm in (P.RPM_NOM, P.RPM_MAX):
        L = hinge_loads(rpm, 0.20, 0.6, 0.25)
        print(rpm, 'F', L['F'], 'M', L['M'])
    json.dump(out, open(os.path.join(HERE, 'loads_base.json'), 'w'), indent=1)
