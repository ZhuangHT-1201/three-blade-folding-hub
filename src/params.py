import numpy as np

NB = 3
R_TIP = 1.40
R_ROOT = 0.28
RHO_AIR = 1.225
RPM_NOM = 982.0
RPM_MAX = 1080.0
SF_ULT = 1.5
SF_FAT = 1.5
GUST_FACTOR = 1.3
ALT_FLAP_RATIO = 0.35

SKIN_T = 1.0e-3
RHO_SKIN = 1550.0
SPAR_FRAC = 0.08
RHO_SPAR = 1600.0
RHO_CORE = 75.0
PAINT_FACTOR = 1.05

MAT = {
    'Al7075-T6': dict(E=71.7e9, nu=0.33, rho=2810., sy=503e6, su=572e6, tau_u=331e6, se=159e6,
                      name='7075-T6 aluminium alloy'),
    'Ti6Al4V':   dict(E=113.8e9, nu=0.34, rho=4430., sy=880e6, su=950e6, tau_u=550e6, se=510e6,
                      name='Ti-6Al-4V titanium alloy'),
    '42CrMo':    dict(E=210e9, nu=0.30, rho=7850., sy=930e6, su=1080e6, tau_u=650e6, se=480e6,
                      name='42CrMo quenched and tempered steel'),
    '17-4PH':    dict(E=196e9, nu=0.27, rho=7800., sy=1000e6, su=1070e6, tau_u=670e6, se=500e6,
                      name='17-4PH (H1025) stainless steel'),
}
MAT_PLATE = 'Al7075-T6'
MAT_GRIP = 'Al7075-T6'
MAT_PIN = '17-4PH'

PARK_AZIMUTH_DEG = 180.0
CLEARANCE_MIN = 0.012


def flip(P):
    Q = np.array(P, dtype=float).copy()
    Q[..., 1] *= -1
    Q[..., 2] *= -1
    return Q
