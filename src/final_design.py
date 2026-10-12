import json, os
import numpy as np

import params as P
import fold_kinematics as K
import sizing as S
import loads as L

HERE = os.path.dirname(os.path.abspath(__file__))

rows = [r for r in json.load(open(os.path.join(HERE, 'pareto.json'))) if r.get('feasible')]
m_min = min(r['mass'] for r in rows)
a_min = min(r['area'] for r in rows)
for r in rows:
    r['J'] = r['mass'] / m_min + r['area'] / a_min
best = min(rows, key=lambda r: r['J'])
r_h, s = best['r_h'], best['s']
print('selected r_h %.2f s %.2f (mass %.3f, area %.3f)' % (r_h, s, best['mass'], best['area']))

sol = S.optimise(r_h=r_h, s=s)
xr = S.round_design(sol['x'])
x = dict(zip(S.NAMES, xr))
gm = S.geometry_mass(xr, r_h, s)
z_h = S.balance_offset(r_h, gm['m_grip'], gm['r_cg_grip'])
z_h = float(np.round(z_h, 3))
LL = S.design_loads(r_h, s, gm['m_grip'], gm['r_cg_grip'], z_h)
ut = S.constraints(xr, r_h, s, LL)

fixed = dict(R_t=x['R_t'], R_p=x['R_p'], a_t=x['a_t'], t_t=x['t_t'], t_p=x['t_p'], z_h=z_h)
(lay, env, clr, sweep), _ = K.optimise(r_list=[r_h], s_list=[s], y_list=np.arange(-0.035, 0.0051, 0.005),
                                       d_list=np.arange(0.0, 4.01, 0.5), **fixed)
dep = K.deployed_envelope(lay)
print('y_h %.3f delta %.1f sweep %.1f mm' % (lay.y_h, lay.delta_deg, sweep * 1e3))
print('folded bbox', env['bbox'], 'plan', env['plan_area'])

b_nom, b_max = L.bemt(P.RPM_NOM), L.bemt(P.RPM_MAX)
design = dict(
    r_h=r_h, s=s, y_h=float(lay.y_h), delta_deg=float(lay.delta_deg), z_h=z_h,
    R_core=lay.R_core, gap=S.GAP, t_cuff=S.T_CUFF, L_cuff=S.L_CUFF, t_bulkhead=S.T_BULKHEAD,
    **{k: float(v) for k, v in x.items()},
    loads={k: float(v) for k, v in LL.items()}, util={k: float(v) for k, v in ut.items()},
    mass_est={k: float(v) for k, v in gm.items()},
    blade_mass=float(L.BLADE_MASS), blade_rcg=float(L.BLADE_RCG),
    rotor=dict(R=P.R_TIP, r_root=P.R_ROOT, rpm_nom=P.RPM_NOM, rpm_max=P.RPM_MAX,
               T_nom=b_nom['T'], Q_nom=b_nom['Q'], P_nom=b_nom['Pw'], T_max=b_max['T'], CT=b_nom['CT']),
    folded=dict(plan_area=env['plan_area'], bbox=env['bbox'].tolist(), bbox_vol=env['bbox_vol'],
                clear={k: float(v) for k, v in clr.items()}, sweep_min=float(sweep)),
    deployed=dict(plan_area=dep['plan_area'], bbox=dep['bbox'].tolist(), bbox_vol=dep['bbox_vol']),
    pareto_pick=best,
)
json.dump(design, open(os.path.join(HERE, 'design.json'), 'w'), indent=1, ensure_ascii=False)
print(json.dumps({k: v for k, v in design.items() if not isinstance(v, dict)}, indent=1))
