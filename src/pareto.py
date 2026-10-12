import json, os, time
import numpy as np

import params as P
import fold_kinematics as K
import sizing as S

HERE = os.path.dirname(os.path.abspath(__file__))

rows = []
t0 = time.time()
for r_h in np.round(np.arange(0.12, 0.2001, 0.01), 3):
    for s in np.round(np.arange(0.04, 0.0801, 0.01), 3):
        sol = S.optimise(r_h=r_h, s=s)
        x = dict(zip(S.NAMES, S.round_design(sol['x'])))
        if max(sol['util'].values()) > 1.02:
            rows.append(dict(r_h=r_h, s=s, feasible=False)); continue
        fixed = dict(R_t=x['R_t'], R_p=x['R_p'], a_t=x['a_t'], t_t=x['t_t'], t_p=x['t_p'], z_h=sol['z_h'])
        try:
            (lay, env, c, w), _ = K.optimise(r_list=[r_h], s_list=[s], **fixed)
        except RuntimeError:
            rows.append(dict(r_h=r_h, s=s, feasible=False, mass=sol['mass']['total'])); continue
        rows.append(dict(r_h=r_h, s=s, feasible=True, mass=sol['mass']['total'], y_h=lay.y_h,
                         delta=lay.delta_deg, area=env['plan_area'], bbox=env['bbox'].tolist(),
                         bbox_vol=env['bbox_vol'], sweep=w, z_h=sol['z_h'], x=x))
        print('r_h %.2f s %.2f  mass %.3f  area %.3f  bbox %s  y_h %.3f  dl %.0f  (%.0fs)'
              % (r_h, s, sol['mass']['total'], env['plan_area'], np.round(env['bbox'], 3), lay.y_h, lay.delta_deg,
                 time.time() - t0), flush=True)
json.dump(rows, open(os.path.join(HERE, 'pareto.json'), 'w'), indent=1)
