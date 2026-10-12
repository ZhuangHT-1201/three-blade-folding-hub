import json, os
import sizing as S

HERE = os.path.dirname(os.path.abspath(__file__))
D = json.load(open(os.path.join(HERE, 'design.json'), encoding='utf-8'))
opt = S.optimise(r_h=D['r_h'], s=D['s'])
ref = S.optimise(z_h=0.0, r_h=D['r_h'], s=D['s'])
out = dict(opt_mass=float(opt['mass']['total']), ref_mass=float(ref['mass']['total']),
           ref_x={n: float(v) for n, v in zip(S.NAMES, ref['x'])}, ref_LL={k: float(v) for k, v in ref['LL'].items()})
json.dump(out, open(os.path.join(HERE, 'sizing_ref.json'), 'w'), indent=1)
print(out)
