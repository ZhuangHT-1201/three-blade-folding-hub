import json, os, itertools
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial import ConvexHull

import params as P
from loads import T_BLADE

HERE = os.path.dirname(os.path.abspath(__file__))


def rotz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])


def blade_cloud():
    tri = T_BLADE
    pts = np.concatenate([tri.reshape(-1, 3), tri.mean(1)])
    return np.unique(np.round(pts, 5), axis=0)


BLADE_PTS = blade_cloud()


def cuff_cloud(t_cuff=0.004, length=0.06):
    sel = BLADE_PTS[:, 0] < P.R_ROOT + length
    p = BLADE_PTS[sel].copy()
    cen = np.array([0, -0.03, -0.012])
    d = p[:, 1:] - cen[1:]
    n = d / (np.linalg.norm(d, axis=1, keepdims=True) + 1e-12)
    p[:, 1:] += n * t_cuff
    return p


def box_cloud(u0, u1, v0, v1, z0, z1, h=0.004):
    us = np.arange(u0, u1 + 1e-9, h); vs = np.arange(v0, v1 + 1e-9, h); zs = np.arange(z0, z1 + 1e-9, h)
    pts = []
    for a, b in ((u0, u1),):
        for uu in (a, b):
            V, Z = np.meshgrid(vs, zs); pts.append(np.stack([np.full(V.size, uu), V.ravel(), Z.ravel()], 1))
    for vv in (v0, v1):
        U, Z = np.meshgrid(us, zs); pts.append(np.stack([U.ravel(), np.full(U.size, vv), Z.ravel()], 1))
    for zz in (z0, z1):
        U, V = np.meshgrid(us, vs); pts.append(np.stack([U.ravel(), V.ravel(), np.full(U.size, zz)], 1))
    return np.concatenate(pts)


def disc_cloud(R, z0, z1, h=0.004):
    pts = []
    for rr in np.arange(h, R + 1e-9, h):
        n = max(6, int(2 * np.pi * rr / h))
        a = np.linspace(0, 2 * np.pi, n, endpoint=False)
        for zz in (z0, z1):
            pts.append(np.stack([rr * np.cos(a), rr * np.sin(a), np.full(n, zz)], 1))
    n = int(2 * np.pi * R / h)
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    for zz in np.arange(z0, z1 + 1e-9, h):
        pts.append(np.stack([R * np.cos(a), R * np.sin(a), np.full(n, zz)], 1))
    return np.concatenate(pts)


class Layout:
    def __init__(self, r_h=0.20, s=0.08, delta_deg=3.0, y_h=0.0, z_h=-0.046,
                 R_t=0.012, R_p=0.014, a_t=0.015, t_t=0.024, t_p=0.004, R_core=0.060, H_core=0.10):
        self.__dict__.update(locals()); del self.__dict__['self']
        self.phi = np.radians(120.0 + delta_deg)

    def local_parts(self):
        blade = BLADE_PTS
        cuff = cuff_cloud()
        L = P.R_ROOT - self.r_h
        tongue = box_cloud(-self.a_t, L, self.y_h - self.s / 2 - self.R_t, self.y_h + self.s / 2 + self.R_t,
                           self.z_h - self.t_t / 2, self.z_h + self.t_t / 2)
        tongue[:, 0] += self.r_h
        plate = box_cloud(self.R_core, self.r_h + self.R_p, self.y_h - self.s / 2 - self.R_p,
                          self.y_h + self.s / 2 + self.R_p,
                          self.z_h - self.t_t / 2 - self.t_p, self.z_h + self.t_t / 2 + self.t_p)
        return dict(blade=blade, cuff=cuff, tongue=tongue), plate

    def core(self):
        z1 = self.z_h + self.t_t / 2 + self.t_p
        return disc_cloud(self.R_core, z1 - self.H_core, z1 + 0.02)

    def pins_local(self):
        return np.array([[self.r_h, self.y_h + self.s / 2, 0], [self.r_h, self.y_h - self.s / 2, 0]])

    def placed(self, frac=1.0):
        loc, plate = self.local_parts()
        out, plates = {}, []
        for name, psi, pin_idx, sgn in (('A', 180, None, 0), ('B', 60, 0, +1), ('C', 300, 1, -1)):
            Rz = rotz(np.radians(psi))
            pins = self.pins_local() @ Rz.T
            moving = np.concatenate(list(loc.values())) @ Rz.T
            plates.append(plate @ Rz.T)
            if pin_idx is not None:
                piv = pins[pin_idx]
                Rf = rotz(sgn * self.phi * frac)
                moving = (moving - piv) @ Rf.T + piv
            out[name] = moving
        return out, plates

    def clearances(self, frac=1.0, pairs=(('A', 'B'), ('A', 'C'), ('B', 'C'))):
        sets, plates = self.placed(frac)
        res = {}
        treeA = cKDTree(sets['A'])
        for a, b in pairs:
            t = treeA if a == 'A' else cKDTree(sets[a])
            res[a + b] = t.query(sets[b], k=1, distance_upper_bound=0.1)[0].min()
        core = self.core()
        for k, (name) in enumerate(('A', 'B', 'C')):
            others = np.concatenate([plates[j] for j in range(3) if j != k] + [core])
            res[name + '-hub'] = cKDTree(others).query(sets[name], k=1, distance_upper_bound=0.1)[0].min()
        return {k: min(v, 0.1) for k, v in res.items()}

    def envelope(self, frac=1.0):
        sets, plates = self.placed(frac)
        allp = np.concatenate(list(sets.values()) + plates)
        hull = ConvexHull(allp[:, :2])
        lo, hi = allp.min(0), allp.max(0)
        return dict(plan_area=hull.volume, bbox=(hi - lo), bbox_vol=float(np.prod(hi - lo)),
                    lo=lo, hi=hi)


def deployed_envelope(lay):
    loc, plate = lay.local_parts()
    zr = np.concatenate([loc['blade'][:, 2], plate[:, 2]])
    h = zr.max() - zr.min()
    D = 2 * P.R_TIP
    return dict(plan_area=np.pi * P.R_TIP ** 2, bbox=np.array([D, D, h]), cyl_vol=np.pi * P.R_TIP ** 2 * h,
                bbox_vol=D * D * h)


def sweep_ok(lay, n=12):
    worst = 1.0
    for f in np.linspace(0, 1, n + 1)[1:]:
        worst = min(worst, min(lay.clearances(f).values()))
    return worst


def optimise(r_list=np.arange(0.16, 0.215, 0.01), s_list=np.arange(0.06, 0.111, 0.01),
             y_list=np.arange(-0.03, 0.0101, 0.01), d_list=np.arange(0.0, 6.01, 1.0), **fixed):
    rows, feas = [], []
    for r_h, s, y_h, dl in itertools.product(r_list, s_list, y_list, d_list):
        lay = Layout(r_h=r_h, s=s, delta_deg=dl, y_h=y_h, **fixed)
        if P.R_ROOT - r_h < lay.R_p + P.CLEARANCE_MIN + 0.02:
            continue
        if r_h - lay.a_t < lay.R_core + P.CLEARANCE_MIN:
            continue
        c = lay.clearances(pairs=(('A', 'B'), ('A', 'C')))
        cmin = min(c.values())
        area = lay.envelope()['plan_area'] if cmin >= P.CLEARANCE_MIN else np.nan
        rows.append((r_h, s, y_h, dl, cmin, area))
        if cmin >= P.CLEARANCE_MIN:
            feas.append((area, lay))
    feas.sort(key=lambda t: t[0])
    for area, lay in feas[:40]:
        w = sweep_ok(lay)
        if w >= P.CLEARANCE_MIN * 0.99:
            return (lay, lay.envelope(), lay.clearances(), w), np.array(rows)
    raise RuntimeError('no feasible layout')


if __name__ == '__main__':
    import time; t0 = time.time()
    (lay, env, c, wsweep), rows = optimise()
    dep = deployed_envelope(lay)
    print('time %.0f s, evaluated %d' % (time.time() - t0, len(rows)))
    print('best r_h %.3f s %.3f y_h %.3f delta %.1f, sweep min clearance %.1f mm'
          % (lay.r_h, lay.s, lay.y_h, lay.delta_deg, wsweep * 1000))
    print('clearances', {k: round(v * 1000, 1) for k, v in c.items()})
    print('folded plan area %.3f m2, bbox %s, bbox vol %.4f' % (env['plan_area'], env['bbox'], env['bbox_vol']))
    print('deployed plan area %.3f m2, bbox %s, bbox vol %.4f' % (dep['plan_area'], dep['bbox'], dep['bbox_vol']))
    np.save(os.path.join(HERE, 'layout_grid.npy'), rows)
    json.dump(dict(r_h=lay.r_h, s=lay.s, y_h=lay.y_h, delta_deg=lay.delta_deg, sweep_min=wsweep,
                   clear={k: float(v) for k, v in c.items()},
                   folded=dict(plan_area=env['plan_area'], bbox=env['bbox'].tolist(), bbox_vol=env['bbox_vol']),
                   deployed=dict(plan_area=dep['plan_area'], bbox=dep['bbox'].tolist(), bbox_vol=dep['bbox_vol'])),
              open(os.path.join(HERE, 'layout_opt.json'), 'w'), indent=1)
