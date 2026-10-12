import glob, json, os
import numpy as np
from scipy.spatial import cKDTree
from render import read_stl

HERE = os.path.dirname(os.path.abspath(__file__))


def sample(tris, h=1.5):
    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
    n = np.maximum(1, np.ceil(area / (h * h)).astype(int))
    idx = np.repeat(np.arange(len(tris)), n)
    r1, r2 = np.random.rand(len(idx)), np.random.rand(len(idx))
    s = np.sqrt(r1)
    p = (1 - s)[:, None] * a[idx] + (s * (1 - r2))[:, None] * b[idx] + (s * r2)[:, None] * c[idx]
    return np.concatenate([p, tris.reshape(-1, 3)])


FC = json.load(open(os.path.join(HERE, 'fc_report.json')))
G = FC['geom']
ARMS = {'A': 180.0, 'B': 60.0, 'C': 300.0}


def own_holes(n):
    a = np.radians(ARMS[n]); R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    return np.array([R @ np.array(G[k]) for k in ('PIN1', 'PIN2', 'PARK1', 'PARK2')])


def load(tag):
    sets = {'A': [], 'B': [], 'C': [], 'hub': []}
    for f in glob.glob(os.path.join(HERE, 'stl', tag + '_*.stl')):
        name = os.path.basename(f)[len(tag) + 1:-4]
        if name.startswith('hub_'):
            sets['hub'].append(read_stl(f))
        elif name.split('_')[1] in ('blade', 'grip', 'sleeve1', 'sleeve2'):
            sets[name[0]].append(read_stl(f))
    return {k: sample(np.concatenate(v)) for k, v in sets.items()}


def outside_own_clevis(pts, n, R=60.0):
    H = own_holes(n)
    dmin = np.min(np.linalg.norm(pts[:, None, :2] - H[None], axis=2), axis=1)
    return pts[dmin > R]


np.random.seed(0)
res = {}
for tag, frac in (('f25', 0.25), ('mid', 0.5), ('f75', 0.75), ('fold', 1.0)):
    S = load(tag)
    T = {k: cKDTree(v) for k, v in S.items()}
    r = {}
    for a, b in (('A', 'B'), ('A', 'C'), ('B', 'C'), ('A', 'hub'), ('B', 'hub'), ('C', 'hub')):
        if b == 'hub':
            pts = outside_own_clevis(S[a], a)
            d = T['hub'].query(pts, distance_upper_bound=200)[0]
            r[a + '-hub'] = float(d.min())
        else:
            r[a + b] = float(T[a].query(S[b], distance_upper_bound=500)[0].min())
    res[str(frac)] = r
    print(tag, {k: round(v, 1) for k, v in r.items()})
json.dump(res, open(os.path.join(HERE, 'clearance.json'), 'w'), indent=1)
