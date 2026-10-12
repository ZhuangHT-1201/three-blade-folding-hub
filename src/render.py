import glob, json, os, struct, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection

plt.rcParams['axes.unicode_minus'] = False

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(ROOT, 'figures')
os.makedirs(FIG, exist_ok=True)


def read_stl(path):
    data = open(path, 'rb').read()
    if data[:5] == b'solid' and b'facet' in data[:300]:
        v = []
        for line in data.decode('latin1').splitlines():
            t = line.split()
            if t and t[0] == 'vertex':
                v.append([float(t[1]), float(t[2]), float(t[3])])
        return np.array(v).reshape(-1, 3, 3)
    n = struct.unpack('<I', data[80:84])[0]
    a = np.frombuffer(data[84:84 + n * 50], dtype=np.dtype([('n', '<3f4'), ('v', '<9f4'), ('a', '<u2')]))
    return a['v'].reshape(-1, 3, 3).astype(float)


COL = {'blade': (0.80, 0.82, 0.86), 'grip': (0.93, 0.55, 0.20), 'plate': (0.35, 0.55, 0.85),
       'spacer': (0.35, 0.55, 0.85), 'motor': (0.35, 0.35, 0.38), 'pin': (0.85, 0.15, 0.15),
       'sleeve': (0.75, 0.65, 0.25)}


def part_color(name):
    for k, c in COL.items():
        if k in name:
            return c
    if 'bolt' in name or 'lock' in name:
        return COL['pin']
    return (0.6, 0.6, 0.6)


def scene(tag):
    items = []
    for f in glob.glob(os.path.join(HERE, 'stl', tag + '_*.stl')):
        name = os.path.basename(f)[len(tag) + 1:-4]
        items.append((name, read_stl(f)))
    return items


def shade(tris, base, light=np.array([0.3, -0.5, 0.8])):
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    l = light / np.linalg.norm(light)
    k = 0.45 + 0.55 * np.abs(n @ l)
    return np.clip(np.outer(k, base), 0, 1)


def render3d(tag, fname, elev=28, azim=-60, title=None, size=(10, 7)):
    items = scene(tag)
    fig = plt.figure(figsize=size, dpi=150)
    ax = fig.add_subplot(111, projection='3d')
    allp = []
    for name, t in items:
        if 'motor' in name:
            continue
        ax.add_collection3d(Poly3DCollection(t, facecolors=shade(t, part_color(name)), edgecolors='none', linewidths=0))
        allp.append(t.reshape(-1, 3))
    P_ = np.concatenate(allp)
    lo, hi = P_.min(0), P_.max(0)
    c = (lo + hi) / 2; r = (hi - lo).max() / 2
    ax.set_xlim(c[0] - r, c[0] + r); ax.set_ylim(c[1] - r, c[1] + r); ax.set_zlim(c[2] - r, c[2] + r)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel('X [mm]'); ax.set_ylabel('Y [mm]'); ax.set_zlabel('Z [mm]')
    if title: ax.set_title(title)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, fname)); plt.close(fig)


def render_zoom(tag, fname, title, R=260, elev=35, azim=-50):
    items = scene(tag)
    fig = plt.figure(figsize=(9, 7), dpi=150)
    ax = fig.add_subplot(111, projection='3d')
    for name, t in items:
        if 'motor' in name:
            continue
        m = (np.abs(t[..., 0]).max(1) < R * 1.3) & (np.abs(t[..., 1]).max(1) < R * 1.3)
        t = t[m]
        if len(t) == 0:
            continue
        ax.add_collection3d(Poly3DCollection(t, facecolors=shade(t, part_color(name)), edgecolors='none'))
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(-R, R)
    ax.set_box_aspect((1, 1, 1)); ax.view_init(elev=elev, azim=azim)
    ax.set_title(title)
    ax.set_xlabel('X [mm]'); ax.set_ylabel('Y [mm]'); ax.set_zlabel('Z [mm]')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fname)); plt.close(fig)


def topview(tags, fname, titles):
    fig, axs = plt.subplots(1, len(tags), figsize=(5.2 * len(tags), 5.4), dpi=150)
    for ax, tag, tt in zip(np.atleast_1d(axs), tags, titles):
        items = scene(tag)
        items.sort(key=lambda it: np.mean(it[1][..., 2]))
        for name, t in items:
            if 'motor' in name:
                continue
            ax.add_collection(PolyCollection(t[..., :2], facecolors=part_color(name), edgecolors='none', alpha=0.9))
        ax.add_patch(plt.Circle((0, 0), 1400, fill=False, ls='--', color='gray', lw=0.8))
        ax.set_xlim(-1500, 1500); ax.set_ylim(-1500, 1500); ax.set_aspect('equal')
        ax.set_title(tt); ax.set_xlabel('X [mm]'); ax.set_ylabel('Y [mm]'); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fname)); plt.close(fig)


def hub_topview(fname):
    rep = json.load(open(os.path.join(HERE, 'fc_report.json')))
    fig, axs = plt.subplots(1, 2, figsize=(12, 6), dpi=150)
    for ax, tag, tt in zip(axs, ('dep', 'fold'), ('Deployed - hub region', 'Folded - hub region')):
        items = scene(tag)
        items.sort(key=lambda it: np.mean(it[1][..., 2]))
        for name, t in items:
            if 'motor' in name:
                continue
            m = (np.abs(t[..., 0]).max(1) < 450) & (np.abs(t[..., 1]).max(1) < 450)
            ax.add_collection(PolyCollection(t[m][..., :2], facecolors=part_color(name), edgecolors='none', alpha=0.85))
        g = rep['geom']
        for psi in (60, 180, 300):
            a = math.radians(psi); Rz = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
            for p in (g['PIN1'], g['PIN2'], g['PARK1'], g['PARK2']):
                q = Rz @ np.array(p)
                ax.plot(q[0], q[1], 'k+', ms=6)
        ax.set_xlim(-420, 420); ax.set_ylim(-420, 420); ax.set_aspect('equal'); ax.grid(alpha=0.3)
        ax.set_title(tt); ax.set_xlabel('X [mm]'); ax.set_ylabel('Y [mm]')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, fname)); plt.close(fig)


def load_plots():
    import loads as L, params as P
    b = L.bemt(P.RPM_NOM)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), dpi=150)
    axs[0].plot(L.r, L.c * 1000); ax2 = axs[0].twinx()
    ax2.plot(L.r, np.degrees(L.theta), 'r'); axs[0].set_xlabel('r [m]')
    axs[0].set_ylabel('Chord [mm]'); ax2.set_ylabel('Pitch [deg]', color='r'); axs[0].set_title('Blade geometry (recovered from SLDPRT)')
    axs[1].plot(L.r, b['dT'], label='dT/dr thrust'); axs[1].plot(L.r, b['dD'], label='dD/dr drag')
    axs[1].set_xlabel('r [m]'); axs[1].set_ylabel('N/m'); axs[1].legend(); axs[1].set_title(f'BEMT load distribution @ {P.RPM_NOM:.0f} rpm')
    Om2 = b['Om'] ** 2
    cf = np.array([np.trapezoid((L.MPL * Om2 * L.r)[k:], L.r[k:]) for k in range(len(L.r))])
    axs[2].plot(L.r, L.MPL); ax3 = axs[2].twinx()
    ax3.plot(L.r, cf / 1000, 'r'); axs[2].set_xlabel('r [m]'); axs[2].set_ylabel('Mass per length [kg/m]')
    ax3.set_ylabel('Centrifugal force [kN]', color='r'); axs[2].set_title('Composite blade mass and centrifugal force')
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'loads.png')); plt.close(fig)


def pareto_plot():
    rows = [r for r in json.load(open(os.path.join(HERE, 'pareto.json'))) if r.get('feasible')]
    D = json.load(open(os.path.join(HERE, 'design.json'), encoding='utf-8'))
    fig, ax = plt.subplots(figsize=(7.5, 5.5), dpi=150)
    sc = ax.scatter([r['area'] for r in rows], [r['mass'] for r in rows], c=[r['r_h'] * 1000 for r in rows],
                    s=[(r['s'] * 1000 - 30) * 4 for r in rows], cmap='viridis')
    pk = D['pareto_pick']
    ax.scatter([pk['area']], [pk['mass']], s=220, facecolors='none', edgecolors='r', lw=2, label='Selected design')
    fig.colorbar(sc, label='Hinge radius r_h [mm] (marker size ~ pin spacing s)')
    ax.set_xlabel('Folded plan area [m$^2$]'); ax.set_ylabel('Mechanism mass [kg]'); ax.legend()
    ax.set_title('Hinge layout Pareto front: folded size vs mass'); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, 'pareto.png')); plt.close(fig)


if __name__ == '__main__':
    import sys
    what = sys.argv[1:] or ['all']
    if 'all' in what or 'asm' in what:
        render3d('dep', 'asm_deployed.png', title='Deployed (flight)')
        render3d('mid', 'asm_midfold.png', title='Folding (50%)')
        render3d('fold', 'asm_folded.png', title='Folded: all three blades stowed in one direction')
        render_zoom('dep', 'hub_deployed_3d.png', 'Folding hub (deployed)')
        render_zoom('fold', 'hub_folded_3d.png', 'Folding hub (folded)')
        topview(['dep', 'mid', 'fold'], 'fold_sequence_top.png', ['Deployed', 'Folding (50%)', 'Folded'])
        hub_topview('hub_top.png')
    if 'all' in what or 'plots' in what:
        load_plots(); pareto_plot()
