import json, math, os, sys, time
import numpy as np
import FreeCAD as App
import Part, Mesh, Import

V = App.Vector
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in dir() else os.getcwd()
ROOT = os.path.dirname(HERE)
CAD = os.path.join(ROOT, 'CAD')
os.makedirs(CAD, exist_ok=True)
os.makedirs(os.path.join(HERE, 'stl'), exist_ok=True)
LOG = open(os.path.join(HERE, 'fc_build_log.txt'), 'w', encoding='utf-8')


def log(*a):
    LOG.write(' '.join(str(x) for x in a) + '\n'); LOG.flush()


D = json.load(open(os.path.join(HERE, 'design.json'), encoding='utf-8'))
mm = lambda k: D[k] * 1000.0
r_h, s, y_h, z_h = mm('r_h'), mm('s'), mm('y_h'), mm('z_h')
d, t_t, R_t, a_t, t_p, R_p, b_a = mm('d'), mm('t_t'), mm('R_t'), mm('a_t'), mm('t_p'), mm('R_p'), mm('b_a')
gap, t_cuff, L_cuff, t_bh = mm('gap'), mm('t_cuff'), mm('L_cuff'), mm('t_bulkhead')
R_core = mm('R_core')
phi = 120.0 + D['delta_deg']
R_ROOT, R_TIP = 280.0, 1400.0
SLEEVE_T = 1.5
RHO = {'plate': 2810, 'grip': 2810, 'spacer': 2810, 'pin': 7800, 'sleeve': 8800, 'blade': None}

z_t0, z_t1 = z_h - t_t / 2, z_h + t_t / 2
z_pu0, z_pu1 = z_t1 + gap, z_t1 + gap + t_p
z_pl0, z_pl1 = z_t0 - gap - t_p, z_t0 - gap
PIN1 = np.array([r_h, y_h + s / 2])
PIN2 = np.array([r_h, y_h - s / 2])


def rot2(p, c, ang):
    a = math.radians(ang)
    R = np.array([[math.cos(a), -math.sin(a)], [math.sin(a), math.cos(a)]])
    return c + R @ (np.asarray(p) - c)


PARK1 = rot2(PIN2, PIN1, +phi)
PARK2 = rot2(PIN1, PIN2, -phi)
log('pins', PIN1, PIN2, 'park', PARK1, PARK2)


def hull2d(P):
    P = sorted(map(tuple, P))
    def cross(o, a, b): return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in P:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(P):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def circle_pts(c, r, n=48):
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([c[0] + r * np.cos(a), c[1] + r * np.sin(a)], 1)


def face_xy(poly, z=0.0):
    pts = [V(float(p[0]), float(p[1]), z) for p in poly]
    return Part.Face(Part.makePolygon(pts + [pts[0]]))


def prism_xy(poly, z0, z1):
    return face_xy(poly, z0).extrude(V(0, 0, z1 - z0))


def cyl(c, r, z0, z1):
    return Part.makeCylinder(r, z1 - z0, V(float(c[0]), float(c[1]), z0), V(0, 0, 1))


def naca0012(c, theta_deg, n=90, offset=0.0):
    beta = np.linspace(0, np.pi, n)
    xi = 0.5 * (1 - np.cos(beta))
    yt = 0.6 * (0.2969 * np.sqrt(xi) - 0.1260 * xi - 0.3516 * xi ** 2 + 0.2843 * xi ** 3 - 0.1036 * xi ** 4)
    up = np.stack([xi, yt], 1); lo = np.stack([xi[::-1], -yt[::-1]], 1)
    P = np.concatenate([up, lo[1:-1]]) * c
    th = math.radians(theta_deg)
    ch = np.array([-math.cos(th), -math.sin(th)]); nr = np.array([-math.sin(th), math.cos(th)])
    le = 0.25 * c * np.array([math.cos(th), math.sin(th)])
    Q = le + P[:, :1] * ch + P[:, 1:] * nr
    if offset:
        nxt = np.roll(Q, -1, 0); prv = np.roll(Q, 1, 0)
        tng = nxt - prv; nrm = np.stack([tng[:, 1], -tng[:, 0]], 1)
        nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
        cen = Q.mean(0)
        sgn = np.sign(((Q - cen) * nrm).sum(1)); sgn[sgn == 0] = 1
        Q = Q + offset * nrm * sgn[:, None]
    return Q


def chord_at(r): return 246.0 + (102.9 - 246.0) * (r - 280.5) / (1399.5 - 280.5)
def pitch_at(r): return 29.30 + (6.80 - 29.30) * (r - 280.5) / (1399.5 - 280.5)


def section_wire(r, offset):
    Q = naca0012(chord_at(r), pitch_at(r), offset=offset)
    pts = [V(r, float(p[0]), float(p[1])) for p in Q]
    return Part.makePolygon(pts + [pts[0]])


def rotz(shape, ang, c=(0, 0)):
    sh = shape.copy()
    sh.rotate(V(c[0], c[1], 0), V(0, 0, 1), ang)
    return sh


brep = os.path.join(HERE, 'blade_design_frame.brep')
if os.path.exists(brep):
    BLADE = Part.Shape(); BLADE.read(brep)
else:
    m = Mesh.Mesh(os.path.join(ROOT, 'blade_rebuilt_mm.stl'))
    sh = Part.Shape(); sh.makeShapeFromMesh(m.Topology, 0.05)
    BLADE = Part.Solid(Part.Shell(sh.Faces))
    BLADE.rotate(V(0, 0, 0), V(1, 0, 0), 180)
    BLADE.exportBrep(brep)
log('blade volume cm3', BLADE.Volume / 1e3, BLADE.BoundBox)

def arm_outline():
    P = [circle_pts(c, R_p) for c in (PIN1, PIN2, PARK1, PARK2)]
    P.append(np.array([[R_core * 0.8, y_h - b_a / 2], [R_core * 0.8, y_h + b_a / 2]]))
    return hull2d(np.concatenate(P))


ARM = arm_outline()
U_PLATE_MAX = max(p[0] for p in ARM)


def plate(z0, z1, lighten=True):
    sh = cyl((0, 0), R_core + 12, z0, z1)
    for psi in (60, 180, 300):
        sh = sh.fuse(rotz(prism_xy(ARM, z0, z1), psi))
    sh = sh.removeSplitter()
    sh = max(sh.Solids, key=lambda x: x.Volume)
    cut = [cyl((0, 0), 16, z0 - 1, z1 + 1)]
    for k in range(6):
        a = math.radians(30 + 60 * k)
        cut.append(cyl((45 * math.cos(a), 45 * math.sin(a)), 3.3, z0 - 1, z1 + 1))
    for psi in (60, 180, 300):
        for c in (PIN1, PIN2, PARK1, PARK2):
            cut.append(rotz(cyl(c, d / 2, z0 - 1, z1 + 1), psi))
        if lighten and LIGHT:
            for (cu, cv, rr) in LIGHT:
                cut.append(rotz(cyl((cu, cv), rr, z0 - 1, z1 + 1), psi))
    for c in cut:
        sh = sh.cut(c)
    return max(sh.removeSplitter().Solids, key=lambda x: x.Volume)


LIGHT = [tuple(x) for x in D.get('lightening_holes', [(118.0, y_h, 15.0)])]
PLATE_U = plate(z_pu0, z_pu1)
PLATE_L = plate(z_pl0, z_pl1)
SPACER = cyl((0, 0), 53.0, z_pl1, z_pu0).cut(cyl((0, 0), 37.0, z_pl1 - 1, z_pu0 + 1))
for k in range(6):
    a = math.radians(30 + 60 * k)
    SPACER = SPACER.cut(cyl((45 * math.cos(a), 45 * math.sin(a)), 3.3, z_pl1 - 1, z_pu0 + 1))
MOTOR = cyl((0, 0), R_core, z_pl0 - 100, z_pl0)

vmin, vmax = y_h - s / 2 - R_t, y_h + s / 2 + R_t
U_S = r_h + s + R_p + 8.0
U_S = max(U_S, U_PLATE_MAX + 8.0)
log('vmin vmax', vmin, vmax, 'U_S', U_S, 'plate max u', U_PLATE_MAX)

tong = prism_xy([[r_h, vmin], [U_S + 2, vmin], [U_S + 2, vmax], [r_h, vmax]], z_t0, z_t1)
tong = tong.fuse(prism_xy([[r_h - a_t, y_h - s / 2], [r_h, y_h - s / 2], [r_h, y_h + s / 2], [r_h - a_t, y_h + s / 2]], z_t0, z_t1))
for c in (PIN1, PIN2):
    tong = tong.fuse(cyl(c, a_t, z_t0, z_t1).common(prism_xy([[r_h - a_t - 1, vmin], [r_h + 1, vmin], [r_h + 1, vmax], [r_h - a_t - 1, vmax]], z_t0 - 1, z_t1 + 1)))

Qr = naca0012(chord_at(R_ROOT), pitch_at(R_ROOT))
low = Qr[len(Qr) // 2:]
order = np.argsort(low[:, 0]); ly, lz = low[order, 0], low[order, 1]
vs = np.linspace(vmin, vmax, 25)
ztop = np.interp(vs, ly, lz) + 2.0
sad2d = [(vmin, z_t0)] + [(vmax, z_t0)] + [(v, z) for v, z in zip(vs[::-1], ztop[::-1])]
sad_face = Part.Face(Part.makePolygon([V(U_S, float(a), float(b)) for a, b in sad2d] + [V(U_S, float(sad2d[0][0]), float(sad2d[0][1]))]))
SADDLE = sad_face.extrude(V(R_ROOT + L_cuff - U_S, 0, 0))

cuff_outer = Part.makeLoft([section_wire(R_ROOT - t_bh, t_cuff), section_wire(R_ROOT + L_cuff, t_cuff)], True, True)
GRIP = tong.fuse(SADDLE).fuse(cuff_outer).removeSplitter()
GRIP = GRIP.cut(BLADE)
for c in (PIN1, PIN2):
    GRIP = GRIP.cut(cyl(c, d / 2 + SLEEVE_T, z_t0 - 1, z_t1 + 1))
POCKET_DEPTH = 16.0
GRIP = GRIP.cut(Part.makeBox(R_ROOT + L_cuff - 6 - (U_S + 6), (vmax - 5) - (vmin + 5), POCKET_DEPTH + 1,
                             V(U_S + 6, vmin + 5, z_t0 - 1)))
GRIP = GRIP.removeSplitter()
GRIP = max(GRIP.Solids, key=lambda x: x.Volume)
log('grip solids', len(GRIP.Solids), 'vol cm3', GRIP.Volume / 1e3, 'valid', GRIP.isValid())

SLEEVE = cyl((0, 0), d / 2 + SLEEVE_T, z_t0 - gap, z_t1 + gap).cut(cyl((0, 0), d / 2, z_t0 - 2, z_t1 + 2))


def pivot_bolt():
    shank = cyl((0, 0), d / 2, z_pl0 - 0.0, z_pu1 + 0.0)
    head = Part.makeCylinder(0.9 * d, 0.65 * d, V(0, 0, z_pu1), V(0, 0, 1))
    nut = Part.makeCylinder(0.9 * d, 0.8 * d, V(0, 0, z_pl0 - 0.8 * d), V(0, 0, 1))
    return shank.fuse(head).fuse(nut).removeSplitter()


def lock_pin():
    shank = cyl((0, 0), d / 2, z_pl0 - 0.5 * d, z_pu1)
    knob = Part.makeCylinder(1.1 * d, 0.6 * d, V(0, 0, z_pu1), V(0, 0, 1))
    stem = Part.makeCylinder(0.35 * d, 1.2 * d, V(0, 0, z_pu1 + 0.6 * d), V(0, 0, 1))
    ring = Part.makeTorus(1.4 * d, 0.18 * d, V(0, 0, z_pu1 + 2.2 * d), V(1, 0, 0))
    return shank.fuse(knob).fuse(stem).fuse(ring).removeSplitter()


BOLT, LPIN = pivot_bolt(), lock_pin()

ARMS = {'A': 180.0, 'B': 60.0, 'C': 300.0}


def blade_set(name, frac):
    psi = ARMS[name]
    parts = {}
    loc = {'blade': BLADE, 'grip': GRIP,
           'sleeve1': SLEEVE.copy().translate(V(PIN1[0], PIN1[1], 0)) if False else SLEEVE.translated(V(PIN1[0], PIN1[1], 0)),
           'sleeve2': SLEEVE.translated(V(PIN2[0], PIN2[1], 0))}
    piv, lock_c = None, None
    if name == 'B':
        piv, ang = PIN1, +phi * frac
    elif name == 'C':
        piv, ang = PIN2, -phi * frac
    for k, sh in loc.items():
        x = sh.copy()
        if piv is not None and frac > 0:
            x.rotate(V(piv[0], piv[1], 0), V(0, 0, 1), ang)
        x.rotate(V(0, 0, 0), V(0, 0, 1), psi)
        parts[k] = x
    if name == 'A' or frac == 0:
        pv, lk = PIN1, PIN2
    elif name == 'B':
        pv, lk = PIN1, (PARK1 if frac >= 1 else None)
    else:
        pv, lk = PIN2, (PARK2 if frac >= 1 else None)
    parts['pivot_bolt'] = rotz(BOLT.translated(V(pv[0], pv[1], 0)), psi)
    if lk is not None:
        parts['lock_pin'] = rotz(LPIN.translated(V(lk[0], lk[1], 0)), psi)
    return parts


def hub_parts():
    return {'plate_upper': PLATE_U, 'plate_lower': PLATE_L, 'spacer': SPACER, 'motor_ref': MOTOR}


def build_doc(frac, tag, doc):
    grp = doc.addObject('App::DocumentObjectGroup', tag)
    objs = []
    for k, sh in hub_parts().items():
        o = doc.addObject('Part::Feature', f'{tag}_{k}'); o.Shape = sh; grp.addObject(o); objs.append(o)
    for n in ARMS:
        for k, sh in blade_set(n, frac).items():
            o = doc.addObject('Part::Feature', f'{tag}_{n}_{k}'); o.Shape = sh; grp.addObject(o); objs.append(o)
    return objs


def interference(frac):
    hub = hub_parts()
    sets = {n: blade_set(n, frac) for n in ARMS}
    res = {}
    def body(n, keys=('blade', 'grip', 'sleeve1', 'sleeve2')):
        return [sets[n][k] for k in keys if k in sets[n]]
    pairs = [('A', 'B'), ('A', 'C'), ('B', 'C')]
    for a, b in pairs:
        dist = min(x.distToShape(y)[0] for x in body(a) for y in body(b))
        vol = sum(x.common(y).Volume for x in body(a) for y in body(b) if x.BoundBox.intersect(y.BoundBox))
        res[a + b] = dict(dist=dist, vol=vol)
    for n in ARMS:
        dist, vol = 1e9, 0.0
        for x in body(n, ('blade', 'grip')):
            for k, y in hub.items():
                if x.BoundBox.intersect(y.BoundBox):
                    vol += x.common(y).Volume
                dist = min(dist, x.distToShape(y)[0])
        res[n + '-hub'] = dict(dist=dist, vol=vol)
    return res


report = {}
t0 = time.time()
for k, sh in dict(PLATE_U=PLATE_U, SPACER=SPACER, GRIP=GRIP, SLEEVE=SLEEVE, BOLT=BOLT, LPIN=LPIN).items():
    log(k, 'valid', sh.isValid(), 'vol', sh.Volume)
masses = dict(plate=PLATE_U.Volume * 1e-9 * RHO['plate'], spacer=SPACER.Volume * 1e-9 * RHO['spacer'],
              grip=GRIP.Volume * 1e-9 * RHO['grip'], sleeve=SLEEVE.Volume * 1e-9 * RHO['sleeve'],
              bolt=BOLT.Volume * 1e-9 * RHO['pin'], lock_pin=LPIN.Volume * 1e-9 * RHO['pin'])
masses['mechanism_total'] = 2 * masses['plate'] + masses['spacer'] + 3 * (masses['grip'] + 2 * masses['sleeve'] + masses['bolt'] + masses['lock_pin'])
gc = GRIP.CenterOfMass
masses['grip_rcg_m'] = gc.x / 1000.0
report['mass'] = masses
log('masses', masses)

report['grip_overlap'] = {}
for frac in (0.25, 0.5, 0.75, 1.0):
    sets = {n: blade_set(n, frac) for n in ARMS}
    vol = {}
    for a, b in (('A', 'B'), ('A', 'C'), ('B', 'C')):
        vol[a + b] = sets[a]['grip'].common(sets[b]['grip']).Volume if sets[a]['grip'].BoundBox.intersect(sets[b]['grip'].BoundBox) else 0.0
    for n in ARMS:
        v = 0.0
        for k, y in hub_parts().items():
            if sets[n]['grip'].BoundBox.intersect(y.BoundBox):
                v += sets[n]['grip'].common(y).Volume
        vol[n + '-hub'] = v
    report['grip_overlap'][str(frac)] = vol
    log('grip overlap frac', frac, vol, '%.0fs' % (time.time() - t0))

def env(frac):
    bb = App.BoundBox()
    for sh in list(hub_parts().values())[:3]:
        bb.add(sh.BoundBox)
    for n in ARMS:
        for k, sh in blade_set(n, frac).items():
            bb.add(sh.BoundBox)
    return [bb.XMin, bb.XMax, bb.YMin, bb.YMax, bb.ZMin, bb.ZMax]
report['bbox_folded'] = env(1.0)
report['bbox_deployed'] = env(0.0)
log('bbox folded', report['bbox_folded'], 'deployed', report['bbox_deployed'])

doc = App.newDocument('FoldingHub')
dep = build_doc(0.0, 'Deployed', doc)
fol = build_doc(1.0, 'Folded', doc)
doc.recompute()
doc.saveAs(os.path.join(CAD, 'folding_hub.FCStd'))
Import.export(dep, os.path.join(CAD, 'assembly_deployed.step'))
Import.export(fol, os.path.join(CAD, 'assembly_folded.step'))
for k, sh in dict(hub_plate_upper=PLATE_U, hub_plate_lower=PLATE_L, hub_spacer=SPACER, blade_grip=GRIP,
                  sleeve_bushing=SLEEVE, pivot_bolt=BOLT, lock_pin=LPIN).items():
    sh.exportStep(os.path.join(CAD, k + '.step'))
BLADE.exportStep(os.path.join(CAD, 'blade_design_frame.step'))

def stl(sh, name, tol=0.3):
    sh.exportStl(os.path.join(HERE, 'stl', name + '.stl'))

FRACS = ((0.0, 'dep'), (1.0, 'fold'), (0.5, 'mid'), (0.25, 'f25'), (0.75, 'f75'))
for frac, tag in FRACS:
    for k, sh in hub_parts().items():
        stl(sh, f'{tag}_hub_{k}')
    for n in ARMS:
        for k, sh in blade_set(n, frac).items():
            stl(sh, f'{tag}_{n}_{k}')
report['blade_transforms'] = {}
for frac, tag in FRACS:
    for n in ARMS:
        b = blade_set(n, frac)['blade']
        report['blade_transforms'][f'{tag}_{n}'] = [list(r) for r in np.array(b.Placement.Matrix.A).reshape(4, 4)]
for k, sh in dict(grip_local=GRIP, plate_upper_local=PLATE_U).items():
    stl(sh, k)
report['geom'] = dict(U_S=U_S, U_PLATE_MAX=U_PLATE_MAX, PIN1=PIN1.tolist(), PIN2=PIN2.tolist(),
                      PARK1=PARK1.tolist(), PARK2=PARK2.tolist(), z=dict(t0=z_t0, t1=z_t1, pu0=z_pu0, pu1=z_pu1,
                                                                            pl0=z_pl0, pl1=z_pl1), ARM=[list(p) for p in ARM])
json.dump(report, open(os.path.join(HERE, 'fc_report.json'), 'w'), indent=1)
log('done %.0fs' % (time.time() - t0))
