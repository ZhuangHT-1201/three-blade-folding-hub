import os, struct
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
P = np.load(os.path.join(here, 'pts_raw.npy')) * 1000.0
T = P.reshape(-1, 3, 3)


def write_stl(path, tris):
    n = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-30
    with open(path, 'wb') as f:
        f.write(b'blade_single rebuilt from SLDPRT, units mm'.ljust(80, b' '))
        f.write(struct.pack('<I', len(tris)))
        for k in range(len(tris)):
            f.write(struct.pack('<12fH', *n[k], *(tris[k].ravel()), 0))


def mesh_props(tris):
    a, b, c = tris[:, 0], tris[:, 1], tris[:, 2]
    v6 = np.einsum('ij,ij->i', a, np.cross(b, c))
    V = v6.sum() / 6
    cen = (v6[:, None] * (a + b + c)).sum(0) / (24 * V)
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum()
    return V, cen, area


def section(tris, x0):
    d = tris[:, :, 0] - x0
    segs = []
    for t, s in zip(tris, d):
        pos = s > 0
        if pos.all() or (~pos).all():
            continue
        pts = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if (s[i] > 0) != (s[j] > 0):
                u = s[i] / (s[i] - s[j])
                pts.append(t[i] + u * (t[j] - t[i]))
        if len(pts) == 2:
            segs.append((pts[0][1:], pts[1][1:]))
    return np.array(segs)


def section_props(segs):
    p = segs[:, 0]; q = segs[:, 1]
    cr = p[:, 0] * q[:, 1] - q[:, 0] * p[:, 1]
    A = cr.sum() / 2
    cy = ((p[:, 0] + q[:, 0]) * cr).sum() / (6 * A)
    cz = ((p[:, 1] + q[:, 1]) * cr).sum() / (6 * A)
    Iyy = ((p[:, 1] ** 2 + p[:, 1] * q[:, 1] + q[:, 1] ** 2) * cr).sum() / 12
    Izz = ((p[:, 0] ** 2 + p[:, 0] * q[:, 0] + q[:, 0] ** 2) * cr).sum() / 12
    Iyz = ((p[:, 0] * q[:, 1] + 2 * p[:, 0] * p[:, 1] + 2 * q[:, 0] * q[:, 1] + q[:, 0] * p[:, 1]) * cr).sum() / 24
    s = np.sign(A); A = abs(A)
    Iyy = s * Iyy - A * cz ** 2; Izz = s * Izz - A * cy ** 2; Iyz = s * Iyz - A * cy * cz
    I1 = (Iyy + Izz) / 2 - np.hypot((Iyy - Izz) / 2, Iyz)
    I2 = (Iyy + Izz) / 2 + np.hypot((Iyy - Izz) / 2, Iyz)
    pts = segs.reshape(-1, 2)
    D = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    i, j = np.unravel_index(D.argmax(), D.shape)
    le, te = pts[i], pts[j]
    if le[0] < te[0]:
        le, te = te, le
    chord = D[i, j]
    ang = np.degrees(np.arctan2(le[1] - te[1], le[0] - te[0]))
    cv = (le - te) / chord
    nv = np.array([-cv[1], cv[0]])
    th = (pts @ nv).max() - (pts @ nv).min()
    return dict(A=A, cy=cy, cz=cz, I1=I1, I2=I2, chord=chord, twist=ang, thick=th,
                ymin=pts[:, 0].min(), ymax=pts[:, 0].max(), zmin=pts[:, 1].min(), zmax=pts[:, 1].max())


if __name__ == '__main__':
    write_stl(os.path.join(root, 'blade_rebuilt_mm.stl'), T * 1000)
    V, cen, S = mesh_props(T)
    print(f'triangles {len(T)}, volume {V*1e6:.1f} cm^3, centroid {cen}, surface {S:.4f} m^2')
    xs = np.linspace(0.2805, 1.3995, 25)
    rows = []
    for x in xs:
        sp = section_props(section(T, x))
        rows.append([x, sp['chord'], sp['twist'], sp['thick'], sp['A'], sp['cy'], sp['cz'], sp['I1'], sp['I2']])
        print(f"r={x:.3f} c={sp['chord']*1000:6.1f}mm tw={sp['twist']:6.2f} t={sp['thick']*1000:5.1f}mm t/c={sp['thick']/sp['chord']:.3f} "
              f"A={sp['A']*1e6:7.1f}mm2 cy={sp['cy']*1000:6.1f} cz={sp['cz']*1000:6.1f}")
    np.save(os.path.join(here, 'sections.npy'), np.array(rows))
