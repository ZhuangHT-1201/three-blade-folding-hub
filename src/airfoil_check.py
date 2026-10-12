import os
import numpy as np
from blade_geom import T, section

for x0 in (0.3, 0.84, 1.38):
    s = section(T, x0)
    pts = s.reshape(-1, 2)
    D = np.linalg.norm(pts[:, None] - pts[None], axis=2)
    i, j = np.unravel_index(D.argmax(), D.shape)
    a, b = pts[i], pts[j]
    c = D[i, j]
    u = (b - a) / c
    n = np.array([-u[1], u[0]])
    xi = (pts - a) @ u / c
    eta = (pts - a) @ n / c
    near_a = eta[xi < 0.05].ptp() if hasattr(eta, 'ptp') else np.ptp(eta[xi < 0.05])
    near_b = np.ptp(eta[xi > 0.95])
    print(f'x={x0} A={a*1000} B={b*1000} chord={c*1000:.1f}  thick near A {near_a:.4f}, near B {near_b:.4f}')
    for xc in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 0.9):
        sel = np.abs(xi - xc) < 0.01
        if sel.any():
            e = eta[sel]
            print(f'   x/c={xc:.2f} upper {e.max():+.4f} lower {e.min():+.4f} camber {(e.max()+e.min())/2:+.4f} t {e.max()-e.min():.4f}')
