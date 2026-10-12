import glob, os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
f = glob.glob(os.path.join(here, '..', '_extract', '018*'))[0]
d = open(f, 'rb').read()

lo = np.array([0.00028, -5.41e-05, -3.29e-05]) - 1e-7
hi = np.array([0.0014, 0.000161, 9.09e-05]) + 1e-7

best = None
for off in range(4):
    n = (len(d) - off) // 4
    a = np.frombuffer(d, dtype='<f4', count=n, offset=off).astype(np.float64)
    a = np.nan_to_num(a, nan=9, posinf=9, neginf=9)
    for ph in range(3):
        m = (len(a) - ph) // 3
        t = a[ph:ph + 3 * m].reshape(m, 3)
        ok = np.all((t >= lo) & (t <= hi), axis=1) & (np.abs(t[:, 0]) > 1e-6)
        idx = np.flatnonzero(ok)
        print(off, ph, len(idx), (idx[:3] * 3 + ph) * 4 + off if len(idx) else None)
        if best is None or len(idx) > len(best[0]):
            best = (idx, t[idx], off, ph)

idx, P, off, ph = best
byte = (idx * 3 + ph) * 4 + off
br = np.flatnonzero(np.diff(idx) != 1)
starts = np.concatenate([[0], br + 1]); ends = np.concatenate([br + 1, [len(idx)]])
for s, e in zip(starts, ends):
    if e - s > 50:
        print('block byte', byte[s], 'n', e - s)
U = np.unique(np.round(P, 9), axis=0)
print('unique pts', len(U), 'min', U.min(0), 'max', U.max(0))
np.save(os.path.join(here, 'pts_raw.npy'), P)
np.save(os.path.join(here, 'pts_bytes.npy'), byte)
