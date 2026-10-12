import os, struct, sys, zlib

MARK = b'\x14\x00\x06\x00\x08\x00'


def unpack(path, out):
    os.makedirs(out, exist_ok=True)
    d = open(path, 'rb').read()
    i, n, k = 0, 0, 0
    while True:
        i = d.find(MARK, i)
        if i < 0:
            break
        _, _, csize, _, nlen = struct.unpack('<IIIII', d[i + 6:i + 26])
        if nlen > 300:
            i += 1
            continue
        name = bytes(((c >> 6) | (c << 2)) & 0xFF for c in d[i + 26:i + 26 + nlen]).decode('latin-1')
        st = i + 26 + nlen
        try:
            raw = zlib.decompress(d[st:st + csize], -15)
        except zlib.error:
            raw = None
        if raw:
            safe = ''.join(c if c.isalnum() or c in '._-' else '_' for c in name)[:60]
            open(os.path.join(out, '%03d_%s' % (n, safe)), 'wb').write(raw)
            k += 1
        n += 1
        i += 1
    return k


if __name__ == '__main__':
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join('..', 'blade_single.SLDPRT')
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join('..', '_extract')
    print(unpack(src, dst))
