#!/usr/bin/env python3
# gap_aaa.py - a-a-a- yapisi E_g: nscf_keep.out (6x6x4 izgara) + bands_s?.dat.gnu (yol); bant 88 dolu (up), 78 dolu (dw)
import re, sys, glob, numpy as np
R = '/arf/scratch/nparan/R1/tilt'; NO = {1: 88, 2: 78}; F = r'-?\d+\.\d+'
def grid(f):
    t = open(f).read(); t = t[t.find('End of band structure calculation'):]
    up, dw = t.split('SPIN DOWN')[0], t.split('SPIN DOWN')[1]
    out = {}
    for s, blk in ((1, up), (2, dw)):
        ks = re.findall(r'k =\s*(' + F + r')\s*(' + F + r')\s*(' + F + r').*?bands \(ev\):\s*\n(.*?)(?=\n\s*(?:k =|occupation|the Fermi|------|Writing))', blk, re.S)
        E = np.array([[float(x) for x in re.findall(F, b[3])][:100] for b in ks]); K = [b[:3] for b in ks]
        n = NO[s]; v, c = E[:, n-1], E[:, n]
        out[s] = (v.max(), c.min(), K[v.argmax()], K[c.argmin()], len(K))
    return out
def path(f, n):
    B = [np.loadtxt(b.splitlines()) for b in re.split(r'\n\s*\n', open(f).read().strip()) if b.strip()]
    v, c = B[n-1], B[n]; return v[:, 1].max(), c[:, 1].min(), v[v[:, 1].argmax(), 0], c[c[:, 1].argmin(), 0]
LB = 'G X S Y G Z U R T Z'.split()
def lab(x, f):                              # yol konumu -> yuksek simetri etiketi (bands.x 'x coordinate' satirlari)
    try: H = [float(l.split()[-1]) for l in open(f) if 'x coordinate' in l]
    except OSError: return f'x={x:.3f}'
    i = int(np.searchsorted(H, x - 1e-3)); i = min(i, len(H) - 1)
    return LB[i] if abs(H[i] - x) < 2e-3 else f'{LB[i-1]}-{LB[i]}({(x-H[i-1])/(H[i]-H[i-1]):.2f})'
for e in (sys.argv[1:] or ['-0.04', '-0.05']):
    G = f'{R}/eps_{e}/gap'; g = grid(f'{G}/nscf_keep.out')
    for s, sp in ((1, 'up'), (2, 'dw')):
        vb, cb, kv, kc, nk = g[s]; pv, pc, xv, xc = path(f'{G}/bands_s{s}.dat.gnu', NO[s])
        hx = '%s/bandsx_s%d.out' % (G, s)
        print(f'{e} {sp}: izgara({nk} k) Eg={cb-vb:.4f} VBM@{kv} CBM@{kc} | yol Eg={pc-pv:.4f} VBM@{lab(xv, hx)} CBM@{lab(xc, hx)} | Eg=min {min(cb-vb, pc-pv):.4f} eV')
