# CMS R1 item 6 - band paths: band-edge table for all strains + E_g = min(6x6x4 mesh, path) + band-structure figure (-5, 0, +5 %)
# Input: R1_bands.tar.gz (bands.in, bandsx_s?.out, bands_s?.dat.gnu of strain_*/FM/R1_bands). Output: band_edges_R1.dat (revision1/data/), fig_bands (Fig. 6)
# Written as a Google Colab cell; runs locally with BASE=<folder>. Column headers: yol = path, izgara = mesh.
import os, re, glob, tarfile, numpy as np, matplotlib.pyplot as plt
try:
    from google.colab import drive; drive.mount('/content/drive'); BASE = '/content/drive/MyDrive/Research/DFT/LNMO_2'
except ImportError: BASE = os.environ.get('BASE', '.')
TGZ = [p for p in (f'{BASE}/R1/R1_bands.tar.gz', f'{BASE}/R1_bands.tar.gz') if os.path.exists(p)][0]
OUT = f'{BASE}/R1'; os.makedirs(OUT, exist_ok=True); X = '/content/r1bands' if os.path.isdir('/content') else './r1bands'
tarfile.open(TGZ).extractall(X); print('input: %s (%.0f kB)' % (TGZ, os.path.getsize(TGZ) / 1024))
assert glob.glob(f'{X}/strain_*/FM/R1_bands/*.gnu'), 'no .gnu files in the archive'
GU = dict(zip(range(-5, 6), [1.0763, 1.1430, 1.1942, 1.2406, 1.2764, 1.2950, 1.2953, 1.2526, 1.1912, 1.1325, 1.0836]))   # data/gap.dat (mesh)
GD = dict(zip(range(-5, 6), [3.1580, 3.2022, 3.2304, 3.2437, 3.2439, 3.2321, 3.1812, 3.1164, 3.0474, 2.9753, 2.8969]))
LAB = ['Γ', 'X', 'S', 'Y', 'Γ', 'Z', 'U', 'R', 'T', 'Z']                 # X(½,0,0) S(½,½,0) Y(0,½,0) Z(0,0,½) U(½,0,½) R(½,½,½) T(0,½,½)
def gnu(f):
    B = [np.array([list(map(float, l.split())) for l in b.splitlines() if l.strip()]) for b in re.split(r'\n\s*\n', open(f).read()) if b.strip()]
    return B[0][:, 0], np.array([b[:, 1] for b in B])
def where(i, ih): return LAB[ih.index(i)] if i in ih else '%s–%s' % (LAB[i // 20], LAB[i // 20 + 1])
D, rows = {}, []
for d in sorted(glob.glob(f'{X}/strain_*/FM/R1_bands'), key=lambda s: float(re.search(r'strain_(-?[\d.]+)', s).group(1))):
    e = round(100 * float(re.search(r'strain_(-?[\d.]+)', d).group(1))); k, up = gnu(f'{d}/bands_s1.dat.gnu'); _, dn = gnu(f'{d}/bands_s2.dat.gnu')
    hs = [float(x) for x in re.findall(r'x coordinate\s+(\S+)', open(f'{d}/bandsx_s1.out').read())]; ih = [int(np.argmin(abs(k - h))) for h in hs]
    gu = up[88].min() - up[87].max(); gd = dn[78].min() - dn[77].max(); D[e] = (k, up, dn, hs, max(up[87].max(), dn[77].max()))
    rows.append((e, gu, GU[e], min(gu, GU[e]), gd, GD[e], min(gd, GD[e]), where(up[87].argmax(), ih), where(up[88].argmin(), ih), where(dn[77].argmax(), ih), where(dn[78].argmin(), ih)))
H = 'eps  Eg_up(yol) Eg_up(izgara) Eg_up(min)  Eg_dn(yol) Eg_dn(izgara) Eg_dn(min)  VBM_up CBM_up VBM_dn CBM_dn'
print(H); [print('%+d  %.4f %.4f %.4f   %.4f %.4f %.4f   %s %s %s %s' % r) for r in rows]
with open(f'{OUT}/band_edges_R1.dat', 'w') as f: f.write('# ' + H + '\n'); [f.write('%+d %.4f %.4f %.4f %.4f %.4f %.4f %s %s %s %s\n' % r) for r in rows]
fig, ax = plt.subplots(1, 3, figsize=(11, 4), sharey=True)
for a, e in zip(ax, (-5, 0, 5)):
    k, up, dn, hs, ev = D[e]
    for b in up: a.plot(k, b - ev, c='#0072B2', lw=.8)
    for b in dn: a.plot(k, b - ev, c='#D55E00', lw=.8, ls=(0, (3, 2)))
    for b, f, m in ((up[87], np.argmax, 'o'), (up[88], np.argmin, 'o'), (dn[77], np.argmax, 's'), (dn[78], np.argmin, 's')):
        i = f(b); a.plot(k[i], b[i] - ev, m, ms=4, mfc='none', mec='k', mew=.8)
    for x in hs: a.axvline(x, c='gray', lw=.4)
    a.axhline(0, c='k', lw=.5, ls=':'); a.set_xlim(k[0], k[-1]); a.set_ylim(-2, 4); a.set_xticks(hs); a.set_xticklabels(LAB[:len(hs)])
    a.set_title(r'$\varepsilon$ = %+d %%' % e if e else r'$\varepsilon$ = 0 %')
ax[0].set_ylabel('$E - E_{\\mathrm{VBM}}$ (eV)')
fig.legend(handles=[plt.Line2D([], [], c='#0072B2', lw=1, label='majority'), plt.Line2D([], [], c='#D55E00', lw=1, ls=(0, (3, 2)), label='minority')],
           loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.02)); fig.tight_layout(rect=(0, 0, 1, 0.95))
for x in ('png', 'pdf'): fig.savefig(f'{OUT}/fig_bands.{x}', dpi=300)
print('saved:', OUT)
