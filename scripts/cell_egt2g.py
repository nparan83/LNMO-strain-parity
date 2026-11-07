# CMS R1 item 5 - e_g/t2g PDOS (pseudocubic frame) + e_g occupations. Input: R1_pdos_d.tar.gz (projwfc.x m-resolved d PDOS of Ni1/Mn1)
# Output: eg_t2g_occ.dat (revision1/data/), fig_egt2g.png/pdf (Fig. 5). Written as a Google Colab cell; runs locally with BASE=<folder>.
import os, re, glob, tarfile, numpy as np, matplotlib.pyplot as plt
try:
    from google.colab import drive; drive.mount('/content/drive'); BASE = '/content/drive/MyDrive/Research/DFT/LNMO_2'
except ImportError: BASE = os.environ.get('BASE', '.')
TGZ = [p for p in (f'{BASE}/R1/R1_pdos_d.tar.gz', f'{BASE}/R1_pdos_d.tar.gz') if os.path.exists(p)][0]
OUT = f'{BASE}/R1'; os.makedirs(OUT, exist_ok=True); X = '/content/r1pdos' if os.path.isdir('/content') else './r1pdos'
tarfile.open(TGZ).extractall(X); print('input:', TGZ)
EPS = sorted({float(re.search(r'strain_(-?\d\.\d+)', p).group(1)) for p in glob.glob(f'{X}/strain_*')})
GAP = {-5: 1.0763, -4: 1.1430, -3: 1.1942, -2: 1.2406, -1: 1.2764, 0: 1.2950, 1: 1.2953, 2: 1.2526, 3: 1.1912, 4: 1.1325, 5: 1.0836}  # data/gap.dat (E_g up, mesh)
def eigs(f, nb=100):
    t = open(f).read(); t = t[t.rfind('End of band structure calculation'):]; up, dn = t.split('SPIN DOWN')[:2]
    def blk(s): return np.array([[float(v) for v in re.findall(r'-?\d+\.\d+', re.split(r'occupation numbers|k =|SPIN|Fermi|highest', s[m.end():])[0])][:nb] for m in re.finditer(r'bands \(ev\):', s)])
    return blk(up), blk(dn)
def struct(f):
    L = open(f).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip().upper().startswith('CELL_PARAMETERS')][0]
    A = np.array([[float(x) for x in L[i + j].split()[:3]] for j in (1, 2, 3)]); j = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')][0]
    S = [l.split()[0] for l in L[j + 1:j + 21]]; C = np.array([[float(v) for v in l.split()[1:4]] for l in L[j + 1:j + 21]]) @ A
    return A, S, C
def Q(n):
    M = np.zeros((3, 3)); d = {'z2': (np.diag([-1, -1, 2]) / 6 ** .5), 'x2y2': np.diag([1, -1, 0]) / 2 ** .5}
    if n in d: return d[n]
    a, b = {'xy': (0, 1), 'xz': (0, 2), 'yz': (1, 2)}[n]; M[a, b] = M[b, a] = 2 ** -.5; return M
F = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 2 ** .5]]) / 2 ** .5
def leak(A, S, C, s):                       # e_g -> t2g mixing caused by the deviation of the real octahedron from the common pseudocubic frame
    c = C[S.index(s)]; v = np.array([C[k] + (np.array(n) - 1) @ A - c for k, t in enumerate(S) if t == 'O' for n in np.ndindex(3, 3, 3)])
    v = v[np.argsort(np.linalg.norm(v, axis=1))][:6]; u = v / np.linalg.norm(v, axis=1)[:, None]
    U, _, Vt = np.linalg.svd(np.array([u[np.argmax(u @ ax)] for ax in F])); R = (U @ Vt) @ F.T
    return sum(np.sum((R.T @ Q(e) @ R) * Q(t)) ** 2 for e in ('z2', 'x2y2') for t in ('xy', 'xz', 'yz')) / 2
TZ = getattr(np, 'trapezoid', None) or np.trapz; rows, P = [], {}
for e in EPS:
    d = f'{X}/strain_{e:.2f}/FM'; up, dn = eigs(f'{d}/nscf_FM.out')
    vbm = max(up[:, 87].max(), dn[:, 77].max()); g_up = up[:, 88].min() - up[:, 87].max(); g_dn = dn[:, 78].min() - dn[:, 77].max()
    mid = vbm + 0.5 * g_up; A, S, C = struct(f'{d}/scf_FM.in'); r = [e * 100, g_up, g_dn, GAP[round(e * 100)] - g_up]
    for s, w in (('Ni1', 3), ('Mn1', 4)):
        D = np.loadtxt(glob.glob(f'{d}/*({s})_wfc#{w}(d)')[0]); E = D[:, 0]; m = D[:, 3:13].reshape(len(E), 5, 2)
        eg = m[:, 0] + m[:, 4]; t2 = m[:, 1] + m[:, 2] + m[:, 3]; P[(e, s)] = (E - vbm, eg, t2); k = E <= mid
        r += [TZ(eg[k, 0], E[k]), TZ(eg[k, 1], E[k]), TZ(t2[k, 0], E[k]), TZ(t2[k, 1], E[k]), leak(A, S, C, s)]
    rows.append(r)
H = 'eps(%) Eg_up Eg_dn dEg(gap.dat-bu) Ni_eg_up Ni_eg_dn Ni_t2g_up Ni_t2g_dn Ni_leak Mn_eg_up Mn_eg_dn Mn_t2g_up Mn_t2g_dn Mn_leak'
np.savetxt(f'{OUT}/eg_t2g_occ.dat', np.array(rows), fmt='%8.4f', header=H); print(H); [print(' '.join('%7.3f' % x for x in r)) for r in rows]
fig, ax = plt.subplots(2, 3, figsize=(10, 6), sharex=True, sharey='row')
for j, e in enumerate([-0.05, 0.0, 0.05]):
    for i, s in enumerate(('Ni1', 'Mn1')):
        E, eg, t2 = P[(e, s)]; a = ax[i, j]
        for y, c, lab in ((eg, '#0072B2', '$e_g$'), (t2, '#D55E00', '$t_{2g}$')):
            a.plot(E, y[:, 0], c=c, lw=1, label=lab); a.plot(E, -y[:, 1], c=c, lw=1)
        a.axhline(0, c='k', lw=.5); a.axvline(0, c='gray', ls=':', lw=.8); a.set_xlim(-8, 5)
        if i == 0: a.set_title(r'$\varepsilon$ = %+d %%' % round(e * 100))
        if j == 0: a.set_ylabel(s[:2] + ' 3$d$ PDOS (states/eV)')
        if i == 1: a.set_xlabel('$E - E_{VBM}$ (eV)')
ax[0, 0].legend(frameon=False); fig.tight_layout()
for x in ('png', 'pdf'): fig.savefig(f'{OUT}/fig_egt2g.{x}', dpi=300)
print('saved:', OUT)
