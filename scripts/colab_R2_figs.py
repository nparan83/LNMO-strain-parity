# ===== Paper-1 R2: Figs. 2, 3, 7, 8 + graphical abstract (one cell). Output: Research/DFT/LNMO_2/R1/Figures =====
import os, io, urllib.request, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from google.colab import drive; drive.mount('/content/drive')

FIGDIR = "/content/drive/MyDrive/Research/DFT/LNMO_2/R1/Figures"; os.makedirs(FIGDIR, exist_ok=True)
GH = "https://raw.githubusercontent.com/nparan83/LNMO-strain-parity/main/data/"
def gh(name, **kw):                       # data file from the public repository (same files as the paper)
    return np.loadtxt(io.StringIO(urllib.request.urlopen(GH + name).read().decode()), **kw)
MM = 1/25.4; CB, CV, CG, CK = '#0072B2', '#D55E00', '#009E73', '#4D4D4D'   # Okabe-Ito (unique names: other cells reuse B/G/K)
def pub_style():
    plt.rcdefaults()
    mpl.rcParams.update({"pdf.fonttype":42, "svg.fonttype":"path", "font.family":"sans-serif", "font.sans-serif":["DejaVu Sans"],
      "font.size":8, "axes.labelsize":8.5, "legend.fontsize":7.5, "xtick.labelsize":7, "ytick.labelsize":7, "axes.linewidth":0.8,
      "xtick.direction":"in", "ytick.direction":"in", "xtick.top":True, "ytick.right":True, "legend.frameon":False, "figure.dpi":150})
def save(fig, name):
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600})): fig.savefig(f"{FIGDIR}/{name}.{ext}", **kw)
    print("saved:", f"{FIGDIR}/{name}.pdf/.png")
# R1 band gaps of the a-a-c+ structure (eV) = min(6x6x4 mesh, high-symmetry path); -1..-4 % mesh only.
# Sources: band_edges_R1.dat (path, Drive LNMO_2/R1) and gap.dat (mesh). Same values as Table 3 of the R2 manuscript.
EPS_R1 = np.arange(-5, 6)
EG_UP = np.array([1.0755,1.1430,1.1942,1.2406,1.2764,1.2947,1.2930,1.2526,1.1911,1.1325,1.0828])
EG_DN = np.array([3.1580,3.2022,3.2304,3.2437,3.2439,3.2321,3.1755,3.1121,3.0398,2.9656,2.8861])
# ===== Fig. 2 — FM/AFM energies; inset: AFM-FM at all strains and FiM-FM at -5, 0, +5 % =====
pub_style()
E = gh("energies.dat"); e = E[:,0]*100     # cols: eps E_FM E_AFM E_FiM dE_AFM-FM dE_FiM-FM E_FM-E_FM(0) E_AFM-E_FM(0) Mtot
fig, a = plt.subplots(figsize=(90*MM, 70*MM))
a.plot(e, E[:,6], 'o-', c=CB, ms=4, lw=1.2, label='FM'); a.plot(e, E[:,7], 's--', c=CV, ms=4, lw=1.2, mfc='white', label='C-type AFM')
a.set_xlabel(r'$\varepsilon$ (%)'); a.set_ylabel(r'$\Delta E$ relative to FM at $\varepsilon=0$ (meV/f.u.)')
a.set_xlim(-5.6, 5.6); a.set_ylim(-30, 930); a.xaxis.set_major_locator(MultipleLocator(1)); a.legend(loc='upper left')
ins = a.inset_axes([0.53, 0.50, 0.33, 0.45])
ins.plot(e, E[:,4], 's-', c=CV, ms=3, lw=1, mfc='white', label='AFM')
m = ~np.isnan(E[:,5]); ins.plot(e[m], E[m,5], 'D', c=CG, ms=4, label='FiM')
ins.set_ylim(110, 205); ins.set_xlim(-5.6, 5.6); ins.tick_params(labelsize=6)
ins.set_xlabel(r'$\varepsilon$ (%)', fontsize=6.5, labelpad=1); ins.set_ylabel(r'$E-E_{\rm FM}$ (meV/f.u.)', fontsize=6.5, labelpad=2)
ins.legend(fontsize=6, loc='center', handletextpad=0.3)
fig.tight_layout(pad=0.4); save(fig, 'Fig2_energies'); plt.show()
print('AFM-FM min/max %.1f / %.1f meV/f.u.; FiM-FM %s' % (np.nanmin(E[:,4]), np.nanmax(E[:,4]), np.round(E[m,5], 1)))
# ===== Fig. 3 — out-of-plane lattice parameter c_perp = V/|a x b| = c sin(beta) =====
# NOTE: the old cell plotted |c|/|a| of the oblique P2_1/c cell (1.85 -> 1.64). That is not c_perp and does not match the caption.
pub_style()
C = gh("cperp.dat"); e = C[:,0]*100; cp = C[:,5]          # cols: eps |a| |b| |c| V c_perp
i0 = int(np.argmin(abs(e)))
fig, a = plt.subplots(figsize=(90*MM, 65*MM))
a.plot(e, cp, 'o-', c=CB, ms=4, lw=1.2)
a.set_xlabel(r'$\varepsilon$ (%)'); a.set_ylabel(r'$c_\perp = c\,\sin\beta$ (Å)'); a.set_xlim(-5.6, 5.6); a.set_ylim(7.50, 8.27); a.xaxis.set_major_locator(MultipleLocator(1))
for j, ha in ((0, 'left'), (-1, 'right')):
    a.annotate('%+.1f %%' % (100*(cp[j]/cp[i0]-1)), (e[j], cp[j]), (8, 0) if j == 0 else (-4, -12), textcoords='offset points', ha=ha, va='center', fontsize=7)
fig.tight_layout(pad=0.4); save(fig, 'Fig3_cperp'); plt.show()
print('c_perp(-5, 0, +5) = %.3f, %.3f, %.3f A; d ln c_perp/d eps (+-1 %%) = %.3f' % (cp[0], cp[i0], cp[-1], (np.log(cp[i0+1]) - np.log(cp[i0-1]))/0.02))
# ===== Fig. 7 — majority gap vs |Q3| of NiO6 (old Fig. 6, R1 gap values) =====
pub_style()
q = gh("q3_geo.dat", usecols=(0, 3, 7)); Q = np.abs(q[:,1]); e = q[:,0]     # Q3 = (2/sqrt6)(d_ap - d_eq)
p = np.polyfit(Q, EG_UP, 1); r = np.corrcoef(Q, EG_UP)[0,1]; rm = 1000*np.sqrt(np.mean((EG_UP - np.polyval(p, Q))**2))
fig, a = plt.subplots(figsize=(90*MM, 68*MM)); xs = np.linspace(0, 0.19, 50)
a.plot(xs, np.polyval(p, xs), '--', c=CK, lw=1)
for msk, col, mk, lab in ((e < 0, CB, 'o', r'$\varepsilon<0$'), (e > 0, CV, 'o', r'$\varepsilon>0$')): a.plot(Q[msk], EG_UP[msk], mk, c=col, ms=5, label=lab)
a.plot(Q[e == 0], EG_UP[e == 0], 'D', c=CK, ms=5, mfc='white', label=r'$\varepsilon=0$')
a.set_xlabel(r'$|Q_3|$ of NiO$_6$ (Å)'); a.set_ylabel(r'$E_g^{\uparrow}$ (eV)'); a.legend(loc='lower left')
a.text(0.97, 0.95, r'$E_g = %.3f - %.3f\,|Q_3|$' '\n' r'$r = %.3f$, RMSE = %.0f meV' % (p[1], -p[0], r, rm), transform=a.transAxes, ha='right', va='top', fontsize=7)
fig.tight_layout(pad=0.4); save(fig, 'Fig7_gap_Q3'); plt.show()
print('fit: Eg = %.3f - %.3f|Q3|, r = %.3f, RMSE = %.1f meV  (manuscript: 1.327 - 1.370|Q3|, r = -0.983, 14 meV)' % (p[1], -p[0], r, rm))
# ===== Fig. 8 — parity of the spin-resolved gaps (old Fig. 7; single-axis panels, R1 gap values) =====
pub_style()
e = EPS_R1.astype(float); i0 = 5; k = np.arange(1, 6)
T = lambda v: np.array([(v[i0+j] - v[i0-j])/2 for j in k])
def peak(v):
    cc = np.polyfit(e, v, 4); xx = np.linspace(-5, 5, 20001); return xx[np.argmax(np.polyval(cc, xx))], cc
fig, ax = plt.subplots(1, 3, figsize=(180*MM, 62*MM)); xx = np.linspace(-5, 5, 400)
for a, v, col, mk, lab in ((ax[0], EG_UP, CB, 'o', 'majority'), (ax[1], EG_DN, CV, 's', 'minority')):
    pk, cc = peak(v); a.plot(xx, np.polyval(cc, xx), '-', c=col, lw=1); a.plot(e, v, mk, c=col, ms=4, mfc='white' if mk == 's' else col)
    a.axvline(pk, c=col, ls=':', lw=0.8); a.axvspan(-5.6, -3.45, color='0.92', lw=0, zorder=0)   # a-a-c+ metastable below eps* (Sec. 3.7)
    a.text(pk + 0.25, v.min() + 0.12*(v.max() - v.min()), 'max at\n%+.2f %%' % pk, color=col, fontsize=7)
    a.set_xlabel(r'$\varepsilon$ (%)'); a.set_ylabel(r'$E_g$ %s (eV)' % lab); a.set_xlim(-5.6, 5.6); a.xaxis.set_major_locator(MultipleLocator(2))
    print(lab, 'max at %+.2f %%' % pk)
b = ax[2]; b.axhline(0, c=CK, lw=0.6)
b.plot(np.r_[0, k], np.r_[0, T(EG_UP)*1000], 'o-', c=CB, ms=4, lw=1.2, label='majority')
b.plot(np.r_[0, k], np.r_[0, T(EG_DN)*1000], 's--', c=CV, ms=4, lw=1.2, mfc='white', label='minority')
b.set_xlabel(r'$|\varepsilon|$ (%)'); b.set_ylabel(r'odd component $T$ (meV)'); b.set_xlim(-0.2, 5.3); b.legend(loc='lower left')
for a, s in zip(ax, 'abc'): a.text(0.03, 0.93, f'({s})', transform=a.transAxes, fontweight='bold')
fig.tight_layout(pad=0.5); save(fig, 'Fig8_parity'); plt.show()
print('T_up (meV):', np.round(T(EG_UP)*1000, 1), ' T_dn (meV):', np.round(T(EG_DN)*1000, 1))
print('manuscript Table 4: T_up 8.3 6.0 -1.5 -5.2 3.7 | T_dn -34.2 -65.8 -95.3 -118.3 -136.0 | maxima +0.03 %, -1.60 %')
# ===== Graphical abstract =====
pub_style()
e = EPS_R1.astype(float)
def peak(v):
    cc = np.polyfit(e, v, 4); xx = np.linspace(-5, 5, 20001); return xx[np.argmax(np.polyval(cc, xx))], cc
fig, ax = plt.subplots(1, 2, figsize=(13.28/2.54*1.2, 5.31/2.54*1.2)); xx = np.linspace(-5, 5, 400)
for a, v, col, mk, ttl in ((ax[0], EG_UP, CB, 'o', 'majority gap: even in strain'), (ax[1], EG_DN, CV, 's', 'minority gap: odd part 136 meV')):
    pk, cc = peak(v); a.plot(xx, np.polyval(cc, xx), '-', c=col, lw=1.2); a.plot(e, v, mk, c=col, ms=5, mfc='white' if mk == 's' else col)
    a.axvline(0, c=CK, lw=0.6); a.axvline(pk, c=col, ls=':', lw=1); a.set_title(ttl, fontsize=8.5)
    a.set_xlabel(r'biaxial strain $\varepsilon$ (%)'); a.set_ylabel(r'$E_g$ (eV)')
fig.suptitle(r'La$_2$NiMnO$_6$, PBE+U: channel-selective parity of the strain response', fontsize=9)
fig.tight_layout(pad=0.5)
for ext in ('png', 'pdf'): fig.savefig(f'{FIGDIR}/Graphical_abstract.{ext}', dpi=300)
plt.show(); print('saved:', f'{FIGDIR}/Graphical_abstract.png/.pdf')
