# make_figs.py - Paper-1 R2 figures: Fig.2 (merged FM/AFM + inset AFM-FM, FiM), Fig.7 (Eg vs |Q3|), Fig.8 (parity, single-axis panels), GA
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
MM = 1/25.4; B, V, G, K = '#0072B2', '#D55E00', '#009E73', '#4D4D4D'
mpl.rcParams.update({"pdf.fonttype":42,"svg.fonttype":"path","font.family":"sans-serif","font.sans-serif":["DejaVu Sans"],
  "font.size":8,"axes.labelsize":8.5,"legend.fontsize":7.5,"xtick.labelsize":7,"ytick.labelsize":7,"axes.linewidth":0.8,
  "xtick.direction":"in","ytick.direction":"in","xtick.top":True,"ytick.right":True,"legend.frameon":False,"figure.dpi":200})
import os
D = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(D, '..', '..', 'data') + '/'          # repository data/ folder
os.makedirs('figs', exist_ok=True)
E = np.loadtxt(R+'energies.dat'); e = E[:,0]*100
q = np.loadtxt(R+'q3_geo.dat', usecols=(0,3,7)); Q = np.abs(q[:,1])
# R1 band gaps: min(6x6x4 mesh, high-symmetry path); -1..-4 mesh only (band_edges_R1.dat, gap.dat)
UP, DN = np.loadtxt(os.path.join(D, '..', 'data', 'gap_R2.dat'), usecols=(1, 2), unpack=True)   # Table 3
def save(fig, n):
    for ext, kw in (("pdf",{}),("png",{"dpi":600})): fig.savefig(f"figs/{n}.{ext}", **kw)
# ---- Fig. 2 ----
fig, a = plt.subplots(figsize=(90*MM, 70*MM))
a.plot(e, E[:,6], 'o-', c=B, ms=4, lw=1.2, label='FM'); a.plot(e, E[:,7], 's--', c=V, ms=4, lw=1.2, mfc='white', label='C-type AFM')
a.set_xlabel(r'$\varepsilon$ (%)'); a.set_ylabel(r'$\Delta E$ relative to FM at $\varepsilon=0$ (meV/f.u.)')
a.set_xlim(-5.6,5.6); a.set_ylim(-30, 930); a.xaxis.set_major_locator(MultipleLocator(1)); a.legend(loc='upper left')
ins = a.inset_axes([0.53, 0.50, 0.33, 0.45])
ins.plot(e, E[:,4], 's-', c=V, ms=3, lw=1, mfc='white', label='AFM')
m = ~np.isnan(E[:,5]); ins.plot(e[m], E[m,5], 'D', c=G, ms=4, label='FiM')
ins.set_ylim(110, 205); ins.set_xlim(-5.6,5.6); ins.tick_params(labelsize=6); ins.set_xlabel(r'$\varepsilon$ (%)', fontsize=6.5, labelpad=1)
ins.set_ylabel(r'$E-E_{\rm FM}$ (meV/f.u.)', fontsize=6.5, labelpad=2); ins.legend(fontsize=6, loc='center', handletextpad=0.3)
fig.tight_layout(pad=0.4); save(fig, 'Fig2_energies')
# ---- Fig. 7 ----
p = np.polyfit(Q, UP, 1); r = np.corrcoef(Q, UP)[0,1]; rm = np.sqrt(np.mean((UP-np.polyval(p,Q))**2))*1000
fig, a = plt.subplots(figsize=(90*MM, 68*MM)); xs = np.linspace(0, 0.19, 50)
a.plot(xs, np.polyval(p, xs), '--', c=K, lw=1)
c = e < 0; t = e > 0; z = e == 0
a.plot(Q[c], UP[c], 'o', c=B, ms=5, label=r'$\varepsilon<0$'); a.plot(Q[t], UP[t], 'o', c=V, ms=5, label=r'$\varepsilon>0$')
a.plot(Q[z], UP[z], 'D', c=K, ms=5, mfc='white', label=r'$\varepsilon=0$')
a.set_xlabel(r'$|Q_3|$ of NiO$_6$ (Å)'); a.set_ylabel(r'$E_g^{\uparrow}$ (eV)'); a.legend(loc='lower left')
a.text(0.97, 0.95, r'$E_g = %.3f - %.3f\,|Q_3|$' '\n' r'$r = %.3f$, RMSE = %.0f meV' % (p[1], -p[0], r, rm), transform=a.transAxes, ha='right', va='top', fontsize=7)
fig.tight_layout(pad=0.4); save(fig, 'Fig7_gap_Q3')
print('Fig7 fit', p, r, rm)
# ---- Fig. 8 ----
i0 = 5; k = np.arange(1,6)
T = lambda v: np.array([(v[i0+j]-v[i0-j])/2 for j in k])
def peak(v):
    cc = np.polyfit(e, v, 4); xx = np.linspace(-5,5,20001); return xx[np.argmax(np.polyval(cc,xx))], cc
fig, ax = plt.subplots(1, 3, figsize=(180*MM, 62*MM)); xx = np.linspace(-5,5,400)
for a, v, col, mk, lab in ((ax[0], UP, B, 'o', 'majority'), (ax[1], DN, V, 's', 'minority')):
    pk, cc = peak(v); a.plot(xx, np.polyval(cc, xx), '-', c=col, lw=1); a.plot(e, v, mk, c=col, ms=4, mfc='white' if mk=='s' else col)
    a.axvline(pk, c=col, ls=':', lw=0.8); a.axvspan(-5.6, -3.45, color='0.92', lw=0, zorder=0)
    a.text(pk+0.25, v.min()+0.12*(v.max()-v.min()), 'max at\n%+.2f %%' % pk, color=col, fontsize=7)
    a.set_xlabel(r'$\varepsilon$ (%)'); a.set_ylabel(r'$E_g$ %s (eV)' % lab); a.set_xlim(-5.6,5.6); a.xaxis.set_major_locator(MultipleLocator(2))
    print(lab, 'peak', pk)
b = ax[2]; b.axhline(0, c=K, lw=0.6)
b.plot(np.r_[0,k], np.r_[0,T(UP)*1000], 'o-', c=B, ms=4, lw=1.2, label='majority'); b.plot(np.r_[0,k], np.r_[0,T(DN)*1000], 's--', c=V, ms=4, lw=1.2, mfc='white', label='minority')
b.set_xlabel(r'$|\varepsilon|$ (%)'); b.set_ylabel(r'odd component $T$ (meV)'); b.set_xlim(-0.2,5.3); b.legend(loc='lower left')
for a, s in zip(ax, 'abc'): a.text(0.03, 0.93, f'({s})', transform=a.transAxes, fontweight='bold')
fig.tight_layout(pad=0.5); save(fig, 'Fig8_parity')
print('T up', np.round(T(UP)*1000,1), 'T dn', np.round(T(DN)*1000,1))
# ---- Graphical abstract (single-axis panels) ----
fig, ax = plt.subplots(1, 2, figsize=(13.28/2.54*1.2, 5.31/2.54*1.2))
for a, v, col, mk, ttl in ((ax[0], UP, B, 'o', 'majority gap: even in strain'), (ax[1], DN, V, 's', 'minority gap: odd part 136 meV')):
    pk, cc = peak(v); a.plot(xx, np.polyval(cc, xx), '-', c=col, lw=1.2); a.plot(e, v, mk, c=col, ms=5, mfc='white' if mk=='s' else col)
    a.axvline(0, c=K, lw=0.6); a.axvline(pk, c=col, ls=':', lw=1); a.set_title(ttl, fontsize=8.5)
    a.set_xlabel(r'biaxial strain $\varepsilon$ (%)'); a.set_ylabel(r'$E_g$ (eV)')
fig.suptitle(r'La$_2$NiMnO$_6$, PBE+U: channel-selective parity of the strain response', fontsize=9)
fig.tight_layout(pad=0.5); fig.savefig('figs/Graphical_abstract.png', dpi=300); fig.savefig('figs/Graphical_abstract.pdf')
