# bl.py - B-O bag uzunluklari (Mn1, Ni1; en kisa 6, A) uretim scf_FM.in'den. Kullanim: python3 bl.py [eps ...]
import sys, numpy as np
P = '/arf/scratch/nparan/lnmo/LNMO_Strain_Production'
for e in (sys.argv[1:] or ['-0.05', '0.00', '0.05']):
    f = e if e.endswith('.in') else f'{P}/strain_{e}/FM/scf_FM.in'
    L = open(f).read().splitlines(); i = [k for k, l in enumerate(L) if l.strip().upper().startswith('CELL_PARAMETERS')][0]
    A = np.array([[float(x) for x in L[i + j].split()[:3]] for j in (1, 2, 3)]); j = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')][0]
    S = [l.split()[0] for l in L[j + 1:j + 21]]; C = np.array([[float(v) for v in l.split()[1:4]] for l in L[j + 1:j + 21]]) @ A
    for s in ('Mn1', 'Ni1'):
        c = C[S.index(s)]; v = np.array([C[k] + (np.array(n) - 1) @ A - c for k, t in enumerate(S) if t.startswith('O') for n in np.ndindex(3, 3, 3)])
        v = v[np.argsort(np.linalg.norm(v, axis=1))][:6]; d = np.linalg.norm(v, axis=1); z = np.abs(v[:, 2]) / d > 0.7
        print('%6s %s  duzlem-ici: %s | apikal: %s | ort %.4f' % (e, s, ' '.join('%.4f' % x for x in sorted(d[~z])), ' '.join('%.4f' % x for x in sorted(d[z])), d.mean()))
