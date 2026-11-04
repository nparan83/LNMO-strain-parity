# el_gen.py - CMS R1 madde 3: C_ij (eps=0), +-0.5 % Voigt zorlanmalari, iyonlar gevsek (relax). Kullanim: python3 el_gen.py scf_FM.in WORKDIR
import sys, os, numpy as np
src, W = sys.argv[1], sys.argv[2]; d = 0.005; L = open(src).read().splitlines()
i = [k for k, l in enumerate(L) if l.strip().upper().startswith('CELL_PARAMETERS')][0]
A = np.array([[float(x) for x in L[i + j].split()[:3]] for j in (1, 2, 3)])
V = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]                     # Voigt: xx yy zz yz xz xy
n = 0
for j, (a, b) in enumerate(V):
    for s, tag in ((1, 'p'), (-1, 'm')):
        E = np.zeros((3, 3)); E[a, b] += s * d / (1 if a == b else 2); E[b, a] = E[a, b]  # muhendislik kaymasi -> tensor d/2
        A2 = A @ (np.eye(3) + E); n += 1; D = f'{W}/{n}'; os.makedirs(D, exist_ok=True)
        out = []
        for k, l in enumerate(L):
            if i < k <= i + 3: out.append('%16.10f%16.10f%16.10f' % tuple(A2[k - i - 1])); continue
            if "'scf'" in l: l = l.replace("'scf'", "'relax'") + "\n    tstress       = .true.\n    tprnfor       = .true.\n    nstep         = 100"
            out.append(l)
        t = '\n'.join(out) + '\n'
        t = t.replace('ATOMIC_SPECIES', "&IONS\n    ion_dynamics    = 'bfgs'\n/\nATOMIC_SPECIES", 1)
        open(f'{D}/pw.in', 'w').write(t); open(f'{D}/tag', 'w').write(f'{j} {s}\n')
print('zorlanma dosyasi:', n, '| delta =', d)
