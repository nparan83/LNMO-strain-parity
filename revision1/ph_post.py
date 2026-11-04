# ph_post.py - Gamma frekanslari. Kullanim: python3 ph_post.py WORKDIR [uretim scf_FM.out]
import sys, re, numpy as np, phonopy
W = sys.argv[1]; RYB = 25.71104309541616  # Ry/bohr -> eV/A
def read(f):
    t = open(f).read(); E = re.findall(r'!\s+total energy\s+=\s+(\S+)', t)
    k = t.rfind('Forces acting on atoms'); ok = 'JOB DONE' in t
    F = [list(map(float, l.split('=')[1].split()[:3])) for l in t[k:].splitlines() if 'force =' in l and 'atom' in l][:20]
    return (float(E[-1]) if E else None), np.array(F) * RYB, ok
n = int(open(f'{W}/ndisp').read()); R = [read(f'{W}/d{i}/pw.out') for i in range(n + 1)]
bad = [i for i, r in enumerate(R) if not r[2] or r[1].shape != (20, 3)]
if bad: sys.exit(f'bitmemis/bozuk: {bad}')
print('E(mukemmel) = %.8f Ry' % R[0][0], ' |F|max mukemmel = %.2e eV/A' % np.abs(R[0][1]).max())
if len(sys.argv) > 2:
    e = re.findall(r'!\s+total energy\s+=\s+(\S+)', open(sys.argv[2]).read())
    print('uretim E = %s Ry, fark = %.2e Ry' % (e[-1], R[0][0] - float(e[-1])))
ph = phonopy.load(f'{W}/phonopy_disp.yaml', produce_fc=False, log_level=0)
ph.forces = np.array([r[1] - R[0][1] for r in R[1:]])
ph.produce_force_constants(); ph.symmetrize_force_constants()
ph.run_qpoints([[0, 0, 0]], with_eigenvectors=False)
f = ph.get_qpoints_dict()['frequencies'][0]; cm = f * 33.35641
np.savetxt(f'{W}/gamma_freq.dat', np.c_[np.arange(1, len(f) + 1), f, cm], fmt='%4d %10.4f %10.2f', header='mod  THz  cm-1')
print('en dusuk 6 (cm-1):', np.round(cm[:6], 2)); print('sanal mod (< -5 cm-1):', int((cm < -5).sum()))
