# ph_gen.py - CMS R1 madde 1: Gamma fononu (phonopy, 1x1x1, 0.01 A). Kullanim: python3 ph_gen.py scf_FM.in WORKDIR
import sys, os, re, numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
src, W = sys.argv[1], sys.argv[2]; L = open(src).read().splitlines()
i = [k for k, l in enumerate(L) if l.strip().upper().startswith('CELL_PARAMETERS')][0]
A = np.array([[float(x) for x in L[i + j].split()[:3]] for j in (1, 2, 3)])
j = [k for k, l in enumerate(L) if l.strip().upper().startswith('ATOMIC_POSITIONS')][0]
S = [re.sub(r'\d', '', l.split()[0]) for l in L[j + 1:j + 21]]
X = [[float(v) for v in l.split()[1:4]] for l in L[j + 1:j + 21]]
ph = Phonopy(PhonopyAtoms(symbols=S, cell=A, scaled_positions=X), supercell_matrix=np.eye(3, dtype=int), symprec=1e-3)
ph.generate_displacements(distance=0.01)
print('uzay grubu', ph.symmetry.get_international_table(), '| yer degistirme', len(ph.supercells_with_displacements))
ph.save(os.path.join(W, 'phonopy_disp.yaml'))
H = """&CONTROL
 calculation='scf', prefix='ph', outdir='./tmp', pseudo_dir='/arf/scratch/nparan/lnmo/pseudo'
 tprnfor=.true., verbosity='low', disk_io='low'
/
&SYSTEM
 ibrav=0, nat=20, ntyp=4, ecutwfc=90.0, ecutrho=900.0, nspin=2
 starting_magnetization(2)=0.5, starting_magnetization(3)=0.5
 occupations='smearing', smearing='mv', degauss=0.01
/
&ELECTRONS
 diagonalization='david', mixing_mode='plain', mixing_beta=0.7, conv_thr=1.0d-10
/
ATOMIC_SPECIES
 La 138.905 La.pbe-spfn-kjpaw_psl.1.0.0.UPF
 Mn  54.938 Mn.pbe-spn-kjpaw_psl.0.3.1.UPF
 Ni  58.693 Ni.pbe-n-kjpaw_psl.0.1.UPF
 O   15.999 O.pbe-n-kjpaw_psl.1.0.0.UPF
K_POINTS automatic
 6 6 4 0 0 0
HUBBARD ortho-atomic
 U Mn-3d 4.0
 U Ni-3d 5.1
"""
cells = [ph.supercell] + list(ph.supercells_with_displacements)
for n, c in enumerate(cells):
    os.makedirs(f'{W}/d{n}', exist_ok=True)
    with open(f'{W}/d{n}/pw.in', 'w') as f:
        f.write((H.replace('ntyp=4,', 'ntyp=4, nosym=.true.,') if n == 0 else H) + 'CELL_PARAMETERS angstrom\n' + ''.join('%16.10f%16.10f%16.10f\n' % tuple(r) for r in c.cell))
        f.write('ATOMIC_POSITIONS crystal\n' + ''.join('%-3s%16.10f%16.10f%16.10f\n' % (s, *p) for s, p in zip(c.symbols, c.scaled_positions)))
open(f'{W}/ndisp', 'w').write(str(len(cells) - 1))
