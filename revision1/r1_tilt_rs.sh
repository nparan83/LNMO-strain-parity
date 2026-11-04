#!/bin/bash
# r1_tilt_rs.sh - nstep=50 sinirina takilan a-a-a- vc-relax'lari SON geometriden surdur (nstep=200) + scf; durum raporu
R=/arf/scratch/nparan/R1; P=/arf/scratch/nparan/lnmo/LNMO_Strain_Production
L=${*:-"0.00 -0.02"}; for e in $L; do
 W=$R/tilt/eps_$e; cd $W || exit 1; T=$P/strain_$e/FM/scf_FM.in
 python3 $R/mkin.py restart vcrelax.in vcrelax.out > vcrelax2.in || exit 1
 echo "== $e entalpi (son 6): $(grep 'enthalpy new' vcrelax.out | tail -6 | awk '{printf "%.6f ", $4}')"
 python3 - vcrelax2.in << 'PY'
import sys; exec(open('/arf/scratch/nparan/R1/tilt_post.py').read().split('E = lambda')[0])
A, S, X = struct(sys.argv[1]); w = rotz(A, S, X); print('   son geometri wz Mn1/Ni1: %+.1f / %+.1f  (zit isaret = c-, ayni = c+)' % (w['Mn1'], w['Ni1']))
PY
 sed -e 's/vcrelax\.in/vcrelax2.in/; s/vcrelax\.out/vcrelax2.out/g; s/-J aaa/-J rs/' run.sh > run2.sh
 [ -z "$DRY" ] && sbatch run2.sh
done
