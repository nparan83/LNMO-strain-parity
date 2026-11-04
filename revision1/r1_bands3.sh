#!/bin/bash
# r1_bands3.sh - ek bant yollari (varsayilan +1..+4 %): save kopyasi + nscf_FM.in sablonu + cg ayarlari + gonderim
P=${P:-/arf/scratch/nparan/lnmo/LNMO_Strain_Production}
L=${*:-"0.01 0.02 0.03 0.04"}; for e in $L; do
 d=$P/strain_$e/FM; B=$d/R1_bands; pf=lnmo_FM_$e; S=$d/tmp/$pf.save
 [ -f $S/data-file-schema.xml ] || { echo "$e: save YOK, atlandi"; continue; }
 mkdir -p $B/tmp/$pf.save; cd $B
 for f in $S/*; do case $(basename $f) in wfc*) ;; *) cp -r $f tmp/$pf.save/ ;; esac; done
 sed -e "s/'nscf'/'bands'/" -e '/K_POINTS/,+1d' $d/nscf_FM.in > bands.in
 grep -q diago_thr_init bands.in || sed -i -e "/diagonalization/a\    diago_cg_maxiter = 100" -e "/diagonalization/a\    diago_thr_init  = 1.0d-6" bands.in
 printf 'K_POINTS crystal_b\n10\n' >> bands.in
 for k in "0 0 0" ".5 0 0" ".5 .5 0" "0 .5 0" "0 0 0" "0 0 .5" ".5 0 .5" ".5 .5 .5" "0 .5 .5"; do echo "$k 20" >> bands.in; done
 echo "0 0 .5 1" >> bands.in
 for s in 1 2; do printf "&BANDS\n prefix='%s', outdir='./tmp', filband='bands_s%d.dat', spin_component=%d, lsym=.false.\n/\n" $pf $s $s > bandsx_s$s.in; done
 sed -e "s/bnd-0.05/bnd$e/" $P/strain_-0.05/FM/R1_bands/run.sh > run.sh
 echo "== $e: $(grep -ciE 'diago' bands.in) diago satiri, $(grep -c ' 20$' bands.in) segment"
 [ -z "$DRY" ] && sbatch run.sh
done
