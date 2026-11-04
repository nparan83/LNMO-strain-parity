#!/bin/bash
# r1_bands2.sh - bands girdisini URETIM nscf_FM.in sablonundan yeniden kur (c_bands yakinsama ayarlari onda) ve yeniden gonder
P=${P:-/arf/scratch/nparan/lnmo/LNMO_Strain_Production}
for e in -0.05 0.00 0.05; do
 d=$P/strain_$e/FM; B=$d/R1_bands; cd $B || exit 1
 rm -f CRASH bands.out bandsx_s*.out slurm_*.out
 sed -e "s/'nscf'/'bands'/" -e '/K_POINTS/,+1d' $d/nscf_FM.in > bands.in
 printf 'K_POINTS crystal_b\n10\n' >> bands.in
 for k in "0 0 0" ".5 0 0" ".5 .5 0" "0 .5 0" "0 0 0" "0 0 .5" ".5 0 .5" ".5 .5 .5" "0 .5 .5"; do echo "$k 20" >> bands.in; done
 echo "0 0 .5 1" >> bands.in
 echo "== $e"; grep -iE "calculation|nbnd|diago|conv_thr|nosym|occupations|degauss|outdir|prefix" bands.in
 [ -z "$DRY" ] && sbatch run.sh
done
