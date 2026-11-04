#!/bin/bash
# r1_aaa_gap.sh - a-a-a- yapisinda bant araligi (eps = -4, -5 %): scf -> nscf (6x6x4, uretim NSCF ayarlari) -> bands (yol) -> bands.x
R=/arf/scratch/nparan/R1; L=${*:-"-0.04 -0.05"}
for e in $L; do
 G=$R/tilt/eps_$e/gap; mkdir -p $G; cd $G || exit 1; cp ../scf.in scf.in
 sed -e "s/'scf'/'nscf'/" -e "/verbosity/d;/diagonalization/d;/diago_/d;/nbnd/d" -e "s/degauss *=.*/degauss       = 0.002/" \
 -e "/calculation/a\    verbosity     = 'high'" -e "/nspin/a\    nbnd          = 100" -e "/^ *\&ELECTRONS/Ia\    diagonalization = 'cg'\n    diago_thr_init  = 1.0d-6\n    diago_cg_maxiter = 100" scf.in > nscf.in
 sed -e "s/'nscf'/'bands'/" -e '/K_POINTS/,+1d' nscf.in > bands.in
 printf 'K_POINTS crystal_b\n10\n' >> bands.in
 for k in "0 0 0" ".5 0 0" ".5 .5 0" "0 .5 0" "0 0 0" "0 0 .5" ".5 0 .5" ".5 .5 .5" "0 .5 .5"; do echo "$k 20" >> bands.in; done; echo "0 0 .5 1" >> bands.in
 pf=$(grep -i prefix scf.in | sed "s/.*=//;s/[\"', ]//g")
 for s in 1 2; do printf "&BANDS\n prefix='%s', outdir='./tmp', filband='bands_s%d.dat', spin_component=%d, lsym=.false.\n/\n" $pf $s $s > bandsx_s$s.in; done
 cat > run.sh << S
#!/bin/bash
#SBATCH -J gap$e
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=14:00:00
#SBATCH --output=slurm_%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023 lib/hdf5/1.14.3-oneapi-2023.0 apps/espresso/7.2-oneapi-2023.0
mpirun pw.x -nk 4 -in scf.in > scf.out && mpirun pw.x -nk 4 -in nscf.in > nscf.out && cp nscf.out nscf_keep.out && \
mpirun pw.x -nk 4 -in bands.in > bands.out && for s in 1 2; do mpirun -np 1 bands.x -in bandsx_s\$s.in > bandsx_s\$s.out; done
rm -rf tmp
S
 echo "== $e: $(grep -cE "'nscf'|nbnd|'cg'|diago_thr|degauss *= *0.002|verbosity" nscf.in)/6 nscf ayari, $(grep -c ' 20$' bands.in) yol satiri (+nat)"
 [ -z "$DRY" ] && sbatch run.sh
done
