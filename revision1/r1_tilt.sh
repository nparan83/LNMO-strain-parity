#!/bin/bash
# r1_tilt.sh - CMS R1 madde 4: a-a-a- (R-3 kokenli) rakip tilt deseni; uretimle ayni hucre/kisit (epitaxial_ab), vc-relax + scf
P=/arf/scratch/nparan/lnmo/LNMO_Strain_Production; R=/arf/scratch/nparan/R1
L=${*:-"0.00 -0.02 -0.04"}; for e in $L; do
 W=$R/tilt/eps_$e; mkdir -p $W; T=$P/strain_$e/FM/scf_FM.in
 python3 $R/tilt_gen.py $T aaa 6 6 > $W/pos.txt || exit 1
 python3 $R/mkin.py relax $T $W/pos.txt > $W/vcrelax.in || exit 1
 cat > $W/run.sh << S
#!/bin/bash
#SBATCH -J aaa$e
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=3-00:00:00
#SBATCH --output=slurm_%j.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023 lib/hdf5/1.14.3-oneapi-2023.0 apps/espresso/7.2-oneapi-2023.0
mpirun pw.x -nk 4 -in vcrelax.in > vcrelax.out
python3 $R/mkin.py scf $T vcrelax.out > scf.in && mpirun pw.x -nk 4 -in scf.in > scf.out && rm -rf tmp
S
 [ -z "$DRY" ] && (cd $W && sbatch run.sh)
done
