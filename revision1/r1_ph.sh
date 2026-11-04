#!/bin/bash
# r1_ph.sh - CMS R1 madde 1: Gamma fononu, eps=-5,0,+5; her gerinim bir SLURM dizisi (d0 = mukemmel hucre)
P=/arf/scratch/nparan/lnmo/LNMO_Strain_Production; R=/arf/scratch/nparan/R1
L=${*:-"-0.05 0.00 0.05"}; for e in $L; do
 W=$R/ph/eps_$e; mkdir -p $W
 python3 $R/ph_gen.py $P/strain_$e/FM/scf_FM.in $W || exit 1
 N=$(cat $W/ndisp)
 cat > $W/array.sh << S
#!/bin/bash
#SBATCH -J ph$e
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=08:00:00
#SBATCH --array=0-$N
#SBATCH --output=slurm_%A_%a.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023 lib/hdf5/1.14.3-oneapi-2023.0 apps/espresso/7.2-oneapi-2023.0
cd d\$SLURM_ARRAY_TASK_ID && mpirun pw.x -nk 4 -in pw.in > pw.out && rm -rf tmp
S
 [ -z "$DRY" ] && (cd $W && sbatch array.sh)
done
