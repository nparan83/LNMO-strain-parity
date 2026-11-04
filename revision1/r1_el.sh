#!/bin/bash
# r1_el.sh - CMS R1 madde 3: elastik sabitler (eps=0), 12 relax isi tek SLURM dizisi
P=/arf/scratch/nparan/lnmo/LNMO_Strain_Production; R=/arf/scratch/nparan/R1; W=$R/el
mkdir -p $W && python3 $R/el_gen.py $P/strain_0.00/FM/scf_FM.in $W || exit 1
grep -E "calculation|tstress|nstep" $W/1/pw.in; diff <(sed -n '/CELL_PARAMETERS/,+3p' $P/strain_0.00/FM/scf_FM.in) <(sed -n '/CELL_PARAMETERS/,+3p' $W/1/pw.in)
cat > $W/array.sh << 'S'
#!/bin/bash
#SBATCH -J elas
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=1-00:00:00
#SBATCH --array=1-12
#SBATCH --output=slurm_%A_%a.out
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023 lib/hdf5/1.14.3-oneapi-2023.0 apps/espresso/7.2-oneapi-2023.0
cd $SLURM_ARRAY_TASK_ID && mpirun pw.x -nk 4 -in pw.in > pw.out && rm -rf tmp
S
[ -z "$DRY" ] && (cd $W && sbatch array.sh)
