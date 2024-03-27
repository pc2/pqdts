#!/bin/bash
#SBATCH -t 1-0
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --exclusive
#SBATCH --cpus-per-task=128


module reset
module load toolchain/foss/2023b
module load lang/Python/3.11.5-GCCcore-13.2.0
source ./pqdts_env/bin/activate

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=true

D=$1
mkdir simulated_D_8200_Dmax_${D}
cd simulated_D_8200_Dmax_${D}
python3 ../pqdts.py -P ../../det-tomo/simulation/data/simulatedP_D_8200.npz -t $SLURM_CPUS_PER_TASK --pqdtspath ../pqdts_omp.x -m 1000 -g 0.0 -e 0.000000001 -D $D -v -b > simulated_D_8200_Dmax_${D}.out
