#!/bin/bash
#SBATCH -t 0:10:00
#SBATCH --mem=230G
#SBATCH --exclusive
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH -x n2dgx01,n2gpu12[01-32],n2hcn[01-05],n2hacc[01-03],n2cn[1101-1196]

#module reset
#module load perf/OSU-Micro-Benchmarks/5.7.1-gompi-2021b
#srun --ntasks-per-node=1 -n $SLURM_JOB_NUM_NODES osu_mbw_mr  > osu.out

module reset
module load toolchain/foss/2023b
module load lang/Python/3.11.5-GCCcore-13.2.0
module load lib/zstd/1.5.5-GCCcore-13.2.0
source ./pqdts_env/bin/activate

make clean
make pqdts
make pqdts_mpi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=true

D=`echo "sqrt(200000000000/8/6/151*$SLURM_JOB_NUM_NODES)" | bc`
#eps=`echo "0.000001/$D" | bc -l`
eps=`echo "0.0000001/$D" | bc -l`
#eps=1
dir="simulated_D_260001_Nodes_${SLURM_JOB_NUM_NODES}_Dmax_${D}_$1"
rm -rf $dir
mkdir $dir
cd $dir

cp ../pqdts_omp_mpi.x .
python ../pqdts.py -P ../simulatedP_D_260001.npz -p ./pqdts_omp_mpi.x -v -b -d -D $D -e $eps > prep.out 2> prep.err

cmd=`tail -n1 prep.out`
wd=`pwd`
srun --ntasks-per-node=1 -n $SLURM_JOB_NUM_NODES cp pqdts_omp_mpi.x  /dev/shm
srun --ntasks-per-node=1 -n $SLURM_JOB_NUM_NODES cp -r data /dev/shm

cd /dev/shm
srun $cmd > $wd/run.out 2> $wd/run.err

#compress output
srun --ntasks-per-node=1 -n $SLURM_JOB_NUM_NODES bash $wd/../compress_and_copy.sh $wd

