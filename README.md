# Parallel Quantum Detector Tomography Solver (pqdts)
## Problem Statement
This program solves the detector tomography problem for a phase-insensitive detector
$$\mathrm{min}\quad \|\|P-F\varPi\|\|^2_{2} + \gamma \sum_{n=0}^{N-1}\sum_{i=0}^{M-1}(\varPi_{i,n}-\varPi_{i+1,n})^2$$
with 

$$ \varPi \in \mathbb{R}^{M\times N} $$

$$ P \in \mathbb{R}^{D\times N} $$

$$ F \in \mathbb{R}^{D\times M} $$

$$\gamma \ge 0$$

and the conditions

$$\varPi_{i,n}\ge 0 \ \forall i \in \{0,...,M-1\}, \ n\in \{0,...,N-1\}$$

$$\sum_{n=0}^{N-1}\varPi_{i,n}=1 \ \forall i \in \{0,...,M-1\}$$

It uses a two-stage approach based on a two-metric projected truncated Newton method and is parallelized with MPI and OpenMP. Scalability to hundreds of compute nodes and more than 100.000 CPU cores has been successfully demonstrated. Details are documented in the corresponding publication.

## How to Build
### Requirements
* a C compiler (e.g. gcc) and a Fortran compiler (e.g. gfortan) supporting OpenMP
* for the MPI version additionally an MPI implementation, e.g., OpenMPI

### Build Instructions
* adjust the build options in the `Makefile` according to your system
* run `make pqdts` for the OpenMP-parallelized version
* run `make pqdts_mpi` for the MPI-and-OpenMP-parallelized version

## How to Use
### Python Wrapper
The most convenient way of using the program is with the included Python wrapper `pqdts.py` which supports input matrices ($P$,$F$) from numpy (dense, npz) and scipy (sparse, npy). The Python wrapper only supports the OpenMP-parallelized version because running the MPI-parallelized version usually requires an adaption to the HPC system where is should run, e.g., the MPI-wrapper of the system like `mpirun`, `mpiexec`, `srun` or others.

This wrapper has the following options:

```
usage: pqdts.py [-h] -P PMATRIX [-F FMATRIX] [-D DMAX] [-M MMAX] [-i INITIAL] [-t THREADS] -p PQDTSPATH [-o OUTPUT] [-e EPSILON] [-g GAMMA] [-m MAXITER] [-T] [-b] [-d] [-v]

options:
  -h, --help            show this help message and exit
  -P, --Pmatrix PMATRIX
                        path to npz file (scipy sparse) or npy file (numpy) of P matrix (dimension D x N)
  -F, --Fmatrix FMATRIX
                        path to npz file (scipy sparse) or npy file (numpy) of F matrix (dimension D x M)
  -D, --Dmax DMAX       truncate to D so that P is a D x N matrix, positive D include D low indices ([:D]), negative D include D high indices [-D:]
  -M, --Mmax MMAX       truncate to M so that F is a D x M matrix, positive M include M low indices ([:M]), negative M include M high indices [-M:]
  -i, --initial INITIAL
                        path to npz file (scipy sparse) or npy file (numpy) of initial povm matrix (dimension M x N)
  -t, --threads THREADS
                        numper of OpenMP threads to use
  -p, --pqdtspath PQDTSPATH
                        path to compiled pqdts_omp.x
  -o, --output OUTPUT   output file for povm as npy
  -e, --epsilon EPSILON
                        convergence parameter of minimization
  -g, --gamma GAMMA     regularization parameter
  -m, --maxiter MAXITER
                        maximal number of iterations
  -T, --timing          measure timing for reconstruction, don't write output POVMs
  -b, --benchmarkops    measure timing for underlying operations
  -d, --dryrun          dry-run: only prepare inputs for pqdts
  -v, --verbose         be more verbose
```

The dependencies of the wrapper can be installed with `pip3 install -r requirements.txt`.

### Command Line Arguments
Without the Python wrapper, the command line arguments of `pqdts_omp.x` and `pqdts_omp_mpi.x` are:
1. size of first dimension of $\varPi$
2. size of second dimension of $\varPi$
3. size of first dimension of $P$
4. number of non-zero elements in $F$ or -1 to calculate $F$ internally
5. number of non-zero elements in $P$
6. computation mode: 2 for two-metric projected truncated Newton method, 1 for projected gradient method 
7. maxiter: maximal number of iterations in stages
8. output: 0 to disable output of the POVMs, 1 to enable output of POVMs at the end, 2 to enable output of POVMs after every minimization stage
9. gamma: regularization parameter $\gamma$
10. epsilon: value of the convergence criterion
11. index of stage to start with, i.e., 1 or 2
12. 0 to start with the initial guess $\varPi=1/N$, 1 to read the output from a previous run as a starting point
13. smoothing distance factor $N_s$
14. benchmark underlying operations: 0 for no, 1 for yes

### Inputs
Without the Python wrapper, the programs `pqdts_omp.x` and `pqdts_omp_mpi.x` expect inputs ($P$ matrix and optionally the $F$ matrix) in the `data` directory in the current working directory. The following files are expected:
* `data/P_row.bin` row indices of elements in $P$ in Fortran binary format
* `data/P_col.bin` column indices of elements in $P$ in Fortran binary format
* `data/P_data.bin` values of elements in $P$ in Fortran binary format
For $F$ in addition:
* `data/F_row.bin` row indices of elements in $F$ in Fortran binary format
* `data/F_col.bin` column indices of elements in $F$ in Fortran binary format
* `data/F_data.bin` values of elements in $F$ in Fortran binary format
For details have a look at the routine `read_sparse_from_python` in `pqdts.f90` and the Python wrapper `pqdts.py`.

For reading in existing POVMs as a starting point for the minimization:
* files of the form `"rank_"+rank+"_oiter"+oiter+".dat"` that were written by pqdts.

## Outputs
* stdout: progress information
* `"rank_"+rank+"_oiter"+oiter+".dat"`: files containing the POVMs in Fortran file format of rank `rank` after stage `oiter`
* `"timing_"+rank+".out"` file containing timing information of the seperate stages and timed routines
* `"status_"+rank+".out"` content of `/proc/self/status` at the end of the program execution to retreive memory usage from high-water mark

# Wigner Functions for High Photon Numbers
The computation of the Wigner function for very high photon numbers is not possible with conventional implementations because it requires the representability of very large floating-point numbers that exceed the range of 64-bit double-precision. To be able to use arbitrary-precision floating-point numbers, we make use of the multiple dispatch feature of the Julia programming language. With this, the Wigner-function implementation of QuantumOptics.jl (https://qojulia.org/, functions `wigner` and `_clenshaw` from https://github.com/qojulia/QuantumOptics.jl/blob/v1.0.15/src/phasespace.jl) can be generalized to arbitrary-precision floating-point numbers. Julia uses the GNU MPFR library for the BigFloat data type.

We have optimized the implementation for phase-insensitive detectors, i.e., diagonal density matrices and have trivially parallelized the computation because computations with arbitrary-precision floating-point numbers are much more cumbersome than with double-precision floating-point numbers.

## How to Build
A Julia installation is required to run the Wigner-programs. See https://julialang.org/downloads/ for details. The Julia packages that the program depends on can be installed with `julia --project=. -e 'using Pkg; Pkg.instantiate()'`.

## How to Use
### Variants
There are three variants:
* `wigner_fp64.jl`: using conventional 64-bit floating-point numbers, useful for exploring the need for larger arbitrary-precision floating-point numbers
* `wigner_mpfr.jl`: using arbitrary-precision floating-point numbers
* `wigner_mpfr_mpi.jl`: using arbitrary-precision floating-point numbers and a trivial MPI-parallelization to scale to many CPU-cores/compute nodes. This variant needs to be run as an MPI program, e.g, `mpirun julia --project=. wigner_mpfr_mpi.jl ...`. 

### Command Line Arguments
The three variants take the following command line arguments
1. input hdf5-file
2. n: name of dataset in hdf5 file 
3. number of points to evualate the wigner function on in positive x-direction
4. maximal x-argument for Wigner function
5. precision: for `wigner_mpfr.jl` and `wigner_mpfr_mpi.jl` the bit size of the mantissa

### Input HDF5-File
The input HDF5 file is expected to have a one-dimensional data set named `n` containing the POVM. It can be created with the following python script from $\varPi$ 

```python
import h5py
#povm (Pi) is a MxN numpy matrix 
mat=np.asarray(povm)
hf = h5py.File(filename, 'w')
for i in range(N):
    hf.create_dataset(str(i), data=mat[:,i])
hf.close()
```

### Output Files 
The scripts plot the Wigner functions and store the plot as `"n_"+(n)*"_M_"+(M)+"_precision_"+(precision)+".png"`. For plotting with Python the Wigner function is written as an HDF5-file `n_"+(n)+"_M_"+(M)+"_precision_"+(precision)+".h5", "w")` with group name `wigner` and data sets `x` for the points on the x-axis and `w` for the values of the Wigner function on these points.

This HDF5-file can be read in Python as follows:

```python
import h5py
hf = h5py.File(filename, 'r')
x=hf['wigner']["x"][:]
w=hf['wigner']["w"][:]
```
