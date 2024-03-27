import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import os
import time
import argparse
import subprocess
import math
from scipy.io import FortranFile
from matplotlib import pyplot
import pickle

def npz_to_fortran(npz,matname):
    try:
        os.mkdir("data")
    except OSError as error:
        pass
    npz.todense().tofile('data/'+matname+'.bin')
    coo=sp.sparse.coo_array(npz)
    coorow=coo.row
    coocol=coo.col
    coodata=coo.data
    coorow.tofile('data/'+matname+'_row.bin')
    coocol.tofile('data/'+matname+'_col.bin')
    coodata.tofile('data/'+matname+'_data.bin')
    return len(coo.data)

def run_pqdts(N,M,D,P,threads,pqdtspath,maxiter=100,tol=1e-6,gamma=0,smo_int=0,smo_dist=0,start_stage=1,read_input=0,smo_fact=0,output=1,verbose=False,F=None):
    #set environment variables for OpenMP
    os.environ["OMP_NUM_THREADS"] = str(threads)
    os.environ["OMP_PROC_BIND"] = "True"
    os.environ["OMP_PLACES"] = "cores"
    #write input for Fortran code
    NP=npz_to_fortran(P,"P")
    NF=-1
    if not (F is None):
        NF=npz_to_fortran(F,"F")
    #call pqdts
    cmd=pqdtspath+" "+str(M)+" "+str(N)+" "+str(D)+" "+str(NF)+" "+str(NP)+" 2 "+str(maxiter)+" "+str(output)+" "+str(gamma)+" "+str(tol)+" "+str(start_stage)+" "+str(read_input)+" "+str(smo_fact)
    if verbose:
        print("executing:",cmd)
    
    start = time.time()
    out=subprocess.run(cmd, shell=True, capture_output=True)
    end = time.time()
    if verbose:
        for o in out.stdout.decode('utf-8').splitlines():
            print(o)
    for o in out.stderr.decode('utf-8').splitlines():
        print(o)

    #get objective value
    O=0
    M=0
    for o in out.stdout.decode('utf-8').splitlines():
        if o.startswith(" iter   "):
            O=float(o.split()[4])
    #get memory usage from status
    with open('status_           0 .out', "r") as f:
        for o in f.readlines():
            if o.startswith(" VmHWM:"):
                M=float(o.split()[1])
    print("memory used=",M,"kB")
    print("wall time=",end-start,"s")
    print("objective value=",O)

    #read output povm
    f='rank_     0_oiter     2.dat'
    fi = FortranFile(f, 'r')
    A=fi.read_ints(np.int32)
    M=A[0]
    A=fi.read_ints(np.int32)
    N=A[0]
    A=fi.read_ints(np.int32)
    D=A[0]
    A=fi.read_ints(np.int32)
    ML=A[0]
    A=fi.read_ints(np.int32)
    rank=A[0]
    X=fi.read_reals(float).reshape((ML,N), order="F")
    #print(M,N,D,rank,ML,X.shape[1])
    return X

parser = argparse.ArgumentParser()
parser.add_argument("-P", "--Pmatrix", help="path to npz file (scipy sparse) or npy file (numpy) of P matrix (dimension D x N)",type=str,required=True)
parser.add_argument("-F", "--Fmatrix", help="path to npz file (scipy sparse) or npy file (numpy) of F matrix (dimension D x M)",type=str,required=False)
parser.add_argument("-D", "--Dmax", help="truncate to D so that P is a D x N matrix",type=int,required=False)
parser.add_argument("-t", "--threads", help="numper of OpenMP threads to use",type=int,required=False,default=1)
parser.add_argument("-p", "--pqdtspath", help="path to compiled pqdts_omp.x",type=str,required=True)
parser.add_argument("-o", "--output", help="output file for povm as pickle",type=str,default="povm.p")
parser.add_argument("-e", "--epsilon", help="convergence parameter of minimization",type=float,default=1e-6)
parser.add_argument("-g", "--gamma", help="regularization parameter",type=float,default=0)
parser.add_argument("-m", "--maxiter", help="maximal number of iterations",type=int,default=200)
parser.add_argument("-b", "--benchmark", help="benchmark mode, don't write output POVMs",action='store_true')
parser.add_argument("-v", "--verbose", help="be more verbose",action='store_true')
args = parser.parse_args()

if args.Pmatrix.endswith(".npz"):
    P = sp.sparse.load_npz(args.Pmatrix)[:,:]
if args.Pmatrix.endswith(".npy"):
    P = np.load(args.Pmatrix)
    P=sp.sparse.csr_matrix(P)


Dmax=P.shape[0]
N=P.shape[1]
print("input P: N=",N,"D=",Dmax)
D=Dmax
if not (args.Dmax is None):
    print("truncating to D=",args.Dmax)
    D=args.Dmax

P=P[-D:,:]
M=0

#read F if given as commandline argument
if not(args.Fmatrix is None):
    if args.Fmatrix.endswith(".npz"):
        F = sp.sparse.load_npz(args.Fmatrix)[:,:]
    if args.Fmatrix.endswith(".npy"):
        F = np.load(args.Fmatrix)
        F=sp.sparse.csr_matrix(F)
    print("shape of read in F",F.shape)
    M=F.shape[1]
    if not (args.Dmax is None):
        F=F[-D:,:]
    if F.shape[0]!=D:
        raise RuntimeError("P and F don't have matching dimensions.")
else:
    print("WARNING: using internal F with quadratic scaling of photon numbers")
    M=(D-1)*(D-1)

print("N=",N,"D=",D,"M=",M,"N*M=",N*M)

#call OpenMP-version of pqdts
output=1
if args.benchmark:
    output=0
if not(args.Fmatrix is None):
    povm=run_pqdts(N,M,D,P,threads=args.threads,pqdtspath=args.pqdtspath,maxiter=args.maxiter,tol=args.epsilon,gamma=args.gamma,verbose=args.verbose,output=output,F=F)
else:
    povm=run_pqdts(N,M,D,P,threads=args.threads,pqdtspath=args.pqdtspath,maxiter=args.maxiter,tol=args.epsilon,gamma=args.gamma,verbose=args.verbose,output=output)

#check constraints
x1=0
x2=0
for i in range(povm.shape[0]):
    x1=max(x1,abs(np.sum(povm[i,:])-1))
    x2=x2+(np.sum(povm[i,:])-1)**2

print("maximal absolute violation of sum-constraint=",x1)
print("l2-norm of violation of sum-constraint=",math.sqrt(x2))
print("minimum povm value=",np.min(povm))

#dump as pickle for convenient further analysis and plotting
pickle.dump(povm,open(args.output,"wb"))


#plot POVM
fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)
ax.set_xscale("linear")
for i in range(povm.shape[1]):
    ax.plot(np.arange(povm.shape[0]),povm[:,i],linewidth=0.25)
fig.savefig("result_lin.png",dpi=300)
pyplot.close(fig)
                          
pyplot.clf()
fig = pyplot.figure()
ax = fig.add_subplot(1,1,1)
ax.set_xscale("log")
for i in range(povm.shape[1]):
    ax.plot(np.arange(povm.shape[0]),povm[:,i],linewidth=0.25)
fig.savefig("result_log.png",dpi=300)
pyplot.close(fig)


