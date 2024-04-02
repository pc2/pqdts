
# Copyright (c) 2024: Paderborn University, Paderborn Center for Parallel Computing, Robert Schade
# 
# Based on the wigner and _clenshaw functions from Quantumoptics.jl
# (https://github.com/qojulia/QuantumOptics.jl/blob/v1.0.15/src/phasespace.jl)
# Copyright (c) 2016: Sebastian Kraemer.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

using LinearAlgebra
using HDF5
using Plots

#wigner and _clenshaw0 are adapted/optimized from https://github.com/qojulia/QuantumOptics.jl/blob/v1.0.15/src/phasespace.jl
function wigner(rho,N, x)
    _2α = x*sqrt(2)
    abs2_2α = abs2(_2α)
    coeff = _clenshaw0( abs2_2α, rho, N)
    return exp(-abs2_2α/2)/pi*coeff
end

function _clenshaw0(abs2_2α::Real, ρ, N::Integer)
    n = N
    if n==0
        return ρ[1]
    elseif n==1
        ϕ1 = -(1-abs2_2α)/sqrt(1)
        return ρ[1] + ρ[2]*ϕ1
    else
        f0 = sqrt(Float64((n-1)*(n-1)))
        f1 = sqrt(Float64((n)*n))
        f0_ = 1/f0
        f1_ = 1/f1
        b2 = 0.
        b1 = ρ[n+1]
        b0 = ρ[n] - (2*n-1-abs2_2α)*f1_*b1
        @inbounds for k=n-2:-1:1
            b1, b2 = b0, b1
            b0 = ρ[k+1] - (2*k+1-abs2_2α)*f0_*b1 - f0*f1_*b2
            f1, f1_ = f0, f0_
            f0 = sqrt((k)*k)
            f0_ = 1/f0
        end
        return ρ[1] - (1-abs2_2α)*f0_*b0 - f0*f1_*b1
    end
end

function main()
  # wigner.jl data.h5 50 20000 2000 64
  d=ARGS[1]
  n=parse(Int,ARGS[2])
  steps=parse(Int,ARGS[3])
  xmax=parse(Int,ARGS[4])
  #prec=parse(Int,ARGS[5])
  prec="fp64"
  
  fid = h5open(d, "r")
  dset = fid[string(n)]
  rho = read(dset)
  
  M=length(rho)
  println("M=",M)


  x = collect(range(0, xmax, length=steps))
  #setprecision(prec)
  w = zeros(Float64,steps)
  wf = zeros(Float64,steps)


  for j in 1:steps
    w[j]=wigner(rho,M-1,Float64(x[j]))
    wf[j]=Float64(w[j])
    #println(j," ",x[j]," ",w[j]," ",precision(w[j]))
  end

  plot(x, wf)
  fid = h5open("n_"*string(n)*"_M_"*string(M)*"_precision_"*string(prec)*".h5", "w")
  gid=create_group(fid, "wigner")
  gid["x"]=x
  gid["w"]=wf
  close(fid)
  savefig("n_"*string(n)*"_M_"*string(M)*"_precision_"*string(prec)*".png")   
end
main()
