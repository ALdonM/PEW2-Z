#!/usr/bin/env python3
"""Phase relaxation for fixed Hopf background and prescribed rho.
Dimensionless diagnostic; orientation remains fixed. Solves div(rho^2(grad chi + a))=0
with chi=0 on boundary using a sparse finite-difference system.
"""
import json, argparse, time
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from model_z_first_test import hopf_background, optimize_amplitude, derivatives

def main():
 p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=25); p.add_argument('--box',type=float,default=5.0); p.add_argument('--L',type=float,default=1.25); p.add_argument('--maxiter',type=int,default=100); p.add_argument('--output',default='/mnt/data/model_z_phase_relaxation_N25.json'); args=p.parse_args()
 x=np.linspace(-args.box,args.box,args.N); h=x[1]-x[0]; N=args.N
 X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); r2=(X*X+Y*Y+Z*Z)/(args.L**2); den=1+r2
 z1=(2*X/args.L+1j*2*Y/args.L)/den; z2=(2*Z/args.L+1j*(r2-1))/den
 norm=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=norm; z2/=norm
 a=[np.imag(np.conj(z1)*derivatives(z1,h)[i]+np.conj(z2)*derivatives(z2,h)[i]) for i in range(3)]
 # Get rho from fixed-background amplitude minimization; then relax chi for this rho.
 _,U2,H2=hopf_background(N,args.box,args.L)
 amp=optimize_amplitude(U2,H2,h,maxiter=args.maxiter)
 # Reconstruct rho by a second solve is intentionally avoided; use rho=1 for this isolated phase-only first diagnostic.
 # This makes the phase PDE unweighted and provides a clean baseline before coupled relaxation.
 rho2=np.ones((N,N,N))
 shape=(N-2,N-2,N-2); m=np.prod(shape); idx=np.arange(m).reshape(shape); M=lil_matrix((m,m)); rhs=np.zeros(m)
 # Discrete div(grad chi + a)=0, chi boundary=0.
 for i in range(1,N-1):
  for j in range(1,N-1):
   for k in range(1,N-1):
    q=idx[i-1,j-1,k-1]; M[q,q]=-6/(h*h)
    for ax,(di,dj,dk) in enumerate([(1,0,0),(0,1,0),(0,0,1)]):
     ip,jp,kp=i+di,j+dj,k+dk; im,jm,km=i-di,j-dj,k-dk
     rhs[q]-=(a[ax][ip,jp,kp]-a[ax][im,jm,km])/(2*h)
     for ii,jj,kk in [(ip,jp,kp),(im,jm,km)]:
      if 0<ii<N-1 and 0<jj<N-1 and 0<kk<N-1: M[q,idx[ii-1,jj-1,kk-1]]=1/(h*h)
 chi=np.zeros((N,N,N)); chi[1:-1,1:-1,1:-1]=spsolve(M.tocsr(),rhs).reshape(shape)
 gradchi=[derivatives(chi,h)[i] for i in range(3)]
 E_before=np.sum(sum(ai*ai for ai in a))*h**3
 E_after=np.sum(sum((gradchi[i]+a[i])**2 for i in range(3)))*h**3
 out={'test':'phase relaxation, fixed Hopf orientation; rho held at vacuum value 1 for this isolated PDE test','N':N,'box':args.box,'dx':h,'L':args.L,'phase_boundary':'chi=0 on boundary','phase_energy_before':float(E_before),'phase_energy_after':float(E_after),'energy_reduction':float(E_before-E_after),'chi_min':float(chi.min()),'chi_max':float(chi.max()),'amplitude_reference_solve':amp,'warning':'This is a phase-only diagnostic with rho=1, not yet the coupled minimization of rho, chi, and n.'}
 open(args.output,'w').write(json.dumps(out,indent=2)); print(json.dumps(out,indent=2))
if __name__=='__main__': main()
