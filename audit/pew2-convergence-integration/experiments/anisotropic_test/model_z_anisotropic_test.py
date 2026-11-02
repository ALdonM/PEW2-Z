#!/usr/bin/env python3
"""Restricted anisotropic deformation test of the prescribed Hopf texture.

The Hopf map is pulled back by x_i -> x_i/s_i, with s=(L exp(t),L exp(-t),L).
This preserves the topology of the smooth map in the continuum. For each t the
amplitude rho is relaxed, while chi_g=0 and the spinor texture itself is not
otherwise varied. This is a diagnostic against a few shape deformations, not
full minimization of n.
"""
import argparse, json, time
import numpy as np
from scipy.optimize import minimize

def d(f,h,axis): return np.gradient(f,h,axis=axis,edge_order=2)

def fields(N,box,L,t):
    x=np.linspace(-box,box,N); h=x[1]-x[0]
    X,Y,Z=np.meshgrid(x,x,x,indexing='ij')
    sx,sy,sz=L*np.exp(t),L*np.exp(-t),L
    xx,yy,zz=X/sx,Y/sy,Z/sz
    r2=xx*xx+yy*yy+zz*zz; den=1+r2
    z1=(2*xx+2j*yy)/den
    z2=(2*zz+1j*(r2-1))/den
    norm=np.sqrt(abs(z1)**2+abs(z2)**2); z1/=norm; z2/=norm
    n=(2*np.real(np.conj(z1)*z2),2*np.imag(np.conj(z1)*z2),abs(z1)**2-abs(z2)**2)
    gradn2=sum(d(n[a],h,i)**2 for a in range(3) for i in range(3))
    a=[np.imag(np.conj(z1)*d(z1,h,i)+np.conj(z2)*d(z2,h,i)) for i in range(3)]
    A=[2*v for v in a]
    H=np.zeros((3,3,N,N,N))
    for i in range(3):
      for j in range(3): H[i,j]=d(A[j],h,i)-d(A[i],h,j)
    return h,0.25*gradn2+sum(v*v for v in a),np.sum(H*H,axis=(0,1))

def relax(U2,H2,h,alpha=1.,lam=1.,mu=1.,e0=1.,maxiter=160):
    N=U2.shape[0]; r0=np.ones((N,N,N)); sl=(slice(1,-1),)*3; shape=r0[sl].shape; dv=h**3
    def full(x):
      r=r0.copy(); r[sl]=x.reshape(shape); return r
    def fg(x):
      r=full(x); Eg=alpha*sum(np.sum(np.diff(r,axis=i)**2)/h**2*dv for i in range(3))
      V=lam/r**2+mu**2*r**2-2*mu*np.sqrt(lam)
      E=Eg+np.sum(U2*r*r+V)*dv
      lp=np.zeros_like(r)
      for ax in range(3):
        diff=np.diff(r,axis=ax); lo=[slice(None)]*3; hi=lo.copy(); lo[ax]=slice(0,-1); hi[ax]=slice(1,None)
        lp[tuple(lo)]-=diff/h**2; lp[tuple(hi)]+=diff/h**2
      g=2*alpha*lp+2*(U2*r+mu**2*r-lam/r**3)
      return E,(g[sl]*dv).ravel()
    res=minimize(fg,np.ones(np.prod(shape)),method='L-BFGS-B',jac=True,bounds=[(.08,3.)]*int(np.prod(shape)),options={'maxiter':maxiter,'ftol':1e-9,'gtol':2e-5,'maxls':20})
    E4=np.sum(H2/(4*e0**2))*dv
    return {'E2_amp':float(res.fun),'E4_Faddeev':float(E4),'E_total':float(res.fun+E4),'rho_min':float(np.min(full(res.x))),'success':bool(res.success),'iterations':int(res.nit),'message':str(res.message),'grad_inf':float(np.max(abs(res.jac)))}

def main():
 p=argparse.ArgumentParser(); p.add_argument('--N',type=int,default=31); p.add_argument('--box',type=float,default=5.); p.add_argument('--L',type=float,default=1.25); p.add_argument('--ts',type=float,nargs='+',default=[-.3,-.2,-.1,0.,.1,.2,.3]); p.add_argument('--maxiter',type=int,default=160); p.add_argument('--output',default='model_z_anisotropic_test_N31.json'); a=p.parse_args()
 rows=[]
 for t in a.ts:
  st=time.time(); h,U,H=fields(a.N,a.box,a.L,t); r=relax(U,H,h,maxiter=a.maxiter); r.update({'t':t,'L':a.L,'N':a.N,'box':a.box,'dx':h,'elapsed_s':round(time.time()-st,2)}); rows.append(r); print(json.dumps(rows[-1]))
 out={'scope':'anisotropic pullback of prescribed Hopf map; rho relaxed; chi_g fixed zero; n not freely varied','parameters':'dimensionless alpha=lambda=mu=e0=1','results':rows,'warning':'Restricted deformation test only; no full n or chi_g relaxation.'}
 open(a.output,'w').write(json.dumps(out,indent=2)); print('WROTE',a.output)
if __name__=='__main__': main()
