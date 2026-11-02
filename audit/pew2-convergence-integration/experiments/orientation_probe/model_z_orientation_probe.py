#!/usr/bin/env python3
"""Restricted local-orientation probe for PEW-C/model Z.

Perturbs the Hopf spinor by smooth, boundary-vanishing local SU(2) rotations.
This preserves |z|=1 and hence |n|=1 exactly. rho is held at 1 and chi_g=0
in this first probe: it is a directional energy test, NOT a full minimization.
"""
import json, time
import numpy as np
from pathlib import Path

def grad(f,h): return np.gradient(f,h,edge_order=2)
def energy(z1,z2,h):
    n1=2*np.real(np.conj(z1)*z2); n2=2*np.imag(np.conj(z1)*z2); n3=np.abs(z1)**2-np.abs(z2)**2
    ns=(n1,n2,n3)
    gradn2=sum(grad(ns[a],h)[i]**2 for a in range(3) for i in range(3))
    avec=[np.imag(np.conj(z1)*grad(z1,h)[i]+np.conj(z2)*grad(z2,h)[i]) for i in range(3)]
    A=[2*a for a in avec]
    H=np.zeros((3,3)+z1.shape)
    for i in range(3):
        for j in range(3): H[i,j]=grad(A[j],h)[i]-grad(A[i],h)[j]
    H2=np.sum(H*H,axis=(0,1))
    terms={'orientation_gradient':float(np.sum(0.25*gradn2)*h**3),
           'phase_connection_chi_zero':float(np.sum(sum(a*a for a in avec))*h**3),
           'faddeev':float(np.sum(H2/4)*h**3)}
    terms['total']=sum(terms.values())
    return terms

def main():
    N=25; box=5.; L=1.25; epsilons=[-0.2,-0.1,-0.05,0.,0.05,0.1,0.2]
    x=np.linspace(-box,box,N); h=x[1]-x[0]
    X,Y,Z=np.meshgrid(x,x,x,indexing='ij'); r2=(X*X+Y*Y+Z*Z)/L**2; den=1+r2
    z1=(2*X/L+1j*2*Y/L)/den; z2=(2*Z/L+1j*(r2-1))/den
    norm=np.sqrt(np.abs(z1)**2+np.abs(z2)**2); z1/=norm; z2/=norm
    base=energy(z1,z2,h); modes=[]
    # smooth local rotations about x, y, z axes; zero on all six faces
    envelope=np.sin(np.pi*(X+box)/(2*box))*np.sin(np.pi*(Y+box)/(2*box))*np.sin(np.pi*(Z+box)/(2*box))
    for axis in range(3):
      rows=[]
      for eps in epsilons:
        wx=np.zeros_like(X); wy=np.zeros_like(X); wz=np.zeros_like(X)
        if axis == 0: wx = eps*envelope
        elif axis == 1: wy = eps*envelope
        else: wz = eps*envelope
        w=np.sqrt(wx*wx+wy*wy+wz*wz); half=w/2
        c=np.cos(half); s=np.zeros_like(w); mask=w>1e-14; s[mask]=np.sin(half[mask])/w[mask]
        # U = c I - i sin(|w|/2) (w.sigma)/|w|
        q1=(wz*z1+(wx-1j*wy)*z2)
        q2=((wx+1j*wy)*z1-wz*z2)
        zz1=c*z1-1j*s*q1; zz2=c*z2-1j*s*q2
        norm2=np.sqrt(np.abs(zz1)**2+np.abs(zz2)**2); zz1/=norm2; zz2/=norm2
        e=energy(zz1,zz2,h)
        rows.append({'epsilon':eps,**e,'delta_total':e['total']-base['total']})
      modes.append({'mode':['local rotation about x','local rotation about y','local rotation about z'][axis], 'results':rows})
    out={'test':'restricted local orientation probe; fixed rho=1, chi_g=0; orientation varied via smooth SU(2) rotations', 'N':N,'box':box,'dx':h,'L':L,'parameters_dimensionless':{'alpha_prime':1,'lambda':1,'mu':1,'e0':1},'base_energy':base,'modes':modes,'interpretation':'Directional probe only. No amplitude/phase relaxation; selected smooth modes only; not a full stability proof. Boundary-vanishing rotations preserve the continuum Hopf class.', 'elapsed_seconds':None}
    path=Path('/mnt/data/model_z_orientation_probe_N25.json'); path.write_text(json.dumps(out,indent=2))
    note='''# PEW-C + modèle Z — sonde de perturbation locale de l’orientation\n\n## Méthode\nLa carte de Hopf est perturbée par des rotations SU(2) locales lisses, nulles au bord. La représentation spinorielle est renormalisée point par point, ce qui maintient exactement |z|=1 et donc |n|=1. Trois familles de perturbations sont testées (rotations autour des axes x, y, z), pour epsilon = -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2.\n\nParamètres exploratoires sans dimension : N=25, boîte [-5,5]^3, L=1.25, alpha'=lambda=mu=e0=1. L’amplitude est fixée à rho=1 et chi_g=0. L’énergie inclut le gradient d’orientation, le terme de connexion de phase avec chi_g=0, et le terme de Faddeev.\n\n## Limites\nCeci est une sonde directionnelle, pas une minimisation complète : rho et chi_g ne sont pas réoptimisés et seules trois familles de modes sont examinées. Une augmentation d’énergie ne démontre pas la stabilité générale ; une baisse persistante identifierait une direction défavorable à examiner ensuite avec relaxation couplée. Les effets de bord et de discrétisation restent à tester.\n\nLes valeurs détaillées sont dans `model_z_orientation_probe_N25.json`.\n'''
    Path('/mnt/data/model_z_orientation_probe_note.md').write_text(note)
    print(json.dumps({'base_energy':base,'modes':[{ 'mode':m['mode'],'deltas':[(r['epsilon'],round(r['delta_total'],6)) for r in m['results']]} for m in modes]},indent=2))
if __name__=='__main__': main()
