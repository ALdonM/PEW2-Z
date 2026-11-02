#!/usr/bin/env python3
"""PEW-C + model Z: first reduced numerical test.

Scope: fixed Q_H=1 Hopf texture, free amplitude rho, and scan of global
texture scale L. The global phase chi_g is fixed to zero and the Hopf
orientation is not yet varied. This is a diagnostic stage, NOT a full 3D
soliton stability calculation.

Static energy in the stated Z convention:
 E = int [ alpha' |grad rho|^2
         + rho^2 (|grad n|^2/4 + |a|^2)
         + H_ij H_ij/(4 e0^2)
         + V(rho)-V_vac ] d^3x,
 a=A/2, A=(1-n_3)dphi, H=dA.
The Hopf connection is evaluated through a smooth spinor lift z of the
Hopf map; a_i = Im(z^dagger partial_i z), avoiding a coordinate-string
singularity. chi_g=0 is an explicit restriction of this first test.
"""
import argparse, json, time
import numpy as np
from scipy.optimize import minimize


def derivatives(f, h):
    return np.gradient(f, h, edge_order=2)


def hopf_background(N=25, box=6.0, scale=1.0):
    x = np.linspace(-box, box, N)
    h = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = (X*X + Y*Y + Z*Z)/(scale*scale)
    den = 1.0 + r2
    # S^3 stereographic coordinates; |z1|^2+|z2|^2=1
    z1 = (2*X/scale + 1j*2*Y/scale) / den
    z2 = (2*Z/scale + 1j*(r2-1.0)) / den
    # Roundoff normalization
    norm = np.sqrt(np.abs(z1)**2 + np.abs(z2)**2)
    z1, z2 = z1/norm, z2/norm
    n1 = 2*np.real(np.conj(z1)*z2)
    n2 = 2*np.imag(np.conj(z1)*z2)
    n3 = np.abs(z1)**2 - np.abs(z2)**2
    n = (n1,n2,n3)
    dn = [[derivatives(n[a],h)[i] for i in range(3)] for a in range(3)]
    gradn2 = sum(dn[a][i]**2 for a in range(3) for i in range(3))
    # Smooth Hopf connection a = Im(z^dagger dz) = A/2
    dz1 = [derivatives(z1,h)[i] for i in range(3)]
    dz2 = [derivatives(z2,h)[i] for i in range(3)]
    a = [np.imag(np.conj(z1)*dz1[i] + np.conj(z2)*dz2[i]) for i in range(3)]
    A = [2*ai for ai in a]
    H = np.zeros((3,3,N,N,N), dtype=float)
    for i in range(3):
        for j in range(3):
            H[i,j] = derivatives(A[j],h)[i] - derivatives(A[i],h)[j]
    H2 = np.sum(H*H, axis=(0,1))
    U2 = 0.25*gradn2 + sum(ai*ai for ai in a)
    return h, U2, H2


def optimize_amplitude(U2, H2, h, alpha=1.0, lam=1.0, mu=1.0, e0=1.0,
                       maxiter=180):
    N = U2.shape[0]
    rho0 = np.ones((N,N,N), dtype=float)
    # Fix boundary to vacuum; optimize interior only.
    sl = (slice(1,-1),)*3
    Ui = U2[sl]
    rho_init = rho0[sl].copy()
    dv = h**3
    Vvac = 2.0*mu*np.sqrt(lam)

    def full(x):
        r = rho0.copy()
        r[sl] = x.reshape(rho_init.shape)
        return r

    def fg(x):
        r = full(x)
        # Forward-edge discretization of |grad rho|^2, with fixed boundary.
        Egrad = 0.0
        for ax in range(3):
            Egrad += alpha*np.sum(np.diff(r, axis=ax)**2)/(h*h)*dv
        V = lam/(r*r) + mu*mu*r*r - Vvac
        Eamp = Egrad + np.sum((U2*r*r + V))*dv
        # Exact discrete derivative of edge-gradient energy.
        lap_pos = np.zeros_like(r)
        for ax in range(3):
            diff = np.diff(r, axis=ax)
            lo = [slice(None)]*3; hi = [slice(None)]*3
            lo[ax] = slice(0,-1); hi[ax] = slice(1,None)
            lap_pos[tuple(lo)] -= diff/(h*h)
            lap_pos[tuple(hi)] += diff/(h*h)
        grad = 2*alpha*lap_pos + 2*(U2*r + mu*mu*r - lam/(r**3))
        gi = grad[sl]*dv
        return Eamp, gi.ravel()

    result = minimize(fg, rho_init.ravel(), method='L-BFGS-B', jac=True,
                      bounds=[(0.08, 3.0)]*rho_init.size,
                      options={'maxiter':maxiter, 'ftol':1e-9, 'gtol':1e-5,
                               'maxls':20})
    rho = full(result.x)
    E4 = np.sum(H2/(4*e0*e0))*dv
    Etotal = result.fun + E4
    return {'success': bool(result.success), 'message': str(result.message),
            'iterations': int(result.nit), 'nfev': int(result.nfev),
            'energy_amplitude_orientation_phase': float(result.fun),
            'energy_faddeev': float(E4), 'energy_total': float(Etotal),
            'rho_min': float(np.min(rho)), 'rho_max': float(np.max(rho)),
            'rho_boundary': float(np.mean(np.concatenate([rho[0].ravel(),rho[-1].ravel(),rho[:,0].ravel(),rho[:,-1].ravel(),rho[:,:,0].ravel(),rho[:,:,-1].ravel()]))),
            'gradient_inf_norm': float(np.max(np.abs(result.jac))),
            'vacuum_energy_density_subtracted': float(2*1.0*np.sqrt(1.0))}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--N', type=int, default=25)
    p.add_argument('--box', type=float, default=6.0)
    p.add_argument('--scales', type=float, nargs='+', default=[0.8,1.0,1.25,1.6])
    p.add_argument('--maxiter', type=int, default=180)
    p.add_argument('--alpha', type=float, default=1.0)
    p.add_argument('--lambda', dest='lam', type=float, default=1.0)
    p.add_argument('--mu', type=float, default=1.0)
    p.add_argument('--e0', type=float, default=1.0)
    p.add_argument('--output', default='model_z_first_test_results.json')
    args = p.parse_args()
    rows=[]
    for L in args.scales:
        t=time.time()
        h,U2,H2=hopf_background(args.N,args.box,L)
        res=optimize_amplitude(U2,H2,h,args.alpha,args.lam,args.mu,args.e0,args.maxiter)
        res.update({'scale_L':float(L),'N':args.N,'box':args.box,'dx':h,
                    'alpha_prime':args.alpha,'lambda':args.lam,'mu':args.mu,'e0':args.e0,
                    'chi_g_fixed':0.0,'orientation_varied':False,
                    'elapsed_seconds':round(time.time()-t,3)})
        rows.append(res)
        print(json.dumps(rows[-1], ensure_ascii=False))
    best=min(rows,key=lambda r:r['energy_total'])
    out={'model':'PEW-C + model Z, reduced first test',
         'scope':'fixed Hopf Q=1 texture; free amplitude; fixed chi_g=0; scale scan; orientation not varied',
         'parameters':'dimensionless exploratory values, not physical predictions',
         'results':rows,'lowest_sampled_scale':best['scale_L'],
         'lowest_sampled_energy':best['energy_total'],
         'warning':'A minimum among sampled scales is not a full 3D variational minimum or stability proof.'}
    with open(args.output,'w') as f: json.dump(out,f,indent=2)
    print('WROTE',args.output)

if __name__=='__main__': main()
