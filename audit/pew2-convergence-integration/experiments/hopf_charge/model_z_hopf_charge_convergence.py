#!/usr/bin/env python3
"""Discrete Hopf-invariant diagnostic for the prescribed PEW/model-Z texture.

Convention: A = 2 a, H = dA, B = curl(A). Then
    Q_H = (1/(16 pi^2)) integral A . B d^3x
        = (1/(4 pi^2)) integral a . curl(a) d^3x.
This checks the background texture only; it does not prove that a relaxed,
fully dynamical solution exists.
"""
import argparse, json
import numpy as np
from model_z_first_test import derivatives


def hopf_charge(N, box, scale=1.0):
    x = np.linspace(-box, box, N)
    h = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = (X*X + Y*Y + Z*Z)/(scale*scale)
    den = 1.0 + r2
    z1 = (2*X/scale + 1j*2*Y/scale)/den
    z2 = (2*Z/scale + 1j*(r2-1.0))/den
    norm = np.sqrt(np.abs(z1)**2 + np.abs(z2)**2)
    z1, z2 = z1/norm, z2/norm
    dz1 = [derivatives(z1, h)[i] for i in range(3)]
    dz2 = [derivatives(z2, h)[i] for i in range(3)]
    a = [np.imag(np.conj(z1)*dz1[i] + np.conj(z2)*dz2[i]) for i in range(3)]
    curl_a = [
        derivatives(a[2], h)[1] - derivatives(a[1], h)[2],
        derivatives(a[0], h)[2] - derivatives(a[2], h)[0],
        derivatives(a[1], h)[0] - derivatives(a[0], h)[1],
    ]
    integrand = sum(a[i]*curl_a[i] for i in range(3))
    integral = float(np.sum(integrand)*h**3)
    return {'N': int(N), 'box': float(box), 'scale': float(scale), 'dx': float(h),
            'integral_a_curl_a': integral,
            'Q_H_discrete': integral/(4*np.pi**2),
            'Q_H_error_vs_1': abs(integral/(4*np.pi**2)-1.0)}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--cases', nargs='*', default=['31,5', '41,5', '61,5', '81,5', '65,4', '97,6', '113,7', '129,8'])
    p.add_argument('--scale', type=float, default=1.0)
    p.add_argument('--output', default='model_z_hopf_charge_convergence.json')
    args = p.parse_args()
    results = [hopf_charge(int(pair.split(',')[0]), float(pair.split(',')[1]), args.scale) for pair in args.cases]
    out = {'diagnostic': 'discrete Hopf invariant of prescribed smooth-spinor background',
           'convention': 'Q_H=(1/(16*pi^2))*integral A dot curl(A) d^3x, A=2a',
           'scope_warning': 'This verifies the background texture only; not a relaxed soliton or full model-Z minimizer.',
           'results': results}
    with open(args.output, 'w') as f:
        json.dump(out, f, indent=2)
    for r in results:
        print(f"N={r['N']:3d} box={r['box']:g} dx={r['dx']:.4f} Q_H={r['Q_H_discrete']:.8f} error={r['Q_H_error_vs_1']:.8f}")
    print('WROTE', args.output)

if __name__ == '__main__':
    main()
