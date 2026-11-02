# PEW-C + modèle Z — sonde de perturbation locale de l’orientation

## Méthode
La carte de Hopf est perturbée par des rotations SU(2) locales lisses, nulles au bord. La représentation spinorielle est renormalisée point par point, ce qui maintient exactement |z|=1 et donc |n|=1. Trois familles de perturbations sont testées (rotations autour des axes x, y, z), pour epsilon = -0.2, -0.1, -0.05, 0, 0.05, 0.1, 0.2.

Paramètres exploratoires sans dimension : N=25, boîte [-5,5]^3, L=1.25, alpha'=lambda=mu=e0=1. L’amplitude est fixée à rho=1 et chi_g=0. L’énergie inclut le gradient d’orientation, le terme de connexion de phase avec chi_g=0, et le terme de Faddeev.

## Limites
Ceci est une sonde directionnelle, pas une minimisation complète : rho et chi_g ne sont pas réoptimisés et seules trois familles de modes sont examinées. Une augmentation d’énergie ne démontre pas la stabilité générale ; une baisse persistante identifierait une direction défavorable à examiner ensuite avec relaxation couplée. Les effets de bord et de discrétisation restent à tester.

Les valeurs détaillées sont dans `model_z_orientation_probe_N25.json`.
