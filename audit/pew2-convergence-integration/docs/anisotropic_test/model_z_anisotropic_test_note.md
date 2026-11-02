# PEW-C + modèle Z — test de déformation anisotrope restreinte
Date : 2026-11-01

## Objet
Tester si le candidat de Hopf prescrit, après relaxation de l'amplitude, résiste à une famille de déformations de forme anisotropes. La texture est obtenue par tiré en arrière de la carte de Hopf standard sous

\[ (x,y,z)\mapsto (x/s_x,y/s_y,z/s_z),\quad (s_x,s_y,s_z)=(Le^t,Le^{-t},L). \]

Le produit \(s_xs_ys_z=L^3\) est conservé : on teste une déformation de forme à volume caractéristique constant. Dans le continuum, cette transformation est un difféomorphisme de l'espace et conserve la classe topologique de la carte.

## Conditions
- \(L=1.25\), boîte cubique \([-5,5]^3\).
- \(\alpha'=\lambda=\mu=e_0=1\), paramètres sans dimension exploratoires.
- \(\chi_g=0\), orientation décrite par la carte de Hopf déformée mais non libre point par point.
- \(\varrho\) relaxée par L-BFGS-B, bord fixé à la valeur du vide 1.
- \(t\in\{-0.2,-0.1,0,0.1,0.2\}\).

## Résultats

| N | dx | E(t=0) | E(t=0.1)-E(0) | E(t=0.2)-E(0) | rho_min à t=0 |
|---:|---:|---:|---:|---:|---:|
| 31 | 0.3333 | 198.4004 | 1.5475 | 6.1713 | 0.6541 |
| 41 | 0.2500 | 207.7118 | 1.8694 | 7.4814 | 0.6369 |
| 51 | 0.2000 | 212.3261 | 2.0434 | 8.1934 | 0.6286 |

Les énergies sont symétriques entre +t et -t à la précision numérique observée. La différence seconde centrée donne :

| N | courbure discrète à t=0.1 | courbure discrète à t=0.2 |
|---:|---:|---:|
| 31 | 309.51 | 308.56 |
| 41 | 373.87 | 374.07 |
| 51 | 408.69 | 409.67 |

## Interprétation
Dans cette famille à un paramètre, \(t=0\) est un minimum local strict : les déformations anisotropes testées augmentent l'énergie après relaxation de l'amplitude. C'est un résultat positif supplémentaire pour le candidat symétrique.

Cependant, la courbure discrète augmente avec N ; les énergies absolues et les profils ne sont donc pas convergés. Le test n'explore qu'une famille particulière de déformations et ne remplace pas la variation générale de \(n(\mathbf x)\). Il ne démontre ni stabilité non axisymétrique générale, ni minimum complet du modèle Z.

## Conclusion
Le candidat résiste à cette déformation anisotrope restreinte pour les trois grilles testées. La prochaine étape doit être une minimisation de l'orientation sur la sphère \(|n|=1\), puis de \(\chi_g\), avec surveillance explicite de la charge de Hopf et des modes instables.
