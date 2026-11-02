# PEW-C + modèle Z — contrôle numérique intermédiaire de l'admission d'un soliton

Date : 2026-11-01

## Question examinée
Le modèle Z admet-il au moins un candidat de type Hopfion à énergie finie et à échelle stationnaire ? Le présent calcul n'examine qu'un sous-problème : texture de Hopf prescrite, amplitude \(\varrho\) relaxée, phase globale fixée \(\chi_g=0\), orientation \(n\) non relaxée, métrique plate.

## Fonctionnelle réduite
\[
E=\int d^3x\left[\alpha'|\nabla\varrho|^2+\varrho^2\left(\frac14|\nabla n|^2+|A/2|^2\right)+\frac{H_{ij}H_{ij}}{4e_0^2}+\frac{\lambda}{\varrho^2}+\mu^2\varrho^2-2\mu\sqrt{\lambda}\right],
\]
avec \(A=2a\), \(H=dA\), et les paramètres exploratoires \(\alpha'=\lambda=\mu=e_0=1\). Le bord impose \(\varrho=\varrho_{\rm vac}=1\). Ces valeurs sont sans dimension et ne sont pas des paramètres physiques dérivés.

## Contrôle de charge topologique
Convention numérique utilisée :
\[
Q_H=\frac{1}{16\pi^2}\int A\cdot(\nabla\times A)\,d^3x
=\frac{1}{4\pi^2}\int a\cdot(\nabla\times a)\,d^3x.
\]
À boîte demi-largeur 5 et échelle L=1 :

| N | dx | Q_H discret |
|---:|---:|---:|
| 31 | 0.3333 | 0.86726 |
| 41 | 0.2500 | 0.91902 |
| 61 | 0.1667 | 0.95881 |
| 81 | 0.1250 | 0.97331 |

À dx=0.125, lorsque la demi-largeur augmente de 4 à 8, Q_H passe de 0.96676 à 0.97885. Cela confirme une convergence vers la charge unitaire, mais la précision intégrale reste limitée par la discrétisation et les dérivées finies.

## Recherche d'une échelle stationnaire
Le balayage initial à boîte fixe était contaminé par une résolution relative variable. Un balayage plus cohérent a ensuite gardé la taille de boîte égale à 5L et le nombre de points constant, puis a été raffiné.

Énergies totales réduites obtenues :

| L | N=61 | N=81 |
|---:|---:|---:|
| 1.0 | 224.876 (optimiseur convergé) | 229.665 (non convergé en 50 itérations) |
| 1.1 | 220.746 | — |
| 1.2 | 218.908 | 223.099 |
| 1.25 | 218.676 | 222.756 |
| 1.3 | 218.832 | 222.812 |
| 1.4 | 220.137 | — |
| 1.5 | 222.546 | — |

Les valeurs d'énergie dérivent notablement avec la résolution : elles ne sont pas encore convergées en valeur absolue. Néanmoins, le minimum du balayage réduit se situe de façon assez stable autour de L≈1.25–1.30 sur N=61 et N=81. Une interpolation quadratique locale des trois points N=81 (L=1.20, 1.25, 1.30) place le minimum près de L≈1.27, à interpréter uniquement comme estimation exploratoire.

Le minimum est un indice d'équilibre contre une dilatation uniforme dans ce sous-problème. Ce n'est ni un test complet du viriel, ni une preuve de stabilité sous déformations générales.

## État de l'admission

**Résultat positif limité :** pour la texture de Hopf imposée, l'amplitude relaxée reste strictement positive sur les grilles testées (minimum environ 0.52–0.63 dans la zone des échelles proches du minimum), et l'énergie réduite présente un minimum d'échelle dans la plage examinée.

**Ce qui n'est pas démontré :**
- que \(n\) et \(\chi_g\) relaxés conservent une solution stationnaire ;
- que la charge de Hopf est préservée lorsque les champs sont entièrement libres, notamment si \(\varrho\) peut s'annuler ;
- que le candidat est stable face aux perturbations non axisymétriques ;
- que le minimum converge en énergie et en profil avec le raffinement ;
- que les paramètres exploratoires ont une signification physique.

## Conclusion
On ne peut pas encore conclure que le modèle Z complet admet un soliton stable. On dispose maintenant d'un **candidat stationnaire réduit** : texture topologique prescrite, amplitude libre, énergie finie sur le domaine numérique, et minimum d'échelle autour de L≈1.27. L'étape suivante décisive est une minimisation où l'orientation \(n\) varie sous la contrainte \(|n|=1\), puis où \(\chi_g\) est également libérée, tout en contrôlant la charge de Hopf, les conditions au bord et la convergence.
