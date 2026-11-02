# PEW-C + modèle Z — premier exécutable réduit

Date : 2026-11-01

## Fonctionnelle statique implémentée

À partir du socle PEW v2.0, en métrique spatiale plate pour ce test (gravité non variée) :

\[
E=\int d^3x\left[\alpha'|\nabla\varrho|^2+\varrho^2\left(\frac14|\nabla n|^2+|\nabla\chi_g+\tfrac12 A|^2\right)+\frac{1}{4e_0^2}H_{ij}H_{ij}+V(\varrho)-V_{\rm vac}\right],
\]
\[
V(\varrho)=\lambda/\varrho^2+\mu^2\varrho^2+\sigma,\quad
\varrho_{\rm vac}=(\lambda/\mu^2)^{1/4},\quad V_{\rm vac}=2\mu\sqrt{\lambda}+\sigma.
\]
La constante \(\sigma\) disparaît de l'énergie localisée après soustraction du vide. La connexion est calculée via un relèvement spinoriel lisse du champ de Hopf, avec \(a=A/2=\operatorname{Im}(z^\dagger\nabla z)\), afin d'éviter la singularité de la carte azimutale.

## Périmètre du premier programme

- texture de Hopf de charge visée \(|Q_H|=1\), prescrite et non relaxée ;
- amplitude \(\varrho\) minimisée numériquement, avec bord fixé au vide ;
- phase globale restreinte à \(\chi_g=0\) ;
- taille globale \(L\) balayée sur une petite liste de valeurs ;
- paramètres exploratoires sans dimension : \(\alpha'=\lambda=\mu=e_0=1\).

Ce programme n'est donc pas encore la minimisation complète sur \((\varrho,n,\chi_g)\), et ne peut pas démontrer la stabilité d'un Hopfion libre. Il teste uniquement si l'amplitude se relaxe de façon régulière dans une texture prescrite, et fournit un premier banc de calcul pour étendre ensuite les degrés de liberté.

## Exécution N=41, boîte [-5,5]^3

| L | E amplitude + termes à deux dérivées | E Faddeev | E totale réduite | min(rho) | convergence de l'optimiseur rho |
|---:|---:|---:|---:|---:|---|
| 0.4 | 27.7579 | 180.9618 | 208.7196 | 0.5158 | oui |
| 0.5 | 36.3386 | 184.7566 | 221.0952 | 0.5164 | oui |
| 0.6 | 44.8334 | 178.4607 | 223.2941 | 0.5263 | oui |
| 0.8 | 61.1678 | 157.0185 | 218.1863 | 0.5580 | oui |
| 1.0 | 76.2827 | 135.9174 | 212.2001 | 0.5938 | oui |

## Interprétation limitée

L'optimiseur d'amplitude converge sur chaque échelle testée et garde \(\varrho_{\min}>0.51\) pour ces paramètres exploratoires. Mais les énergies ne présentent pas encore de minimum d'échelle fiable : la texture de petite taille est sous-résolue (à \(L=0.4\), le pas de grille vaut 0.25, donc le cœur ne couvre que quelques cellules), et le terme de Faddeev montre des erreurs de discrétisation visibles. Le plus bas point échantillonné, \(L=0.4\), ne doit pas être interprété comme un rayon optimal.

## Étapes indispensables suivantes

1. Vérifier numériquement la charge de Hopf de la texture discrète et sa convergence en résolution.
2. Faire converger les profils sur des grilles emboîtées, avec \(\Delta x/L\) assez petit et boîte agrandie.
3. Libérer \(\chi_g\) et minimiser conjointement la phase et l'amplitude.
4. Libérer l'orientation \(n\) avec une paramétrisation qui conserve \(|n|=1\), puis vérifier si la charge est préservée ou si un passage par \(\varrho=0\) la détruit.
5. Seulement après ces contrôles, tester le viriel et la stabilité sous perturbations non axisymétriques.
