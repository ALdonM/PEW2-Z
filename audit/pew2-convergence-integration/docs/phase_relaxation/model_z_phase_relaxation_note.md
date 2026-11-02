# PEW-C + modèle Z — première relaxation de phase
Date : 2026-11-02

## Test exécuté
Orientation de Hopf standard fixée, amplitude tenue à la valeur du vide rho=1 pour isoler l'équation de phase. Résolution discrète de div(grad chi_g + a)=0 avec chi_g=0 au bord, sur N=25, boîte [-5,5]^3, L=1.25. Paramètres sans dimension.

## Résultats
- Énergie du terme de phase avant relaxation : 35.4334646
- Énergie du terme de phase après relaxation : 29.9661136
- Réduction : 5.4673510, soit environ 15.43 %
- min/max chi_g : -0.3509116 / +0.3509116

## Interprétation et limites
La phase nulle n'est pas le minimum de ce terme pour cette discrétisation et ces conditions aux limites : laisser varier chi_g réduit l'énergie de phase d'environ 15.4 %. C'est un diagnostic cohérent avec l'équation elliptique de phase.

Ce test n'est PAS encore la minimisation couplée : rho a été maintenue à 1 pour le calcul de phase et l'orientation n reste prescrite. Le solveur d'amplitude exécuté en parallèle dans le script est un diagnostic de référence séparé, pas utilisé dans le calcul de la phase. La prochaine étape est de résoudre le problème pondéré avec rho variable, puis de tester des perturbations locales de n en réoptimisant les champs auxiliaires.
