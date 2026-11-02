# PEW-C + modèle Z — archive de travail

Date: 2026-11-01

Cette archive regroupe les scripts, résultats JSON et notes produits durant les premiers tests numériques du modèle Z.

## Statut scientifique
- Test réduit : texture de Hopf prescrite, amplitude rho minimisée, phase globale chi_g fixée à zéro, orientation n non relaxée.
- Les paramètres sont exploratoires et sans dimension.
- La charge de Hopf discrète converge vers 1 avec le raffinement, mais reste affectée par les erreurs de discrétisation et de domaine fini.
- Le balayage réduit suggère un minimum en échelle vers L ~ 1.25–1.30, mais les énergies absolues ne sont pas convergées.
- Cela ne démontre pas l’existence ni la stabilité d’un soliton du modèle Z complet.

## Contenu
Tous les fichiers `model_z_*` présents au moment de la création de cette archive sont inclus. Les résultats ne doivent pas être interprétés comme des prédictions physiques.

## Étapes suivantes
1. Contrôler les conventions et la convergence topologique.
2. Améliorer les tests de résolution et de taille du domaine.
3. Libérer l’orientation n sous contrainte |n|=1.
4. Libérer chi_g, puis tester stabilité non axisymétrique et viriel.

Cette archive est un instantané de travail, pas une publication finale.
