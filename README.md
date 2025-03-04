# Solveur parallèle modulaire pour l'équation de la chaleur en 1D

Ce code résout l'équation de la chaleur en une dimension spatiale avec la méthode d'Euler implicite, 
discrétisation spatiale par différences finies centrées d'ordre 2 et un solveur gradient conjugué
parallélisé avec MPI.

## Structure modulaire

Le code est organisé en plusieurs modules pour améliorer la lisibilité et la maintenance :

- **heat.f90** : Programme principal qui orchestre la résolution
- **parameters_mod.f90** : Définition des constantes et paramètres
- **domain_mod.f90** : Gestion du domaine spatial et décomposition parallèle
- **solver_mod.f90** : Méthodes numériques (Euler implicite, gradient conjugué)
- **output_mod.f90** : Fonctions d'entrée/sortie et visualisation

## Équation résolue

L'équation de la chaleur:
```
∂u/∂t = α ∂²u/∂x²
```

où:
- u(x,t) est la température
- α est la diffusivité thermique
- x ∈ [0,L] est la variable spatiale
- t est le temps

## Méthode numérique

- **Discrétisation temporelle**: Schéma d'Euler implicite
  ```
  (u^{n+1} - u^n)/Δt = α ∂²u^{n+1}/∂x²
  ```

- **Discrétisation spatiale**: Différences finies centrées d'ordre 2
  ```
  ∂²u/∂x² ≈ (u_{i+1} - 2u_i + u_{i-1})/(Δx)²
  ```

- **Système linéaire résultant**: Tridiagonal de la forme
  ```
  -λ u_{i-1}^{n+1} + (1+2λ) u_i^{n+1} - λ u_{i+1}^{n+1} = u_i^n
  ```
  où λ = α·Δt/(Δx)²

- **Résolution du système linéaire**: Méthode du gradient conjugué
  - Algorithme itératif pour systèmes symétriques définis positifs
  - Convergence rapide pour les matrices tridiagonales
  - Facilement parallélisable avec MPI

- **Parallélisation**: 
  - Décomposition du domaine spatial
  - Chaque processus gère une portion des points intérieurs
  - Communications MPI aux frontières des sous-domaines

## Compilation

```bash
make
```

Le compilateur et les flags utilisés sont définis dans le Makefile:
- Compilateur: mpif90 (wrapper Fortran pour MPI)
- Flags: -O3 -Wall -Wextra (optimisation et avertissements)

## Exécution

```bash
make run
```
ou
```bash
mpirun -np <nombre_de_processus> ./heat
```

Par défaut, la commande `make run` utilise 4 processus. Vous pouvez ajuster ce nombre en modifiant le Makefile.

## Paramètres du modèle

Les paramètres sont définis dans `parameters_mod.f90`:

- N = 100: nombre de points de la grille spatiale
- L = 1.0: longueur du domaine
- α = 1.0: diffusivité thermique
- dt = 0.001: pas de temps
- t_final = 1.0: temps de simulation final
- epsilon = 1.0e-6: tolérance pour la convergence du gradient conjugué

## Visualisation

Le script Python `plot_solution.py` offre plusieurs options pour visualiser les résultats:

```bash
# Afficher la dernière solution calculée
python3 plot_solution.py

# Afficher la solution à un pas de temps spécifique
python3 plot_solution.py -s 50

# Afficher plusieurs pas de temps sur le même graphique 
python3 plot_solution.py -m -n 100

# Créer une animation de l'évolution de la solution
python3 plot_solution.py -a

# Créer une visualisation en carte de chaleur spatio-temporelle
python3 plot_solution.py -c


```

## Conditions initiales et aux limites

- **Conditions aux limites**: Dirichlet (valeurs fixes aux extrémités)
  - u(0,t) = 0
  - u(L,t) = 0

- **Condition initiale**: 
  - u(x,0) = sin(πx) 

Ces conditions peuvent être modifiées dans les fichiers `domain_mod.f90` et `solver_mod.f90`.

## Performance

Le code est optimisé pour les architectures HPC avec:
- Décomposition de domaine équilibrée, même quand N-1 n'est pas divisible par le nombre de processus
- Communications MPI non bloquantes pour maximiser l'efficacité
- Stockage minimal pour réduire l'empreinte mémoire
- Gradient conjugué qui minimise le nombre d'itérations nécessaires

## Organisation des fichiers

Le code génère deux types de fichiers de sortie:
- **Données brutes**: Stockées dans le dossier `results/`
  - Format: `results/heat_solution_XXXXX.dat`
  - Contenu: position (x) et température (u) à un pas de temps donné

- **Visualisations**: Stockées dans le dossier `viz/`
  - Graphiques statiques (PNG)
  - Animations (GIF, MP4)
  - Cartes de chaleur (heatmaps)

## Limitations et extensions possibles

- Le code actuel utilise des conditions aux limites de Dirichlet fixes. Il pourrait être étendu pour supporter d'autres types (Neumann, Robin).
- Pour améliorer la performance, un préconditionnement du gradient conjugué pourrait être implémenté.
- Le domaine est actuellement 1D; l'extension à 2D ou 3D nécessiterait une refonte des schémas de communication.

## Installation

### Prérequis
- Python 3.6+
- Compilateur Fortran compatible MPI (mpif90)
- Une implémentation MPI (comme OpenMPI ou MPICH)

### Configuration de l'environnement Python
```bash
# Créer un environnement virtuel (recommandé)
python3 -m venv venv

# Activer l'environnement
source venv/bin/activate  # Sur Linux/Mac
# venv\Scripts\activate  # Sur Windows

# Installer les dépendances
pip install numpy matplotlib
```

## Contribuer

Les contributions sont les bienvenues ! Si vous souhaitez améliorer le code ou la documentation, n'hésitez pas à soumettre une pull request.

## Auteurs

- [AÏNA Yanel, AHOUANYE Elonm ()] - Auteurs principaux
- Autres contributeurs : ...

## Licence

Ce projet est sous la licence MIT - voir le fichier [LICENSE](LICENSE) pour plus de détails.




