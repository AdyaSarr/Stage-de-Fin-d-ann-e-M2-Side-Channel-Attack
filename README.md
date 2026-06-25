# Stage M2 — Attaque par Canal Auxiliaire contre Classic McEliece

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![galois v0.4.11](https://img.shields.io/badge/galois-v0.4.11-blue.svg)](https://mhostetter.github.io/galois/latest/)

Ce dépôt contient l'implémentation complète d'une attaque par canal auxiliaire (*side-channel attack*) contre le cryptosystème post-quantique **Classic McEliece**. Le travail est mené dans le cadre d'un stage de fin d'études M2 au **Laboratoire Hubert Curien** (Saint-Étienne), sous l'encadrement de **Pierre-Louis Cayrel(Maître de conférences HDR UJM - LabHC), Brice Colombier(chercheur et professeur en electronique et programmation a IUT Saint-Etienne), Michaël Bulois(Maitre de Conférences en Mathématiques  au laboratoire Institut Camille Jordan), Vincent GROSSO (CR Researcher, CNRS) et Vlad-Florin Dragoi (Professor (Associate) at Aurel Vlaicu University of Ara in Roumania)**

## Contexte

Classic McEliece est un cryptosystème post-quantique basé sur les **codes de Goppa binaires**, finaliste du standard NIST PQC. Sa sécurité repose sur la difficulté de décoder un code aléatoire et la dissimulation de la structure de Goppa via une permutation secrète.

Ce travail étudie une **attaque par canal auxiliaire** exploitant des fuites lors de la décapsulation. À partir d'acquisitions bruitées de poids de Hamming :

$$y_{i,j} = \mathrm{wt}(\alpha_i^j \cdot G(\alpha_i)^{-2}) + \varepsilon_{i,j}, \qquad \varepsilon_{i,j} \sim \mathcal{N}(0, \sigma^2)$$

l'attaquant reconstruit intégralement la clé secrète $(G, L)$ : le polynôme de Goppa $G$ et le support $L = (\alpha_0, \ldots, \alpha_{n-1})$.

## Pipeline d'attaque

<table>
<tr>
<td><b>1. Décodage MLE par FFT</b></td>
<td>Acquisitions bruitées (n×N) → (αᵢ, βᵢ = G(αᵢ)⁻²)</td>
<td>O(P² log P), P = 2ᵐ − 1</td>
</tr>
<tr>
<td><b>2. Reconstruction algébrique de G</b></td>
<td>Bernstein (2024), interpolation Reed–Solomon avec correction d'erreurs</td>
<td>O(n²m²)</td>
</tr>
<tr>
<td><b>3. Reconstruction du support L</b></td>
<td>Pivot de Gauss sur la matrice de parité publique</td>
<td>O((mt)³ + n(mt)²)</td>
</tr>
</table>
## Résultats expérimentaux

### Validation sur paramètres officiels Classic McEliece

| Paramètre        | m  | t   | n    | mt   | σ testé | Attaque |
|------------------|----|-----|------|------|---------|---------|
| mceliece348864   | 12 | 64  | 3488 | 768  | 0 – 2.3 | ✓       |
| mceliece348864   | 12 | 64  | 3488 | 768  | 2.4+    | ✗       |
| mceliece460896   | 13 | 96  | 4608 | 1248 | 0 – 3.0+| ✓       |

L'attaque réussit complètement (reconstruction de $G$ **et** de $L$) tant que le taux d'erreur du décodeur FFT reste sous environ 5 % pour `mceliece348864` et sous 3 % pour `mceliece460896`.


### Balayage de σ sur mceliece348864

| σ    | p_succès | erreurs | rec. G | rec. L |
|------|----------|---------|--------|--------|
| 0.5  | 0.9969   | 3       | ✓      | ✓      |
| 1.0  | 0.9979   | 2       | ✓      | ✓      |
| 1.5  | 0.9886   | 11      | ✓      | ✓      |
| 2.0  | 0.9329   | 65      | ✓      | ✓      |
| 2.3  | 0.9132   | 84      | ✓      | ✓      |

### Balayage de σ sur mceliece460896

| σ    | p_succès | erreurs | rec. G | rec. L |
|------|----------|---------|--------|--------|
| 0.5  | 1.0000   | 0       | ✓      | ✓      |
| 1.0  | 1.0000   | 0       | ✓      | ✓      |
| 1.5  | 1.0000   | 0       | ✓      | ✓      |
| 2.0  | 1.0000   | 0       | ✓      | ✓      |
| 2.5  | 0.9979   | 3       | ✓      | ✓      |
| 3.0  | 0.9827   | 25      | ✓      | ✓      |

### Balayage de σ sur mceliece8192128

| σ    | p_succès | erreurs | rec. G | rec. L |
|------|----------|---------|--------|--------|
| 0.5  | 1.0000   | 0       | ✓      | ✓      |
| 1.0  | 1.0000   | 0       | ✓      | ✓      |
| 1.5  | 1.0000   | 0       | ✓      | ✓      |
| 2.0  | 1.0000   | 0       | ✓      | ✓      |
| 2.5  | 0.9979   | 3       | ✓      | ✓      |
| 3.0  | 0.9827   | 25      | ✓      | ✓      |
## Structure du dépôt
```text
.
├── 📁 src/
│   ├── Classe_Classic_McEliece.py    # Implémentation du cryptosystème
│   ├── Classe_Decodeur_FFT.py        # Décodeur MLE par FFT
│   └── GoppaReconstructor.py         # Reconstruction algébrique
│
├── 📁 tests/
│   └── tests.py                      # Pipeline de tests complet
|
├── 📁 theoremes/
│   └── theoreme1_identifiabilite.py  # Theoreme 1
│
├── 📁 docs/                          # Mémoire de stage (à venir)
├── 📁 results/                       # Résultats expérimentaux
├── 📄 requirements.txt
├── 📄 .gitignore
└── 📄 README.md
```
## Installation et utilisation

### Prérequis

- Python ≥ 3.9
- Bibliothèques : `numpy`, `galois`, `joblib`

### Installation

```bash
git clone https://github.com/AdyaSarr/Stage-de-Fin-d-annee-M2-Side-Channel-Attack.git
cd Stage-de-Fin-d-annee-M2-Side-Channel-Attack
pip install -r requirements.txt
```

### Exécution

Pour lancer le pipeline complet de test :

```bash
python tests/tests.py mceliece348864       
python tests/tests.py mceliece460896 
```

### Exemple d'utilisation

```python
import sys
sys.path.insert(0, 'src')

from Classe_Classic_McEliece import ClassicMcEliece
from Classe_Decodeur_FFT import DecodorFFT
from GoppaReconstructor import GoppaAttack
import numpy as np

# 1. Génération de clé
mc = ClassicMcEliece('mceliece348864')
mc.key_generation()

# 2. Acquisitions bruitées
decoder = DecodorFFT(mc)
acq = decoder.noise_acquisitions(sigma=2.0)

# 3. Décodage FFT
alphas, betas, scores = decoder.found_all_alpha_beta(
    acq, mode='alpha_any', Nb_couple=968
)

# 4. Reconstruction de la clé
attacker = GoppaAttack(mc)
positions = np.arange(968, dtype=np.int64)
G_hat, L_hat = attacker.full_attack(
    alphas, betas, positions,
    e_max=70, r_for_G=300
)

print(f"G_hat == G ? {G_hat == mc.G}")
print(f"L_hat == L ? {np.array_equal(L_hat, [int(a) for a in mc.L])}")
```

## Cadre théorique

L'attaque repose sur les outils suivants :

- **Décodage par maximum de vraisemblance** (MLE) avec accélération FFT via décimation circulaire et changement de variable par logarithme discret.
- **Algorithme d'interpolation Reed-Solomon avec erreurs** de Daniel J. Bernstein (*Understanding binary-Goppa decoding*, IACR ePrint 2024), avec calcul d'approximants par résolution d'un système linéaire.
- **Reconstruction du support de Goppa** via la matrice de parité publique mise en forme systématique, en exploitant l'injectivité de l'application $\Phi_G$.

Une étude détaillée du seuil de bruit toléré ($\sigma_c$), du nombre minimum de couples requis $n_{min}$(p), ainsi qu'une analyse du cas pathologique des bases normales engendrées par $X$ (premiers d'Artin) sont présentées dans le mémoire.

## Comparaison

Ce modele d'attque presente des avantages comparés aux autres decrit sur les articles en refferences.
Voici un tableau comparatif:

## Tableau comparatif synthétique

| Critère                      | Guo et al. 2022           | Dragoi et al. 2025         | Vallet et al. 2025        | Mon Sujet                                          |
| ---------------------------- | ------------------------- | -------------------------- | ------------------------- | ---------------------------------------------------- |
| **Cible d'implémentation**   | Additive FFT (FPGA, M4)   | Syndrome (matrix-vector)   | Syndrome bruité           | L'algorithme de **decapsulation**                         |
| **Modèle de fuite**          | Power analysis            | Hamming weight             | HW bruité (Gaussien)      | HW bruité (Gaussien, Laplacien)                                 |
| **Méthode du décodeur**      | ML classifier (templates) | Distinguisher déterministe | ML + correction           | MLE par FFT                                          |
| **Bruit toléré (σ)**         | Bas (templates parfaits)  | σ = 0 (modèle idéalisé)    | σ ≤ ~0.38<===>accuracy=0.85 (mceliece348864) | σ ≤ ~2.0<==>accuracy=0.2 (mceliece348864), σ ≤ ~3.0<=====>accuracy=0.1 (mceliece460896) |
| **Reconstruction de G**      | Indirect via support      | Berlekamp-Massey + LFSR    | LFSR avec erreurs         | Bernstein RS interpolation                           |
| **Reconstruction de L**      | Itératif                  | Pivot de Gauss             | Pivot de Gauss            | Pivot de Gauss (identique aux autres)                |
| **Complexité totale**        | O(n)                      | O(n) traces                | O(n³)                     | O(n³)                                                |
|                              | O(n³) algorithme          | O(n³)                      | O(n²m² + (mt)³ + n(mt)²)  | O(n²m² + (mt)³ + n(mt)²)                             |
| **Validation expérimentale** | ChipWhisperer (réel)      | Cortex-M4 (réel)           | Simulation + réel         | Simulation seule (à ce stade)                        |
| **Code source public**       | Non                       | Oui (TCHES 2025)           | Oui (CIC 2025)            | Oui                                |

## Verification des theoremes et propositions
### Théorème 1: Identifiabilite
Ce théorème donne des **conditions necessaires et suffisantes(C1, C2, C1)** pour lesquelles la **séquence de poids** determine de manière **unique** le couple **(a,b)**.

Soient $(a, b), (a', b')$ $\in$ $\mathbb{Z}/P\mathbb{Z}$ $\times$ $\mathbb{Z}/P\mathbb{Z}$ deux couples vérifiant les trois hypothèses suivantes :

-**(C1)** $a$ et $a'$ sont **primitifs** modulo $P$, c'est-à-dire $\gcd(a, P) = \gcd(a', P) = 1$.

-**(C2)** $N \geq N_0(m, \pi)$, où $N_0(m, \pi)$ est le seuil de la séquence.

-**(C3 étendue)** $Stab_{Aff}(f) = \{(1, 0)\}$ : le seul élément du groupe affine préservant $f$ est l'identité.

|   **Parametres**                |   **C1**                                   |**C2**                                        |   **C3**  |
|---------------------------------|--------------------------------------------|----------------------------------------------|-----------|
|   mceliece348864                |  Regarder les 42.19% des elements du corps |   $N_0$ = 18 suffit pour distinguer 42.19% des elements de maniere unique |  Le **stablisateur affine** est **trivial**|
|   mceliece460896                |  Regarder les 99.99% des elements du corps |   $N_0$ = 21 suffit pour distinguer 99.99% des elements de maniere unique |  Le **stablisateur affine** est **trivial**|
|   mceliece6688128               |  Regarder les 99.99% des elements du corps |   $N_0$ = 21 suffit pour distinguer 99.99% des elements de maniere unique |  Le **stablisateur affine** est **trivial**|
|   mceliece6960119               |  Regarder les 99.99% des elements du corps |   $N_0$ = 21 suffit pour distinguer 99.99% des elements de maniere unique |  Le **stablisateur affine** est **trivial**|
|   mceliece8192128               |  Regarder les 99.99% des elements du corps |   $N_0$ = 21 suffit pour distinguer 99.99% des elements de maniere unique |  Le **stablisateur affine** est **trivial**|
|   mceliece348864_poly_cyclo     |  Regarder les 42.19% des elements du corps |   Même a $N_0$ = 256 ne suffit pas pour distinguer les elements |  Le **stablisateur affine** n'est pas **trivial** car $X$ engendre une **base normale** avec le polynome cyclotomique|

### Proposition 3: Critere d'un element soit un element engendre une base normale

Cette proposition permet de caractériser les éléments qui engendrent une base normale du corps $\mathbb{F}_{2^m}$.

Soit $\pi(X) \in \mathbb{F}*2[X]$ un polynôme irréductible de degré $m$, et soit
[
\mathbb{F}*{2^m} \simeq \mathbb{F}_2[X]/(\pi(X)).
]

Notons $M_X \in \mathbb{F}*2^{m \times m}$ la matrice dont la $k$-ième colonne est la représentation binaire de $X^{2^k}$ dans la base canonique
[
(1, X, X^2, \ldots, X^{m-1})
]
de $\mathbb{F}*{2^m}$, c'est-à-dire
[
M_X[i,k]
========

\text{coefficient de } X^i \text{ dans } X^{2^k}
\in \mathbb{F}_2.
]

Alors $X$ engendre une base normale de $\mathbb{F}_{2^m}$ sur $\mathbb{F}*2$ si et seulement si
[
\det*{\mathbb{F}_2}(M_X) \neq 0.
]


## Références

1. **Classic McEliece team.** *Classic McEliece: conservative code-based cryptography*. NIST PQC Round 4 submission (2022). [https://classic.mceliece.org/](https://classic.mceliece.org/)

2. **Daniel J. Bernstein.** *Understanding binary-Goppa decoding*. IACR Communications in Cryptology, 2024. [https://cic.iacr.org/p/1/1/14/](https://cic.iacr.org/p/1/1/14/)

3. **Nicolas Vallet, Brice Colombie, Vlad-Florin Drăgoi, Vincent Grosso, Pierre-Louis Cayrel**[https://cic.iacr.org/p/2/2/26/pdf](https://cic.iacr.org/p/2/2/26/pdf).

4. **Vlad Dragoi, Brice Colombier, Nicolas Vallet, Pierre-Louis Cayrel, Vincent Grosso.**[https://hal.science/hal-04835914/document](https://hal.science/hal-04835914/document)

5. **Michaël Bulois , Pierre-Louis Cayrel , Vlad-Florin Drăgoi , and Vincent Grosso**.[https://hal.science/hal-05621977v1/document](https://hal.science/hal-05621977v1/document)

6. **Annelie Heuser⋆, Olivier Rioul, and Sylvain Guilley**.[https://eprint.iacr.org/2014/527.pdf](https://eprint.iacr.org/2014/527.pdf)

7. **D. PEl, C. WANG AND J.OMURA.**[https://ieeexplore.ieee.org/document/1057152](https://ieeexplore.ieee.org/document/1057152)

8. **Ian F. Blake, XuHong Gao, Ronald C. Mullin, Scott A. Vanstone & Tomik Yaghoobian.** [https://link.springer.com/chapter/10.1007/978-1-4757-2226-0_4](https://link.springer.com/chapter/10.1007/978-1-4757-2226-0_4)

9. **Deking.**[https://perso.eleves.ens-rennes.fr/people/david.michel/agreg/Dedekind.pdf](https://perso.eleves.ens-rennes.fr/people/david.michel/agreg/Dedekind.pdf)
## Auteur

**Adya SARR**  
M2 Mathématiques, Informatique et Applications à la Cryptologie (MIC)  
Université Paris Cité

Stagiaire au Laboratoire Hubert Curien — Université Jean-Monnet, Saint-Étienne

## Licence

Ce code est distribué sous licence MIT. Voir le fichier [LICENSE](LICENSE) pour plus de détails.
