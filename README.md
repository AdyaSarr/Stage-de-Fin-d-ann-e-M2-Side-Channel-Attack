# Stage M2 — Attaque par Canal Auxiliaire contre Classic McEliece

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Ce dépôt contient l'implémentation complète d'une attaque par canal auxiliaire (*side-channel attack*) contre le cryptosystème post-quantique **Classic McEliece**. Le travail est mené dans le cadre d'un stage de fin d'études M2 au **Laboratoire Hubert Curien** (Saint-Étienne), sous l'encadrement de **Pierre-Louis Cayrel(Maître de conférences HDR UJM - LabHC), Brice Colombier(chercheur et professeur en electronique et programmation a IUT Saint-Etienne), Michaël Bulois(Maitre de Conférences en Mathématiques  au laboratoire Institut Camille Jordan), Vincent GROSSO (CR Researcher, CNRS) et Vlad-Florin Dragoi (Professor (Associate) at Aurel Vlaicu University of Ara in Roumania)**

## Contexte

Classic McEliece est un cryptosystème post-quantique basé sur les **codes de Goppa binaires**, finaliste du standard NIST PQC. Sa sécurité repose sur la difficulté de décoder un code aléatoire et la dissimulation de la structure de Goppa via une permutation secrète.

Ce travail étudie une **attaque par canal auxiliaire** exploitant des fuites lors de la décapsulation. À partir d'acquisitions bruitées de poids de Hamming :

$$y_{i,j} = \mathrm{wt}(\alpha_i^j \cdot G(\alpha_i)^{-2}) + \varepsilon_{i,j}, \qquad \varepsilon_{i,j} \sim \mathcal{N}(0, \sigma^2)$$

l'attaquant reconstruit intégralement la clé secrète $(G, L)$ : le polynôme de Goppa $G$ et le support $L = (\alpha_0, \ldots, \alpha_{n-1})$.

## Pipeline d'attaque

┌─────────────────────────────────────────────────────────────────┐

───────────────────┐

│1. Décodage MLE par FFT│$

│Acquisitions bruitées (n×N) → couples (α_i, β_i = G(α_i)^-2)│

│ Coût : O(P² log P) où P = 2^m - 1│
├─────────────────────────────────────────────────────────────────
─────────┤
│2. Reconstruction algébrique de G│

│Algorithme de Bernstein (2024) — interpolation Reed-Solomon│

│avec correction d'erreurs│

│Coût : O(n² m²) opérations binaires│
├─────────────────────────────────────────────────────────────────
│3. Reconstruction du support L│

│Pivot de Gauss sur la matrice de parité publique│

│Coût : O((mt)³ + n(mt)²) opérations binaires│
└─────────────────────────────────────────────────────────────────
─────────┘
## Résultats expérimentaux

### Validation sur paramètres officiels Classic McEliece

| Paramètre        | m  | t   | n    | mt   | σ testé | Attaque |
|------------------|----|-----|------|------|---------|---------|
| mceliece348864   | 12 | 64  | 3488 | 768  | 0 – 2.3 | ✓       |
| mceliece348864   | 12 | 64  | 3488 | 2.5+ | ✗ (seuil atteint) |
| mceliece460896   | 13 | 96  | 4608 | 1248 | 0 – 3.0+| ✓       |

L'attaque réussit complètement (reconstruction de $G$ **et** de $L$) tant que le taux d'erreur du décodeur FFT reste sous environ 5 % pour `mceliece348864` et sous 3 % pour `mceliece460896`.

### Balayage de σ sur mceliece460896

| σ    | p_succès | erreurs | rec. G | rec. L |
|------|----------|---------|--------|--------|
| 0.5  | 1.0000   | 0       | ✓      | ✓      |
| 1.0  | 1.0000   | 0       | ✓      | ✓      |
| 1.5  | 1.0000   | 0       | ✓      | ✓      |
| 2.0  | 1.0000   | 0       | ✓      | ✓      |
| 2.5  | 0.9979   | 3       | ✓      | ✓      |
| 3.0  | 0.9827   | 25      | ✓      | ✓      |

## Structure du dépôt
.
├──────src/
│   ├── Classe_Classic_McEliece.py    # Implémentation du cryptosystème
│   ├── Classe_Decodeur_FFT.py        # Décodeur MLE par FFT
│   └── GoppaReconstructor.py         # Reconstruction algébrique

├──────tests/
│   └── tests.py                      # Pipeline de tests complet
├──────docs/                             # Mémoire de stage (à venir)
├──────results/                          # Résultats expérimentaux
├──────requirements.txt
├──────.gitignore
└──────README.md
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

### Exemple d'utilisation programmatique

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

Une étude détaillée du seuil de bruit toléré ($\sigma_c$), du nombre minimum de couples requis ($n_{\min}(p)$), ainsi qu'une analyse du cas pathologique des bases normales engendrées par $X$ (premiers d'Artin) sont présentées dans le mémoire.

## Références

1. **Classic McEliece team.** *Classic McEliece: conservative code-based cryptography*. NIST PQC Round 4 submission (2022). [https://classic.mceliece.org/](https://classic.mceliece.org/)

2. **Daniel J. Bernstein.** *Understanding binary-Goppa decoding*. IACR Communications in Cryptology, 2024. [https://cic.iacr.org/p/1/1/14/](https://cic.iacr.org/p/1/1/14/)

3. **Nicolas Vallet, Brice Colombie, Vlad-Florin Drăgoi, Vincent Grosso, Pierre-Louis Cayrel**[https://cic.iacr.org/p/2/2/26/pdf](https://cic.iacr.org/p/2/2/26/pdf).

4. **Vlad Dragoi, Brice Colombier, Nicolas Vallet, Pierre-Louis Cayrel, Vincent Grosso.**[https://hal.science/hal-04835914/document](https://hal.science/hal-04835914/document)

5. **Michaël Bulois , Pierre-Louis Cayrel , Vlad-Florin Drăgoi , and Vincent Grosso**.**<Algebraic Key-Recovery Side-Channel Attack on Classic McEliece>**

6. **Annelie Heuser⋆, Olivier Rioul, and Sylvain Guilley**.[https://eprint.iacr.org/2014/527.pdf](https://eprint.iacr.org/2014/527.pdf)
## Auteur

**Adya SARR**  
M2 Mathématiques, Informatique et Applications à la Cryptologie (MIC)  
Université Paris Cité

Stagiaire au Laboratoire Hubert Curien — Université Jean-Monnet, Saint-Étienne

## Licence

Ce code est distribué sous licence MIT. Voir le fichier [LICENSE](LICENSE) pour plus de détails.
