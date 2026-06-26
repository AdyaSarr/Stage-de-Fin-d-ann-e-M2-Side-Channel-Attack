"""Vérification de la Proposition 3 : critère matriciel effectif pour les bases normales."""

import sys
sys.path.insert(0, 'src')

import numpy as np
import galois

from Classe_Classic_McEliece import CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL


def proposition3(param, theta_int=None, verbose=True):
    """Vérifie si theta engendre une base normale de F_{2^m} sur F_2.
    
    Convention des coefficients dans le dictionnaire spec :
        coefficients HIGH-TO-LOW : 
        coeffs[0] = coeff de X^m, coeffs[m] = coeff constant.
    
    Args:
        param (str): paramètre Classic McEliece.
        theta_int (int or None): élément à tester. None → X (=2).
        verbose (bool): affichage.
    
    Returns:
        dict avec 'engendre_base_normale', 'M_theta', 'rank', 'puissances_Frobenius'.
    """
    m = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['m']
    spec_coeffs = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['irreducible_polynomial']
    
    # ATTENTION : tes coefficients sont en convention HIGH-TO-LOW.
    # galois.Poly attend HIGH-TO-LOW par défaut. Donc PAS d'inversion.
    irred_poly = galois.Poly(spec_coeffs, field=galois.GF(2))
    
    # Vérification de l'irréductibilité
    if not irred_poly.is_irreducible():
        raise ValueError(f"Le polynôme {irred_poly} n'est pas irréductible !")
    
    GF2m = galois.GF(2**m, irreducible_poly=irred_poly)
    GF2 = galois.GF(2)
    
    if verbose:
        print(f"[Prop3] Paramètre : {param}")
        print(f"[Prop3] m = {m}")
        print(f"[Prop3] Polynôme irréductible : {irred_poly}")
    
    # Choix de theta
    if theta_int is None:
        theta_int = 2                                # represente X
    theta = GF2m(theta_int)
    
    if verbose:
        print(f"[Prop3] theta = {theta_int}, "
              f"représentation polynomiale : {galois.Poly.Int(theta_int, field=GF2)}")
    
    # Étape 1 : calculer les conjugués de Frobenius theta^{2^k}
    if verbose:
        print(f"[Prop3] Calcul des conjugués de Frobenius theta^(2^k) pour k=0,...,{m-1}")
    
    puissances_Frob = []
    cur = theta
    for k in range(m):
        puissances_Frob.append(cur)
        if verbose and k < 5:
            print(f"  k={k} : theta^(2^{k}) = {int(cur)}")
        cur = cur ** 2
    
    # Étape 2 : construire M_theta
    if verbose:
        print(f"[Prop3] Construction de M_theta de taille ({m}, {m})")
    
    M_theta = np.zeros((m, m), dtype=np.int8)
    for k in range(m):
        q = int(puissances_Frob[k])
        for i in range(m):
            M_theta[i, k] = (q >> i) & 1
    
    if verbose and m <= 13:
        print(f"[Prop3] M_theta =")
        for i in range(m):
            print(f"  {M_theta[i]}")
    
    # Étape 3 : calculer le rang
    rank = int(np.linalg.matrix_rank(GF2(M_theta)))
    engendre = (rank == m)
    
    if verbose:
        print(f"[Prop3] rang(M_theta) = {rank} (sur {m})")
        if engendre:
            print(f"[Prop3] ✓ theta = {theta_int} ENGENDRE une base normale.")
        else:
            print(f"[Prop3] ✗ theta = {theta_int} N'ENGENDRE PAS une base normale "
                  f"(rang {rank} < {m}).")
    
    return {
        'engendre_base_normale': engendre,
        'M_theta': M_theta,
        'rank': rank,
        'puissances_Frobenius': puissances_Frob,
    }


if __name__ == '__main__':
    print(f"{'Paramètre':<28} {'m':>3} {'Polynôme':<40} {'rang':>5} {'Base normale ?'}")
    print("-" * 140)
    for param in CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL:
        m = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['m']
        coeffs = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['irreducible_polynomial']
        irred = galois.Poly(coeffs, field=galois.GF(2))
        
        result = proposition3(param, theta_int=2, verbose=True)
        verdict = "✓ OUI" if result['engendre_base_normale'] else "✗ NON"
        print(f"{param:<28} {m:>3} {str(irred):<40} {result['rank']:>5} {verdict}")