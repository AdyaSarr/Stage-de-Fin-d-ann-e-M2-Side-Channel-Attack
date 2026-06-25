"""""Verification du Théorème 1:
        Soient (a, b), (a', b') dans Z/PZxZ/PZ deux couples vérifiant les trois hypothèses suivantes :
        -(C1): a et a' sont primitifs modulo P, c'est-à-dire PGCD(a, P) = PGCD(a', P) = 1.
        -(C2): N >= N_0(m, pi), où N_0(m, pi) est le seuil.
        -(C3 étendu): Stab_{Aff}(f) = {(1, 0)} : le seul élément du groupe affine préservant f est l'identité.
        Alors:
        Quelque soit  i dans {0, 1,......., N-1}, w^{(a,b)}_i = w^{(a',b')}_i <===========> (a, b) = (a', b').
"""""

import sys
import numpy as np
from math import gcd
import time
sys.path.insert(0, 'src')

from Classe_Classic_McEliece import CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL

#===========================================================================
#                           Arithmetique sur GF(2^m)
#===========================================================================
def gf_mul(a, b, m, reduction):
    result = 0
    while b:
        if b & 1: result ^= a
        b >>= 1
        a <<= 1
        if a & (1 << m): a ^= (1 << m); a ^= reduction
    return result


def gf_pow(base, e, m, reduction):
    result = 1
    while e > 0:
        if e & 1: result = gf_mul(result, base, m, reduction)
        e >>= 1
        base = gf_mul(base, base, m, reduction)
    return result

#===========================================================================
#                           Construction du Corps Fini
#===========================================================================
def build_field(param):
    m = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['m']
    poly_coeffs = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['irreducible_polynomial']
    P = 2**m - 1
    if poly_coeffs is None:
        raise ValueError(f"param={param} n'est par supporté")
    poly_int = 0
    deg = len(poly_coeffs) - 1
    for i, c in enumerate(poly_coeffs):
        if c: poly_int |= (1 << (deg - i))
    reduction = poly_int ^ (1 << m)
    # Detecter si X (=2) est primitif
    a = 1
    pows = [1]
    order_X = None
    for k in range(1, P + 1):
        a <<= 1
        if a & (1 << m): a ^= (1 << m); a ^= reduction
        pows.append(a)
        if a == 1: order_X = k; break

    if order_X == P:
        g = 2
        Tab_exp = np.array(pows[:P], dtype=np.int64)
    else:
        def prime_factors(n):
            factors = set(); d = 2
            while d * d <= n:
                while n % d == 0: factors.add(d); n //= d
                d += 1
            if n > 1: factors.add(n)
            return factors

        primes = prime_factors(P)
        g = None
        for cand in range(2, 1 << m):
            if all(gf_pow(cand, P // p, m, reduction) != 1 for p in primes):
                g = cand; break

        Tab_exp = np.zeros(P, dtype=np.int64)
        x = 1
        for k in range(P):
            Tab_exp[k] = x
            x = gf_mul(x, g, m, reduction)

    Tab_log = np.full(2 ** m, -1, dtype=np.int64)
    for k in range(P):
        Tab_log[Tab_exp[k]] = k

    return Tab_exp, Tab_log, g, reduction

#===========================================================================
#                           Verification C3 étendu
#===========================================================================
def verifier_C3_etendu(param, verbose=True, early_exit=True):
    """
        Stab_Aff(f) = {(c,e) dans Aff(Z/PZ): f(cj+e) = f(j) pour j dans Z/PZ}
        Args:
            param(str): Nom du parametre Classic McEliece
            verbose(bool): affichage
            early_exit(bool): retourne dés qu'un element non-trivial est trouvé
    """
    Tab_exp, Tab_log, g, reduction = build_field(param=param)
    m = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['m']
    P = 2**m - 1
    
    if verbose:
        print(f"[C3 étendu sur] {param}: m={m}, P = {P}")
    #Precalculer tout f(k) = wt(g^k)
    f = np.array([bin(int(x)).count('1') for x in Tab_exp], dtype=np.int8)
    #Prendre les elements premiers avec P
    units = np.array([u for u in range(1, P) if gcd(u, P)==1], dtype=np.int16)
    
    if verbose:
        print(f"[C3 étendu] Card(units) = {len(units)}, candidates = {len(units)*P:,}")
    
    k_range = np.arange(P, dtype=np.int64)
    stab_elements = []
    
    t0 = time.time()
    for ui, u in enumerate(units):
        uk_mod = (u*k_range)%P
        f_uk = f[uk_mod]
        g_u = np.empty_like(f)
        g_u[uk_mod] = f
        for v in range(P):
            f_shifted = np.roll(f, v)
            if np.array_equal(g_u, f_shifted):
                if not (u == 1 and v == 0):
                    stab_elements.append((int(u), int(v)))
                    if early_exit:
                        if verbose:
                            print(f"[C3 etendu] le premier element non-trivial du stablisateur: (u={u}, v={v})")
                        return{
                            'C3_satisfied': False,
                            'stab_elements': stab_elements,
                            'm': m,
                            'P': P
                        }
    elapsed = time.time() - t0
    if verbose:
        print(f"[C3 étendu] Done in {elapsed:.1f}s. "
              f"Non-trivial stab : {len(stab_elements)}")
    return {
        'C3_satisfied': len(stab_elements) == 0,
        'stab_elements': stab_elements,
        'm': m,
        'P': P,
    }
    
#===========================================================================
#                           Verification de l'hypothese C2
#===========================================================================
def verifier_C2(param, N, verbose=True):
    """Vérifie la condition (C2) : N >= N_0(m, π).
    
    Le seuil N_0 dépend du polynôme irréductible π. 
    Pour simplifier, on utilise le critère pratique :
        N_0 ≈ taille minimale pour discriminer toutes les paires (a, b)
        avec a primitif.
    
    Args:
        param (str): nom du paramètre Classic McEliece
        N (int): valeur de N à tester
    
    Returns:
        dict avec 'C2_satisfied', 'N', 'N_0_estimated'.
    """
    Tab_exp, Tab_log, g, reduction = build_field(param)
    m = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['m']
    P = 2**m - 1
    
    # Precalculer f
    f = np.array([bin(int(x)).count('1') for x in Tab_exp], dtype=np.int8)
    #Tout les elements premiers avec P
    units = [u for u in range(1, P) if gcd(u, P) == 1]
    
    # on calcule la sequence (f(j*a + b mod P))_{j=0,...,N-1} pour chaqque couple (a, b)
    # et on regarde la distinguabilité.
    
    if verbose:
        print(f"[C2] Test pour N={N} et pour {param} (m={m}, P={P})")
    
    sigs = {}
    n_collisions = 0
    
    import time
    t0 = time.time()
    
    for a in units:
        # Precalculer (j*a) mod P for j=0,...,N-1
        ja_mod = (np.arange(N, dtype=np.int64) * a) % P
        
        for b in range(P):
            idx = (ja_mod + b) % P
            sig = f[idx].tobytes()
            
            if sig in sigs:
                n_collisions += 1
            else:
                sigs[sig] = (a, b)
    
    elapsed = time.time() - t0
    
    n_total = len(units) * P
    C2_satisfied = (n_collisions == 0)
    
    if verbose:
        print(f"[C2] Testés {n_total:,} (a, b) pairs en {elapsed:.1f}s")
        print(f"[C2] Sequences distinctes: {len(sigs):,}")
        print(f"[C2] Collisions: {n_collisions:,}")
        if C2_satisfied:
            print(f"[C2] ✓ N={N} est suffisant")
        else:
            print(f"[C2] ✗ N={N} est insuffisant ({n_collisions} collisions)")
    
    return {
        'C2_satisfied': C2_satisfied,
        'N': N,
        'n_collisions': n_collisions,
    }

#===========================================================================
#                         Vérification complete du theoreme
#===========================================================================
def verifier_theoreme(param, N=None, verbose=True):
    """Vérifie globalement le théorème d'identifiabilité pour un paramètre.
    
    Stratégie :
        - (C1) : on quantifie combien de couples (a, b) satisfont gcd(a, P) = 1.
                 Cette condition filtre les (a, b) admissibles.
        - (C2) : on teste pour le N pratique (par defaut 2t) si tous les
                 couples (a, b) avec a primitif donnent des signatures distinctes.
        - (C3) : on vérifie que Stab_Aff(f) = {(1, 0)}.
    
    Le théorème est dit "vérifié pour ces paramètres" si (C2) et (C3) le sont.
    
    Args:
        param (str): nom du paramètre Classic McEliece (ex: 'mceliece348864').
        N (int or None): valeur de N à tester pour (C2). Par défaut, on prend N=2t.
        verbose (bool): affichage détaillé.
    
    Returns:
        dict avec les clés :
            'theoreme_verifie' (bool) : True si C2 et C3 satisfaites
            'C1_info' (dict) : statistiques sur (C1)
            'C2_result' (dict) : sortie de verifier_C2
            'C3_result' (dict) : sortie de verifier_C3_etendu
            'condition_violee' (str or None) : nom de la première condition non vérifiée
    """
    # --- Récupération des paramètres ---
    m = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['m']
    t = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[param]['t']
    P = 2**m - 1
    
    if N is None:
        N = 2 * t                                    # valeur par défaut
    
    if verbose:
        print()
        print("#" * 72)
        print(f"#  VÉRIFICATION DU THÉORÈME D'IDENTIFIABILITÉ")
        print(f"#  Paramètre : {param}")
        print(f"#  m = {m}, t = {t}, P = 2^m - 1 = {P}, N = {N}")
        print("#" * 72)
        print()
    
    # ===Pour C1==
    if verbose:
        print("-" * 72)
        print("(C1) : a et a' primitifs modulo P, i.e. gcd(a, P) = 1")
        print("-" * 72)
    
    n_units = sum(1 for u in range(1, P) if gcd(u, P) == 1)
    n_couples_C1 = n_units * P                       
    n_couples_total = P * P
    proportion_C1 = n_couples_C1 / n_couples_total
    
    if verbose:
        print(f"  Nombre d'éléments a primitifs                 : {n_units:,}")
        print(f"  Couples (a, b) satisfaisant (C1)              : {n_couples_C1:,}")
        print(f"  Couples (a, b) totaux dans Z/PZ x Z/PZ        : {n_couples_total:,}")
        print(f"  Proportion satisfaisant (C1)                  : {proportion_C1:.4%}")
        print(f"  Note : (C1) est une hypothèse sur les couples ; elle filtre")
        print(f"            elements du corps qu'on peut distinguer")
        print()
    
    C1_info = {
        'n_units': n_units,
        'n_couples_satisfying_C1': n_couples_C1,
        'n_couples_total': n_couples_total,
        'proportion': proportion_C1,
    }
    
    # === Pour C2 ===
    if verbose:
        print("-" * 72)
        print(f"(C2) : N = {N} >= N_0(m, pi) ?")
        print("-" * 72)
    
    result_C2 = verifier_C2(param, N=N, verbose=verbose)
    
    if verbose:
        print()
    
    # Si C2 violée, on arrête : le théorème n'est pas vérifié pour ce N
    if not result_C2['C2_satisfied']:
        if verbose:
            print("=" * 72)
            print(f"VERDICT : THÉORÈME NON VÉRIFIÉ pour {param} avec N = {N}")
            print(f"  Condition violée : (C2) — N = {N} insuffisant pour discriminer")
            print(f"                     {result_C2['n_collisions']} couples produisent")
            print(f"                     de séquences identiques.")
            print(f"  Augmenter N et retester.")
            print("=" * 72)
        return {
            'theoreme_verifie': False,
            'C1_info': C1_info,
            'C2_result': result_C2,
            'C3_result': None,
            'condition_violee': 'C2',
        }
    
    # === Pour C3 étendu ===
    if verbose:
        print("-" * 72)
        print("(C3 étendue) : Stab_Aff(f) = {(1, 0)} ?")
        print("-" * 72)
    
    result_C3 = verifier_C3_etendu(param, verbose=verbose, early_exit=False)
    
    if verbose:
        print()
    
    # === Verdict final ===
    if verbose:
        print("=" * 72)
        print("VERDICT FINAL")
        print("=" * 72)
    
    if result_C3['C3_satisfied']:
        if verbose:
            print(f"  THÉORÈME VÉRIFIÉ pour {param} avec N = {N}")
            print(f"    (C1) : {n_couples_C1:,} couples admissibles "
                  f"({proportion_C1:.2%} du total)")
            print(f"    (C2) : N = {N} suffisant (aucune collision)")
            print(f"    (C3) : Stab_Aff(f) = {{(1, 0)}} (trivial)")
            print(f"  Conclusion : pour tout couple (a, b) avec a primitif,")
            print(f"               la séquenece (w_i^(a,b))_{{i=0..{N-1}}} détermine")
            print(f"               de manière unique le couple (a, b).")
            print("=" * 72)
        return {
            'theoreme_verifie': True,
            'C1_info': C1_info,
            'C2_result': result_C2,
            'C3_result': result_C3,
            'condition_violee': None,
        }
    else:
        if verbose:
            print(f"  THÉORÈME NON VÉRIFIÉ pour {param}")
            print(f"  Condition violée : (C3 étendue) — Stab_Aff(f) contient")
            print(f"                     {len(result_C3['stab_elements'])} éléments")
            print(f"                     non triviaux.")
            print(f"  Premiers éléments du stabilisateur :")
            for (u, v) in result_C3['stab_elements'][:5]:
                print(f"      (u = {u}, v = {v})")
            print(f"  Conséquence : il existe des couples (a, b) ≠ (a', b')")
            print(f"                produisant la même séquence, donc le décodeur")
            print(f"                FFT a une ambiguïté inhérente.")
            print("=" * 72)
        return {
            'theoreme_verifie': False,
            'C1_info': C1_info,
            'C2_result': result_C2,
            'C3_result': result_C3,
            'condition_violee': 'C3',
        }


#===========================================================================
#                           Point d'entree
#===========================================================================
if __name__ == '__main__':
    if len(sys.argv) > 1:
        param = sys.argv[1]
        if param.isdigit():
            param = 'mceliece' + param
    else:
        param = 'mceliece348864'
    N = 18
    result = verifier_theoreme(param, N, verbose=True)
    
    # Code de sortie pour usage en script
    sys.exit(0 if result['theoreme_verifie'] else 1)