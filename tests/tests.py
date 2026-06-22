"""COMPLETE TESTS SIDE-CHANEL ATTACK OF CLASSIC MCELIECE.
Autor: Adya SARR

Steps Tests :
    [1] Generation of Keys of the Classic McEliece (G, L) + clé publique T
    [2] Noise Acquisitions : y_{i,j} = wt(alpha_i^j · G(alpha_i)^{-2}) + N(0, sigma^2)
    [3] Decoding FFT : (alpha_i, beta_i) = argmax S(alpha, beta)
    [4] Algebric Reconstruction : (alphas, betas) → G_hat → L_hat
    [5] Verification : (G_hat, L_hat) == (G, L) ?

Tests performed :
    - Test 1 : determinist pipeline (sigma=0).
    - Test 2 : noise pipeline (sigma given).
    - Test 3 : scanning sigma.

Auteur : Adya SARR
"""

import numpy as np
import time
import sys

from Classe_Classic_McEliece import ClassicMcEliece
from Classe_Decodeur_FFT import DecodorFFT
from GoppaReconstructor import GoppaAttack


# =========================================================================
#                   Prinyter functions
# =========================================================================

def banner(title, width=70):
    """Print the title"""
    print()
    print("=" * width)
    print(f"  {title}")
    print("=" * width)


def subsection(title, width=70):
    """Print the sub-title"""
    print()
    print("-" * width)
    print(f"  {title}")
    print("-" * width)


def compare_alpha_beta(decoded_alphas, decoded_betas, mc, n_to_check=None):
    """Compare the (alpha_i, beta_i) decoded into tje real values.
    
    Args:
        decoded_alphas, decoded_betas : output of found_all_alpha_beta
        mc : instance of ClassicMcEliece after key_generation
        n_to_check : number of positions to compare
    
    Returns:
        (n_correct, n_total, error_positions) : statistiques decoding
    """
    L = mc.L
    G = mc.G
    
    if n_to_check is None:
        n_to_check = len(decoded_alphas)
    
    true_alphas = np.array([int(a) for a in L[:n_to_check]], dtype=np.int64)
    true_betas = np.array([int(G(L[i])**(-2)) for i in range(n_to_check)],
                          dtype=np.int64)
    
    correct_alpha = (decoded_alphas[:n_to_check] == true_alphas)
    correct_beta = (decoded_betas[:n_to_check] == true_betas)
    correct = correct_alpha & correct_beta
    
    n_correct = int(np.sum(correct))
    error_positions = np.where(~correct)[0]
    
    return n_correct, n_to_check, error_positions


def compute_attack_params(n_pairs, n_errors, t, margin=20, extra_pairs=50):
    """Compute consistent (e_max, r_for_G) parameters for full_attack.
    
    The constraint is : r_for_G ≥ t + 2*e_max + 1.
    
    Strategy :
        1. e_max = n_errors + margin (safety margin against under-estimation)
        2. r_for_G = t + 2*e_max + extra_pairs (extra margin)
        3. Clamp r_for_G ≤ n_pairs
        4. Re-derive e_max from final r_for_G if r_for_G was clamped
    
    Args:
        n_pairs (int) : total available pairs
        n_errors (int) : observed/estimated number of errors
        t (int) : Goppa polynomial degree
        margin (int) : safety margin on e_max above n_errors
        extra_pairs (int) : extra pairs above the minimum r = t + 2*e_max + 1
    
    Returns:
        (e_max, r_for_G) : consistent parameters guaranteeing r ≥ t + 2*e_max + 1.
    """
    # Step 1 : target e_max
    e_max = n_errors + margin
    
    # Step 2 : target r_for_G
    r_for_G = t + 2 * e_max + extra_pairs
    
    # Step 3 : clamp to available pairs
    if r_for_G > n_pairs:
        r_for_G = n_pairs
        # Step 4 : recompute e_max from r_for_G
        e_max = (r_for_G - t - 1) // 2     # ensures r ≥ t + 2*e_max + 1
    
    # Final safety check
    assert r_for_G >= t + 2 * e_max + 1, \
        f"Inconsistent params : r={r_for_G}, t={t}, e_max={e_max}"
    
    return e_max, r_for_G


# =========================================================================
# TEST 1 : Deterministic pipeline (sigma = 0)
# =========================================================================

def test_pipeline_deterministe(param_name='mceliece348864', n_pairs=None):
    """Test of pipeline complet en l'absence de bruit (σ = 0).
    
    Sans bruit, le décodeur FFT doit retrouver TOUS les (α_i, β_i) exacts,
    puis la reconstruction algébrique doit donner (G_hat, L_hat) = (G, L).
    
    Ce test isole les bugs de la reconstruction algébrique (G, L) sans
    interférence du bruit.
    """
    banner(f"TEST 1 : Pipeline déterministe (σ=0) sur '{param_name}'")
    
    # --- [1] Génération de clé ---
    subsection("[1] Génération de clé")
    mc = ClassicMcEliece(param_name)
    print(f"  Paramètres : m={mc.m}, t={mc.t}, n={mc.n}, mt={mc.mt}")
    
    t0 = time.time()
    public_key, secret_key = mc.key_generation()
    print(f"  key_generation : {time.time()-t0:.2f}s")
    
    G, L = secret_key
    n, mt, t = mc.n, mc.mt, mc.t
    
    if n_pairs is None:
        n_pairs = n                          # tout le support
    n_pairs = min(n_pairs, n)
    print(f"  Nombre de couples à acquérir : {n_pairs}")
    
    # --- [2] Acquisitions sans bruit ---
    subsection("[2] Acquisitions (σ=0)")
    decoder = DecodorFFT(mc)
    
    t0 = time.time()
    acq = decoder.noise_acquisitions(sigma=0.0)
    print(f"  noise_acquisitions : shape={acq.shape}, "
          f"{time.time()-t0:.2f}s")
    
    # --- [3] Décodage FFT (mode 'alpha_any') ---
    subsection("[3] Décodage FFT")
    t0 = time.time()
    alphas, betas, scores = decoder.found_all_alpha_beta(
        acq, mode='alpha_any', Nb_couple=n_pairs
    )
    elapsed_fft = time.time() - t0
    print(f"  found_all_alpha_beta : {elapsed_fft:.2f}s "
          f"({elapsed_fft/n_pairs*1000:.1f} ms/couple)")
    
    n_correct, n_total, err_pos = compare_alpha_beta(alphas, betas, mc, n_pairs)
    p_success = n_correct / n_total
    n_errors = n_total - n_correct
    print(f"  Décodage : {n_correct}/{n_total} corrects "
          f"(p_success = {p_success:.4f}, {n_errors} erreurs)")
    if len(err_pos) > 0:
        print(f"  Erreurs aux positions : {err_pos[:10]}{'...' if len(err_pos)>10 else ''}")
    
    # --- [4] Reconstruction algébrique ---
    subsection("[4] Reconstruction algébrique")
    attacker = GoppaAttack(mc)
    positions = np.arange(n_pairs, dtype=np.int64)
    
    # Auto-calcul (e_max, r_for_G) cohérents
    e_max, r_for_G = compute_attack_params(n_pairs, n_errors, t)
    print(f"  r_for_G = {r_for_G}, e_max = {e_max} (n_errors observés = {n_errors})")
    
    t0 = time.time()
    G_hat, L_hat = attacker.full_attack(
        alphas, betas, positions,
        e_max=e_max,
        r_for_G=r_for_G
    )
    elapsed_reco = time.time() - t0
    print(f"  full_attack : {elapsed_reco:.2f}s")
    
    # --- [5] Vérification ---
    subsection("[5] Vérification")
    
    if G_hat is None:
        print("  ÉCHEC : reconstruct_G a échoué")
        return False
    ok_G = (G_hat == G)
    print(f"  G_hat == G ? {ok_G}")
    if not ok_G:
        print(f"    G     = {G}")
        print(f"    G_hat = {G_hat}")
    
    if L_hat is None:
        print("  ÉCHEC : reconstruct_L a échoué")
        return False
    L_int = np.array([int(a) for a in L], dtype=np.int64)
    ok_L = np.array_equal(L_hat, L_int)
    print(f"  L_hat == L ? {ok_L}")
    if not ok_L:
        n_diff = int(np.sum(L_hat != L_int))
        print(f"    Différences : {n_diff} positions")
    
    success = ok_G and ok_L
    subsection("RÉSULTAT")
    if success:
        print(f"  ✓ ATTAQUE COMPLÈTE RÉUSSIE")
        print(f"    Temps total : FFT {elapsed_fft:.1f}s + "
              f"reconstruction {elapsed_reco:.1f}s")
    else:
        print(f"  ✗ ATTAQUE ÉCHOUÉE")
    return success


# =========================================================================
# TEST 2 : Pipeline complet avec bruit
# =========================================================================

def test_pipeline_bruit(param_name='mceliece348864', sigma=1.0, n_pairs=None):
    """Test du pipeline complet avec bruit gaussien d'écart-type σ.
    
    Le bruit introduit des erreurs dans le décodage FFT, qui sont
    ensuite corrigées par l'algorithme de Bernstein dans reconstruct_G.
    
    Args:
        param_name : nom du paramètre Classic McEliece
        sigma : écart-type du bruit gaussien
        n_pairs : nombre de couples à acquérir (default : tout le support)
    """
    banner(f"TEST 2 : Pipeline avec bruit σ={sigma} sur '{param_name}'")
    
    # --- [1] Génération de clé ---
    subsection("[1] Génération de clé")
    mc = ClassicMcEliece(param_name)
    print(f"  Paramètres : m={mc.m}, t={mc.t}, n={mc.n}, mt={mc.mt}")
    
    t0 = time.time()
    public_key, secret_key = mc.key_generation()
    print(f"  key_generation : {time.time()-t0:.2f}s")
    
    G, L = secret_key
    n, mt, t = mc.n, mc.mt, mc.t
    
    # Pour reconstruire L on a besoin de ≥ mt couples corrects.
    # On prend une marge confortable pour absorber les erreurs.
    if n_pairs is None:
        n_pairs = min(mt + 200, n)           # marge de 200 couples
    n_pairs = min(n_pairs, n)
    print(f"  Nombre de couples à acquérir : {n_pairs}")
    
    # --- [2] Acquisitions bruitées ---
    subsection(f"[2] Acquisitions (σ={sigma})")
    decoder = DecodorFFT(mc)
    
    t0 = time.time()
    acq = decoder.noise_acquisitions(sigma=sigma)
    print(f"  noise_acquisitions : shape={acq.shape}, {time.time()-t0:.2f}s")
    
    # --- [3] Décodage FFT ---
    subsection("[3] Décodage FFT (mode 'alpha_any')")
    t0 = time.time()
    alphas, betas, scores = decoder.found_all_alpha_beta(
        acq, mode='alpha_any', Nb_couple=n_pairs
    )
    elapsed_fft = time.time() - t0
    print(f"  found_all_alpha_beta : {elapsed_fft:.2f}s "
          f"({elapsed_fft/n_pairs*1000:.1f} ms/couple)")
    
    n_correct, n_total, err_pos = compare_alpha_beta(alphas, betas, mc, n_pairs)
    p_success = n_correct / n_total
    n_errors = n_total - n_correct
    print(f"  Décodage FFT : {n_correct}/{n_total} corrects "
          f"(p_success = {p_success:.4f}, {n_errors} erreurs)")
    
    # Capacité théorique de correction (avec tous les couples)
    e_max_theorique = (n_pairs - t) // 2
    print(f"  e_max théorique (avec tous les couples) = {e_max_theorique}")
    if n_errors > e_max_theorique:
        print(f"  ⚠ ATTENTION : nombre d'erreurs ({n_errors}) > e_max théorique "
              f"({e_max_theorique})")
        print(f"  L'attaque va probablement échouer.")
    
    # --- [4] Reconstruction algébrique ---
    subsection("[4] Reconstruction algébrique")
    attacker = GoppaAttack(mc)
    positions = np.arange(n_pairs, dtype=np.int64)
    
    # Auto-calcul (e_max, r_for_G) cohérents
    e_max, r_for_G = compute_attack_params(n_pairs, n_errors, t)
    print(f"  r_for_G = {r_for_G}, e_max = {e_max} (n_errors observés = {n_errors})")
    
    t0 = time.time()
    G_hat, L_hat = attacker.full_attack(
        alphas, betas, positions,
        e_max=e_max,
        r_for_G=r_for_G
    )
    elapsed_reco = time.time() - t0
    print(f"  full_attack : {elapsed_reco:.2f}s")
    
    # --- [5] Vérification ---
    subsection("[5] Vérification")
    
    if G_hat is None:
        print("  ÉCHEC : reconstruct_G a échoué (trop d'erreurs ?)")
        return False
    ok_G = (G_hat == G)
    print(f"  G_hat == G ? {ok_G}")
    
    if L_hat is None:
        print("  ÉCHEC : reconstruct_L a échoué")
        return False
    L_int = np.array([int(a) for a in L], dtype=np.int64)
    ok_L = np.array_equal(L_hat, L_int)
    print(f"  L_hat == L ? {ok_L}")
    
    success = ok_G and ok_L
    subsection("RÉSULTAT")
    if success:
        print(f"  ✓ ATTAQUE COMPLÈTE RÉUSSIE (σ={sigma})")
        print(f"    Temps total : FFT {elapsed_fft:.1f}s + "
              f"reconstruction {elapsed_reco:.1f}s")
        print(f"    p_success décodeur : {p_success:.4f}")
    else:
        print(f"  ✗ ATTAQUE ÉCHOUÉE (σ={sigma})")
    return success


# =========================================================================
# TEST 3 : Balayage de sigma
# =========================================================================

def test_sigma_sweep(param_name='mceliece348864', sigma_list=None,
                     n_pairs=None):
    """Test du pipeline pour plusieurs valeurs de σ, pour observer
    le seuil pratique de bruit toléré.
    
    Args:
        param_name : nom du paramètre
        sigma_list : liste des σ à tester
        n_pairs : nombre de couples (constant à travers les σ)
    """
    if sigma_list is None:
        sigma_list = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6]
    
    banner(f"TEST 3 : Balayage σ sur '{param_name}'")
    
    # --- Préparation (génération de clé une seule fois) ---
    mc = ClassicMcEliece(param_name)
    mc.key_generation()
    G, L = mc.G, mc.L
    n, mt, t = mc.n, mc.mt, mc.t
    
    if n_pairs is None:
        n_pairs = min(mt + 200, n)
    n_pairs = min(n_pairs, n)
    
    print(f"  Paramètres : m={mc.m}, t={t}, n={n}, mt={mt}")
    print(f"  Nombre de couples : {n_pairs}, "
          f"e_max théorique max = {(n_pairs-t)//2}")
    
    decoder = DecodorFFT(mc)
    attacker = GoppaAttack(mc)
    positions = np.arange(n_pairs, dtype=np.int64)
    
    # Tableau récapitulatif
    print()
    print(f"  {'σ':>5} | {'p_succès':>10} | {'erreurs':>8} | "
          f"{'r_for_G':>8} | {'e_max':>6} | {'rec. G':>8} | "
          f"{'rec. L':>8} | {'temps':>8}")
    print(f"  {'-'*5}-+-{'-'*10}-+-{'-'*8}-+-{'-'*8}-+-{'-'*6}-+-"
          f"{'-'*8}-+-{'-'*8}-+-{'-'*8}")
    
    results = []
    for sigma in sigma_list:
        acq = decoder.noise_acquisitions(sigma=sigma)
        
        t0 = time.time()
        alphas, betas, _ = decoder.found_all_alpha_beta(
            acq, mode='alpha_any', Nb_couple=n_pairs
        )
        elapsed_fft = time.time() - t0
        
        n_correct, _, _ = compare_alpha_beta(alphas, betas, mc, n_pairs)
        n_errors = n_pairs - n_correct
        p_success = n_correct / n_pairs
        
        # Auto-calcul (e_max, r_for_G) cohérents
        e_max, r_for_G = compute_attack_params(n_pairs, n_errors, t)
        
        t0 = time.time()
        G_hat, L_hat = attacker.full_attack(
            alphas, betas, positions,
            e_max=e_max,
            r_for_G=r_for_G
        )
        elapsed_reco = time.time() - t0
        
        ok_G = (G_hat is not None and G_hat == G)
        if L_hat is not None:
            L_int = np.array([int(a) for a in L], dtype=np.int64)
            ok_L = np.array_equal(L_hat, L_int)
        else:
            ok_L = False
        
        results.append((sigma, p_success, n_errors, r_for_G, e_max,
                        ok_G, ok_L, elapsed_fft + elapsed_reco))
        
        ok_G_str = "✓" if ok_G else "✗"
        ok_L_str = "✓" if ok_L else "✗"
        print(f"  {sigma:>5.2f} | {p_success:>10.4f} | {n_errors:>8d} | "
              f"{r_for_G:>8d} | {e_max:>6d} | "
              f"{ok_G_str:>8s} | {ok_L_str:>8s} | "
              f"{elapsed_fft+elapsed_reco:>7.1f}s")
    
    return results


# =========================================================================
# POINT D'ENTRÉE PRINCIPAL
# =========================================================================

if __name__ == '__main__':
    # Paramètres par défaut. Pour des tests plus grands, passer un argument
    # en ligne de commande : python tests.py 8192128
    
    if len(sys.argv) > 1:
        param = sys.argv[1]
        if param.isdigit():
            param = 'mceliece' + param
    else:
        #param = 'mceliece348864'
        param = 'mceliece460896'
        #param = 'mceliece6688128'
        #param = 'mceliece6960119'
        #param = 'mceliece8192128'
    
    print()
    print("#" * 70)
    print(f"#  PIPELINE COMPLET D'ATTAQUE PAR CANAL AUXILIAIRE")
    print(f"#  contre Classic McEliece — {param}")
    print("#" * 70)
    
    # Sigma utilisé pour le Test 2 (modifiable selon les besoins)
    sigma_test2 = 4.5
    
    # --- Test 1 : déterministe ---
    success_1 = test_pipeline_deterministe(param)
    
    # --- Test 2 : avec bruit ---
    success_2 = test_pipeline_bruit(param, sigma=sigma_test2)
    
    #--- Test 3 : balayage σ (décommenter pour activer) ---
    results = test_sigma_sweep(param,
                               sigma_list=[0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5, 5.5, 6, 6.5])
    
    # --- Résumé final ---
    banner("RÉSUMÉ DES TESTS")
    print(f"  Test 1 (σ=0)          : {'PASSED' if success_1 else 'FAILED'}")
    print(f"  Test 2 (σ={sigma_test2})       : {'PASSED' if success_2 else 'FAILED'}")
    
    if success_1 and success_2:
        print()
        print("  ✓ TOUS LES TESTS SONT PASSÉS")
        sys.exit(0)
    else:
        print()
        print("  ✗ AU MOINS UN TEST A ÉCHOUÉ")
        sys.exit(1)