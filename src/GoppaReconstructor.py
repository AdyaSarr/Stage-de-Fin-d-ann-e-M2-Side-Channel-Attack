"""This class implements the complete attack described in the paper.
It allows to find the complete secret key (G, L) of Classic McEliece
from pairs (alpha_i, beta_i) recovered by side-channel decoding.

Author: Adya SARR
"""
import numpy as np
import galois


class GoppaAttack:
    """Full algebraic reconstruction of the Goppa secret key (G, L)
    from pairs (alpha_i, beta_i) recovered by the FFT side-channel decoder.
    
    Implements:
      - Step 1 (reconstruct_G) : Goppa polynomial reconstruction via
                 Reed-Solomon interpolation with errors (Bernstein 2024).
      - Step 2 (reconstruct_L) : Support reconstruction via Gaussian
                 elimination on the public matrix H_pub.
      - full_attack : pipeline combining the two steps with filtering.
    """
    
    def __init__(self, mc_instance):
        if not hasattr(mc_instance, 'T') or mc_instance.T is None:
            raise RuntimeError(
                "Public key not available. Call mc_instance.key_generation() first."
            )
        
        self.m = mc_instance.m
        self.t = mc_instance.t
        self.n = mc_instance.n
        self.mt = mc_instance.mt
        self.GF2m = mc_instance.GF2m
        self.GF2 = mc_instance.GF2
        
        self.T = mc_instance.T
        I_mt = self.GF2.Identity(self.mt)
        self.H_pub = np.hstack([I_mt, self.T])
        
        self.mc = mc_instance
        self.G_hat = None
        self.L_hat = None
    
    # =====================================================================
    #           STEP 1 : Goppa polynomial reconstruction
    # =====================================================================
    
    def reconstruct_G(self, alphas, betas, e_max=None):
        """Reconstruct the Goppa polynomial G from r ≥ t + 2*e_max pairs with (degP = t and r=n in the paper)
        (alpha_i, beta_tilde_i) where beta_tilde_i = G(alph_i)^{-2} for correct pairs.
        
        Tolerates at least e_max errors among beta_tilde_i.
        
        Args:
            alphas (np.ndarray, dtype int64): recovered alpha_i, shape (r,)
            betas (np.ndarray, dtype int64): recovered beta_tilde_i = G(alpha_i)^{-2}, shape (r,)
            e_max (int or None): maximum number of errors to tolerate.
                If None, defaults to floor((r - t) / 2).
        
        Returns:
            G_hat : galois.Poly over GF2m, of degree t, monic, or None on failure.
        """
        Poly = galois.Poly
        GF = self.GF2m
        t = self.t
        m = self.m
        r = len(alphas)
        
        #verify if I have enough pairs
        if r <= t:
            raise ValueError(f"Not enough pairs: r = {r}, need r > t = {t}.")
        
        if e_max is None:
            e_max = (r - t) // 2
        
        if r < t + 2 * e_max + 1:
            raise ValueError(
                f"With r = {r} pairs, cannot tolerate e_max = {e_max} errors. "
                f"Need r ≥ t + 2*e_max + 1 = {t + 2*e_max + 1}."
            )
        
        # Convert the elements of alphas and betas to elements in F2m
        alpha_gf = GF(alphas.astype(int).tolist())
        beta_gf  = GF(betas.astype(int).tolist())
        
        # === STEP 0 : Transition from values y_tilde===
        exp_sqrt_inv = (2 ** m - 1) - (2 ** (m - 1)) # to have an efficient algorithme
        y_tilde = beta_gf ** exp_sqrt_inv
        
        # === STEP 1 : Reformulation Reed-Solomon ===
        alpha_t = alpha_gf ** t
        r_vec = y_tilde - alpha_t
        
        # === STEP 2a : Lagrange interpolation of B ===
        A = Poly([GF(1)], field=GF)
        for ai in alpha_gf:
            A = A * Poly([GF(1), ai], field=GF)
        
        A_coeffs = A.coeffs
        deg_A = A.degree
        Ap_coeffs = []
        for k in range(deg_A, 0, -1):
            ck = A_coeffs[deg_A - k]
            Ap_coeffs.append(ck * GF(k % 2))
        A_prime = Poly(Ap_coeffs, field=GF)
        #To take the general formule for Lagrange interpolation on page 10 CIC 
        B = Poly([GF(0)], field=GF)
        for i in range(r):
            ai = alpha_gf[i]
            ri = r_vec[i]
            Ap_ai = A_prime(ai)
            if Ap_ai == GF(0):
                return None
            Li = ri / Ap_ai
            A_div = A // Poly([GF(1), ai], field=GF)
            B = B + Poly([Li], field=GF) * A_div
        
        # === STEP 2b : Find an approximation (a, b) of  B/A: Algorithme 3.1.1 on CIC===
        if e_max == 0:
            a = Poly([GF(1)], field=GF)
            b = Poly([GF(0)], field=GF)
        else:
            n_eq = 2 * e_max
            n_var = 2 * e_max + 1
            
            B_low = list(B.coeffs[::-1])
            while len(B_low) < r:
                B_low.append(GF(0))
            
            A_low = list(A.coeffs[::-1])
            while len(A_low) < r + 1:
                A_low.append(GF(0))
            
            M_mat = np.zeros((n_eq, n_var), dtype=int)
            for j in range(n_eq):
                for i in range(e_max + 1):
                    idx = e_max + r - 1 - i - j
                    if 0 <= idx < len(B_low):
                        M_mat[j, i] = int(B_low[idx])
                for i in range(e_max):
                    idx = e_max + r - 1 - i - j
                    if 0 <= idx < len(A_low):
                        M_mat[j, e_max + 1 + i] = int(A_low[idx])
            
            M_gf = GF(M_mat)
            kernel = self._right_kernel(M_gf)
            if kernel is None or len(kernel) == 0:
                return None
            
            ab_vec = kernel[0]
            a_coeffs_low = [ab_vec[i] for i in range(e_max + 1)]
            b_coeffs_low = [ab_vec[e_max + 1 + i] for i in range(e_max)]
            
            a = Poly(a_coeffs_low[::-1], field=GF)
            b = Poly(b_coeffs_low[::-1], field=GF) if any(c != GF(0) for c in b_coeffs_low) \
                else Poly([GF(0)], field=GF)
            
            g_ab = self._poly_gcd(a, b)
            if g_ab.degree > 0:
                a = a // g_ab
                b = b // g_ab
        
        # === STEP 2c : Verify if a divides A ===
        aB_bA = a * B - b * A
        
        if a.degree > A.degree:
            return None
        if A % a != Poly([GF(0)], field=GF):
            return None
        
        bound = r - 2 * e_max + a.degree
        is_zero = (aB_bA == Poly([GF(0)], field=GF))
        if not (is_zero or aB_bA.degree < bound):
            return None
        
        # === STEP 2d : Recover P = B - b·A/a ===
        A_div_a = A // a
        P_hat = B - b * A_div_a
        
        if P_hat.degree >= t:
            return None
        
        # === STEP 3 : Reconstruct G(X) = X^t + P(X) ===
        Xt_coeffs = [GF(0)] * (t + 1)
        Xt_coeffs[0] = GF(1)
        X_to_t = Poly(Xt_coeffs, field=GF)
        
        G_hat = X_to_t + P_hat
        
        if G_hat.degree != t:
            return None
        if G_hat.coeffs[0] != GF(1):
            return None
        
        self.G_hat = G_hat
        return G_hat
    
    # =====================================================================
    #               STEP 2 : Support reconstruction L
    # =====================================================================
    
    def reconstruct_L(self, G, alphas, positions):
        """Reconstruct the full Goppa support L from G, H_pub and CORRECT pairs only.
        
        IMPORTANT : This function ASSUMES all (alphas[k], positions[k]) are correct,
        i.e. that alpha_{positions[k]} = alphas[k] in the true support L.
        
        Returns:
            L_hat (np.ndarray of shape (n,), dtype int64), or None on failure.
        """
        GF = self.GF2m
        n = self.n
        mt = self.mt
        H_pub = self.H_pub
        
        r = len(alphas)
        if r < mt:
            raise ValueError(
                f"Not enough pairs: r = {r}, need r ≥ mt = {mt} to reconstruct L."
            )
        if len(positions) != r:
            raise ValueError(
                f"alphas and positions must have the same length: {r} vs {len(positions)}."
            )
        
        # === STEP 1 : Precompute lookup table T_G ===
        lookup_TG = self._build_lookup_TG(G)
        
        # === STEP 2 : Build the candidate matrix M̃ ===
        alpha_gf = GF(alphas.astype(int).tolist())
        M_tilde = self._build_phi_G_matrix(G, alpha_gf)
        
        # === STEP 3 : Find mt independent columns ===
        J_local = self._find_independent_columns(M_tilde, mt)
        if J_local is None or len(J_local) < mt:
            return None
        
        J_global = positions[J_local]
        
        # === STEP 4 : Compute M = M_J · H_pub[:, J_global]^{-1} ===
        M_J = M_tilde[:, J_local]
        H_pub_J = H_pub[:, J_global]
        
        try:
            H_pub_J_inv = np.linalg.inv(H_pub_J)
        except np.linalg.LinAlgError:
            return None
        
        M = M_J @ H_pub_J_inv
        
        # === STEP 5 + STEP 6 : Recover H[:, j] and α_j for each j ===
        L_hat = np.zeros(n, dtype=np.int64)
        H_full = M @ H_pub
        
        for j in range(n):
            col = np.array(H_full[:, j], dtype=int)
            key = col.tobytes()
            if key not in lookup_TG:
                return None
            L_hat[j] = lookup_TG[key]
        
        if len(set(L_hat.tolist())) != n:
            return None
        
        self.L_hat = L_hat
        return L_hat
    
    # =====================================================================
    # PART 3 : Full attack pipeline
    # =====================================================================
    
    def full_attack(self, alphas, betas, positions, e_max=None, r_for_G=None,
                    verbose=True):
        """Run the complete attack: reconstruct G, filter errors, reconstruct L.
        
        Pipeline:
            [1] reconstruct_G on a subset of r_for_G pairs.
            [2] Filter ALL pairs by direct verification : G_hat(alpha_i)^{-2} == beta_i.
            [3] reconstruct_L using only the verified correct pairs.
        
        Args:
            alphas, betas (np.ndarray): all available pairs
            positions (np.ndarray): support positions of each pair
            e_max (int or None): max errors for reconstruct_G
            r_for_G (int or None): how many pairs to use for reconstruct_G.
                If None, defaults to min(2*t + 60, len(alphas)).
            verbose (bool): print intermediate diagnostics
        
        Returns:
            (G_hat, L_hat) or (None, None) on failure.
        """
        GF = self.GF2m
        r_total = len(alphas)
        
        # --- [1] Reconstruct G on a small subset ---
        if r_for_G is None:
            r_for_G = min(2 * self.t + 60, r_total)
        r_for_G = min(r_for_G, r_total)
        
        alphas_G = alphas[:r_for_G]
        betas_G = betas[:r_for_G]
        
        if e_max is None:
            e_max = (r_for_G - self.t) // 2
        # Make sure e_max is feasible for this subset
        e_max = min(e_max, (r_for_G - self.t) // 2)
        
        if verbose:
            print(f"  [full_attack] reconstruct_G : r={r_for_G}, e_max={e_max}")
        
        G_hat = self.reconstruct_G(alphas_G, betas_G, e_max=e_max)
        if G_hat is None:
            if verbose:
                print("  [full_attack] reconstruct_G FAILED")
            return None, None
        
        # --- [2] Filter erroneous pairs by direct verification ---
        # A pair (alpha_i, beta_i) is correct iff G_hat(alpah_i)^{-2} == beta_i.
        # This works for ANY pair (alpha_i, beta_i) regardless of whether it was used
        # to reconstruct G or not.
        alpha_gf_all = GF(alphas.astype(int).tolist())
        beta_gf_all  = GF(betas.astype(int).tolist())
        
        # Compute G_hat(alpha_i)^{-2} for all i
        G_at_alpha = G_hat(alpha_gf_all)
        
        # Detect alpha_i that are roots of G_hat
        zero_mask = np.array(G_at_alpha == GF(0))
        
        G_at_alpha_safe = G_at_alpha.copy()
        if np.any(zero_mask):
            ones = GF.Ones(int(np.sum(zero_mask)))
            G_at_alpha_safe[zero_mask] = ones
        
        expected_beta = G_at_alpha_safe ** (-2)
        is_correct = np.array(expected_beta == beta_gf_all) & (~zero_mask)
        
        n_correct = int(np.sum(is_correct))
        n_errors = r_total - n_correct
        
        if verbose:
            print(f"  [full_attack] after filtering : {n_correct} correct / "
                  f"{n_errors} erroneous (out of {r_total})")
        
        # Keep only correct pairs
        alphas_clean = alphas[is_correct]
        positions_clean = positions[is_correct]
        
        if verbose:
            print(f"  [full_attack] clean pairs for reconstruct_L : "
                  f"{len(alphas_clean)} (need ≥ {self.mt})")
        
        if len(alphas_clean) < self.mt:
            if verbose:
                print(f"  [full_attack] not enough clean pairs for reconstruct_L")
            return G_hat, None
        
        # --- [3] Reconstruct L on clean pairs ---
        L_hat = self.reconstruct_L(G_hat, alphas_clean, positions_clean)
        if L_hat is None and verbose:
            print(f"  [full_attack] reconstruct_L FAILED")
        
        return G_hat, L_hat
    
    # =====================================================================
    #                               Intermediates methods
    # =====================================================================
    
    def _right_kernel(self, M):
        GF = self.GF2m
        rows, cols = M.shape
        
        M_work = M.copy()
        pivot_cols = []
        row = 0
        for col in range(cols):
            pivot_row = None
            for r in range(row, rows):
                if M_work[r, col] != GF(0):
                    pivot_row = r
                    break
            if pivot_row is None:
                continue
            if pivot_row != row:
                M_work[[row, pivot_row]] = M_work[[pivot_row, row]]
            pivot_val = M_work[row, col]
            M_work[row] = M_work[row] / pivot_val
            for r in range(rows):
                if r != row and M_work[r, col] != GF(0):
                    factor = M_work[r, col]
                    M_work[r] = M_work[r] - factor * M_work[row]
            pivot_cols.append(col)
            row += 1
        
        free_cols = [c for c in range(cols) if c not in pivot_cols]
        kernel = []
        for fc in free_cols:
            vec = GF.Zeros(cols)
            vec[fc] = GF(1)
            for i, pc in enumerate(pivot_cols):
                vec[pc] = -M_work[i, fc]
            kernel.append(vec)
        
        return kernel
    
    def _poly_gcd(self, a, b):
        Poly = galois.Poly
        GF = self.GF2m
        zero = Poly([GF(0)], field=GF)
        
        while b != zero:
            a, b = b, a % b
        
        if a.degree >= 0 and a.coeffs[0] != GF(1):
            leading = a.coeffs[0]
            new_coeffs = [c / leading for c in a.coeffs]
            a = Poly(new_coeffs, field=GF)
        return a
    
    def _build_phi_G_matrix(self, G, alpha_gf):
        GF = self.GF2m
        GF2 = self.GF2
        t = self.t
        mt = self.mt
        
        alpha_arr = np.atleast_1d(alpha_gf)
        r = len(alpha_arr)
        
        G_alpha = G(alpha_arr)
        G_alpha_inv = G_alpha ** (-1)
        
        exponents = np.arange(t)
        V_t = alpha_arr[:, None] ** exponents[None, :]
        H_field = V_t * G_alpha_inv[:, None]
        H_field = H_field.T
        H_vec = H_field.vector()# shape (t, r, m)
        M = H_vec.transpose(0, 2, 1).reshape(mt, r)
        M = GF2(M)
        return M
    
    def _find_independent_columns(self, M_tilde, mt):
        GF2 = self.GF2
        rows, cols = M_tilde.shape
        if cols < mt:
            return None
        
        M_work = M_tilde.copy()
        pivot_cols = []
        row = 0
        for col in range(cols):
            pivot_row = None
            for r in range(row, rows):
                if M_work[r, col] != GF2(0):
                    pivot_row = r
                    break
            if pivot_row is None:
                continue
            if pivot_row != row:
                M_work[[row, pivot_row]] = M_work[[pivot_row, row]]
            for r in range(rows):
                if r != row and M_work[r, col] != GF2(0):
                    M_work[r] = M_work[r] - M_work[row]
            pivot_cols.append(col)
            row += 1
            if len(pivot_cols) == mt:
                break
        
        if len(pivot_cols) < mt:
            return None
        return np.array(pivot_cols, dtype=np.int64)
    
    def _build_lookup_TG(self, G):
        GF = self.GF2m
        all_elements = GF.elements
        
        G_vals = G(all_elements)
        mask = (G_vals != GF(0))
        valid_alphas = all_elements[mask]
        
        M = self._build_phi_G_matrix(G, valid_alphas)
        M_int = np.array(M, dtype=int)
        
        #I construct a dictionary 
        lookup = {}
        for k, alpha in enumerate(valid_alphas):
            key = M_int[:, k].tobytes()# Array into a hashable bytes string
            lookup[key] = int(alpha)
        return lookup