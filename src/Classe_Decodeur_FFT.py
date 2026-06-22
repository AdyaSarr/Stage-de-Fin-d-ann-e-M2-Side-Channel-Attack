""""
    This class allows to built a decodor using FFT
    Autor: Adya SARR
"""

import numpy as np
from math import gcd
from joblib import Parallel, delayed

#====================================================================================================
#                                        Hamming Weight
#====================================================================================================
HW = np.array([bin(i).count("1") for i in range(2 ** 13)], dtype=np.int32)

#====================================================================================================
#                                  Parallel Worker: a candidate a => one or two FFTs
#====================================================================================================
def _score_for_a(a, tilde_y, one_N, F_f, F_h, P, use_quadratic):
    """ For a candidat a in [0, P-1](corresponding an elment alpha=g^a):
        z_l = sum u_j such that a*j=l[P]: Decimation;
        And we compute:
            - S_lin(b) = correlation(z_y, f)(b)
            - S_squad(b) = correlation(z_1, h)(b)

    Args:
        a (int): candidat a : alpha = g^a
        tilde_y (np.array): noise acquisition and padding of lenght P
        one_N (np.array): map identity in [0, N-1]
        F_f (np.array): FFT of f
        F_h (np.array): FFT of h
        P (int): period field
        use_quadratic (bool): _description_
    
    Return: (b_max, best_score)
    """
    j_arr = np.arange(P, dtype=np.int64)
    l_arr = (a * j_arr) % P
    # Decimation of tilde_y of a and correlation with f
    z_y = np.zeros(P, dtype=np.float64)
    #For j, we add tilde_y[j] a z_y[l_arr[j]], so a z_y[a*j mod P] and they are collisions 
    # this instruction accumulate the contribution in z_y.
    np.add.at(z_y, l_arr, tilde_y)
    F_zy = np.fft.fft(z_y)
    #S_lin = (S_1(g^a, g^b))_b
    S_lin = np.fft.ifft(np.conj(F_zy) * F_f).real

    if use_quadratic:
        # Decimation of 1_N per a,  and h's correlation
        z_1 = np.zeros(P, dtype=np.float64)
        np.add.at(z_1, l_arr, one_N)
        F_z1 = np.fft.fft(z_1)
        #S_quad = (S_2(g^a, g^b))_b
        S_quad = np.fft.ifft(np.conj(F_z1) * F_h).real
        S = S_lin - 0.5 * S_quad
    else:
        S = S_lin

    b_max = int(np.argmax(S))
    return b_max, float(S[b_max])

# ============================================================
# Worker for one full acquisition row
# ============================================================
def _decode_one_row(acq_row, candidates, one_N, F_f, F_h, g_pows, P, N, use_quadratic):
    """Decode one acquisition row : find (alpha, beta) maximizing the score.
        
    This worker pads the row internally and iterates over all candidate exponents.
    It is called in parallel by found_all_alpha_beta.
        
    Args:
        acq_row (np.ndarray): one acquisition of length N
        candidates (np.ndarray): array of candidate exponents a ∈ Z/PZ
        one_N, F_f, F_h, g_pows, P, N, use_quadratic : decoder constants
        
    Returns:
        (alpha_int, beta_int, best_score) : decoded pair and its score
    """
    # Pad the acquisition to length P
    tilde_y = np.zeros(P, dtype=np.float64)
    tilde_y[:N] = acq_row
        
    # Iterate over candidates (sequential loop INSIDE this worker)
    best_score, best_a, best_b = -np.inf, -1, -1
    for a in candidates:
        b_max, score = _score_for_a(int(a), tilde_y, one_N, F_f, F_h, P, use_quadratic)
        if score > best_score:
            best_score = score
            best_a = int(a)
            best_b = int(b_max)
        
    alpha_int = int(g_pows[best_a % P])
    beta_int = int(g_pows[best_b])
    return alpha_int, beta_int, best_score

class DecodorFFT:
    """FFT-based decoder for recovering pairs (alpha, beta) from noisy
        Hamming-weight acquisitions, as part of a side-channel attack on Classic McEliece.
    
        Usage:
            mc = ClassicMcEliece('mceliece8192128')
            mc.key_generation()
            decoder = DecodorFFT(mc, N=256, sigma=1.0)
            alpha_hat, beta_hat, score = decoder.decode(tilde_y)
    """
    def __init__(self, mc_instance):
        """Constructor of the class.

            Args:
                mc_instance: an instance of ClassicMcEliece, from which all
                            field parameters (m, t, GF2m, etc.) are inherited.
        """
        #Parameters
        self.mc = mc_instance
        self.m = mc_instance.m
        self.t = mc_instance.t
        self.N = 2*mc_instance.t # Acquisition lenght
        self.P = 2**mc_instance.m - 1
        self.irreducible_poly = mc_instance.irreducible_polynomial
        self.GF2m = mc_instance.GF2m
        
        #Generator over the field
        self.g_int = int(mc_instance.GF2m.primitive_element)
        
        #Tables logs and exponential: g_pows[k] = g^k and log_table[x] = k
        self.g_pows = np.empty(self.P, dtype=np.int64)
        self.log_table = np.full(2**self.m, -1, np.int64)#log_table[0] = -1
        
        a = self.GF2m(1)
        g = self.GF2m(self.g_int)
        
        for k in range(self.P):
            v = int(a)
            self.g_pows[k] = v
            self.log_table[v] = k
            a = a*g
        
        #Function P-periods: f(k) = wt(g^k) and h(k) = f(k)^2
        self.f = HW[self.g_pows].astype(np.float64)
        self.h = self.f * self.f
        
        #FFT for each function f and h
        self.F_f = np.fft.fft(self.f)
        self.F_h = np.fft.fft(self.h)
        
        #indicatrice map of lenght P: 1's N followed by zero's P-N
        self.one_N = np.zeros(self.P, dtype=np.float64)
        self.one_N[:self.N] = 1.0
        
        #all no nulls elements
        self.all_nonzero = np.arange(1, 2**self.m, dtype=np.int64)
        
        #results
        self.results = []


    def noise_acquisitions(self, sigma):
        """Generate noisy acquisitions for all n positions of the McEliece support.
        
        For each position i in {0, ..., n-1}, the acquisition is:
            y_{i,j} = wt(aloha_i^j · G(alpha_i)^{-2}) + N(0, σ²)   for j = 0, ..., N-1
        
        where (G, L) are the secret key of self.mc.
        
        Args:
            sigma (float): noise standard deviation
        
         Returns:
            np.ndarray of shape (n, N), float64 :
                acquisitions[i, j] = wt(alpha_i^j · G(alpha_i)^{-2}) + noise
                - rows : positions of the support (i = 0, ..., n-1)
                - cols : time samples (j = 0, ..., N-1)
                
                Note: this is the TRANSPOSE of the HW_matrix returned by
                ClassicMcEliece.decapsulation_modified with shape (2t, n).
        """
        # verify if the keys exist
        if not hasattr(self.mc, 'G') or self.mc.G is None:
            raise RuntimeError("McEliece key not generated. Call mc.key_generation() first.")
        
        # cumpute H_wt directly on G and L
        # H_wt[j, i] = wt(α_i^j · G(α_i)^{-2})
        G = self.mc.G
        L = self.mc.L
        
        gL_inv_2 = G(L) ** (-2)                            # shape (n,)
        exponents = np.arange(self.N)                       # 0, 1, ..., N-1
        V_N = L[:, None] ** exponents[None, :]              # shape (n, N)
        H_field = (V_N * gL_inv_2[:, None])                 # shape (n, N), in F_{2^m}
        # For each element compute the hamming weigth
        HW_theoretical = np.array(H_field.vector(), dtype=int).sum(axis=2)  # shape (n, N)
        # Add the noise Gaussian
        noise = np.random.normal(0, sigma, size=HW_theoretical.shape)
        return HW_theoretical.astype(np.float64) + noise

    def _pad(self, acq_np):
        """Pad an acquisition of length N to length P with zeros.

        Args:
            acq_np (np.ndarray): acquisition of length N

        Returns:
            np.ndarray of shape (P,) and dtype float64 : zero-padded acquisition
        """
        acq_np = np.asarray(acq_np, dtype=np.float64)
        if acq_np.shape != (self.N,):
            raise ValueError(f"acquisition shape: {acq_np.shape}, expected ({self.N},)")
        tilde_y = np.zeros(self.P, dtype=np.float64)
        tilde_y[:self.N] = acq_np
        return tilde_y
    
    def _candidates_a(self, mode):
        """list of candidates a in Z/PZ

        Args:
            mode (str): gives which mode we can use depending of a
        Return:
            np.array
        """
    
        if mode == "alpha_primitive_fixed":
            #alpha = g
            return np.array([1], dtype=np.int64)
        if mode == "alpha_primitive_any":
            #Take elements that are prime with P
            return np.array(
                [a for a in range(1, self.P) if gcd(a, self.P)==1],
                dtype=np.int64
            )
        if mode == "alpha_any":
            #Take all elements
            return np.arange(0, self.P, dtype=np.int64)
        raise ValueError(f"Unkwon mode:{mode}")

    

    # ============================================================
    # Méthode found_all_alpha_beta (corrigée et optimisée)
    # ============================================================
    def found_all_alpha_beta(self, acq_matrix, mode="alpha_primitive_fixed",
                            Nb_couple=None, n_jobs=-1, use_quadratic=True):
        """For each row of acq_matrix (truncated to Nb_couple first rows), find
        the pair (alpha, beta) ∈ F_{2^m}^* × F_{2^m}^* that maximizes the FFT-based score:
            S(alpha, beta) = <y, w_t(alpha, beta)> - 0.5 ||w_t(alpha, beta)||²
        using parallel FFT.
        
        Args:
            acq_matrix (np.ndarray): noisy acquisitions of shape (n, N), where
                each row i is the acquisition for position alpha_i of the support.
            mode (str): one of:
                - 'alpha_primitive_fixed' : alpha = g (fixed primitive element)
                - 'alpha_primitive_any'   : alpha is any primitive element (a coprime with P)
                - 'alpha_any'             : alpha is any element of F_{2^m}^* (no constraint)
            Nb_couple (int or None, default None): number of rows to decode. If None, decode all rows of acq_matrix.
            n_jobs (int, default -1): number of parallel jobs (-1 = all cores).
            use_quadratic (bool, default True): whether to use the quadratic
                correction term -0.5 ||w_t||² in the score.
        
        Returns:
            alphas (np.ndarray of shape (k,), dtype int64) where k = min(Nb_couple, rows):
                recovered alpha_i (integer representation) for i = 0, ..., k-1.
            betas (np.ndarray of shape (k,), dtype int64):
                recovered beta_i (integer representation).
            scores (np.ndarray of shape (k,), dtype float64):
                score of the best (alpha_i, beta_i) for each acquisition.
        
        Side effects:
            self.results : list of tuples (alpha_int, beta_int, score), updated.
        """
        rows = acq_matrix.shape[0]
        if Nb_couple is None:
            Nb_couple = rows
        else:
            Nb_couple = min(Nb_couple, rows)
        
        #initialise again
        self.results = []
        
        # Constants captured for the workers (read-only)
        candidates = self._candidates_a(mode=mode)
        g_pows = self.g_pows
        P = self.P
        N = self.N
        one_N = self.one_N
        F_f, F_h = self.F_f, self.F_h
        
        acq_matrix_trunc = acq_matrix[:Nb_couple, :]
        
        
        if len(candidates) == 1 or n_jobs == 1 or Nb_couple == 1:
            # For this case it is not necessairy to do parrallelism
            results = [
                _decode_one_row(
                    acq_matrix_trunc[i], candidates,
                    one_N, F_f, F_h, g_pows, P, N, use_quadratic
                )
                for i in range(Nb_couple)
            ]
        else:
            # Parallel over rows
            results = Parallel(n_jobs=n_jobs, backend="loky")(
                delayed(_decode_one_row)(
                    acq_matrix_trunc[i], candidates,
                    one_N, F_f, F_h, g_pows, P, N, use_quadratic
                )
                for i in range(Nb_couple)
            )
        
        # Store and return
        self.results = list(results)
        alphas = np.array([r[0] for r in self.results], dtype=np.int64)
        betas  = np.array([r[1] for r in self.results], dtype=np.int64)
        scores = np.array([r[2] for r in self.results], dtype=np.float64)
        return alphas, betas, scores
        