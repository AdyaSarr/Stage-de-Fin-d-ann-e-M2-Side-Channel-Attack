"""Classic McEliece Cryptosystem: this class implemente the  cryptosystem ClassicMcEliece
    From key genaration to the encapsulation and the decapsulation.
    Autor: Adya SARR
"""

import numpy as np
import galois as gf 
from random import choice

#====================================================================================================
#                               The official irreducible polynomial of Classic McEliece
#====================================================================================================
CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL = {
    'test': {'m':4, 't': 3, 'n': 16, 'irreducible_polynomial': [1, 0, 0, 1, 1]},
    'mceliece348864':  {'m': 12, 't': 64,  'n': 3488, 'irreducible_polynomial':   [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1]},
    'mceliece460896':  {'m': 13, 't': 96,  'n': 4608, 'irreducible_polynomial':   [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]},
    'mceliece6688128': {'m': 13, 't': 128, 'n': 6688, 'irreducible_polynomial':   [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]},
    'mceliece6960119': {'m': 13, 't': 119, 'n': 6960, 'irreducible_polynomial':   [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]},
    'mceliece8192128': {'m': 13, 't': 128, 'n': 8192, 'irreducible_polynomial':   [1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1]},
}


#====================================================================================================
#                               The class Classic McEliece
#====================================================================================================
class ClassicMcEliece:
    """This class implements the cryptosystem Classic McEliece, from key genaration to the encapsulation and the decapsulation.
    The parameters of the cryptosystem are defined in the dictionary CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL.
    The constructor takes as input the name of the cryptosystem, which is one of the keys of the dictionary CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL.
    The class has three main methods: key_generation, encapsulation and decapsulation."""
    
    def __init__(self, name):
        """Constructor of the class ClassicMcEliece. It takes as input the name of the cryptosystem, which is one of the keys of the dictionary 
        CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL."""
        if name not in CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL:
            raise ValueError("Invalid name for Classic McEliece. Valid names are: " + ", ".join(CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL.keys()))
        self.name = name
        self.m = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[name]['m']
        self.t = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[name]['t']
        self.n = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[name]['n']
        self.irreducible_polynomial = CLASSIC_MCELIECE_IRREDUCIBLE_POLYNOMIAL[name]['irreducible_polynomial']
        self.GF2m = gf.GF(2**self.m, irreducible_poly=self.irreducible_polynomial)
        self.GF2 = gf.GF(2)
        self.mt = self.m * self.t
    def key_generation(self):
        """Key generation of the Classic McEliece Cryptosystem.
        It generates the public key and the secret key of the cryptosystem, 
        following the parameters given by the class.
        
        Returns:
            public_key: the public key T of the cryptosystem (matrix of shape (mt, n-mt) in F_2)
            secret_key: tuple (G, L) where G is the Goppa polynomial and L is the support
        
        Side effects:
            self.G     : Goppa polynomial (irreducible, monic, degree t)
            self.L     : Goppa support (n distinct elements of F_{2^m})
            self.H_F2m : parity-check matrix in F_{2^m}, shape (t, n)
            self.H     : binary parity-check matrix, shape (mt, n)
            self.T     : public key, shape (mt, n - mt)
            self.S_inv : matrix S^{-1} used to obtain systematic form, shape (mt, mt)
        """
        mt = self.m * self.t
        GF2 = gf.GF(2)
        
        while True:
            # === Step 1 : Generate the Goppa polynomial G ===
            while True:
                g = gf.Poly.Random(self.t, field=self.GF2m)
                if g.is_irreducible:
                    break
            if not g.is_monic:
                g = g * (g.coeffs[0] ** -1)
            assert g.is_irreducible and g.is_monic and g.degree == self.t
            self.G = g

            # === Step 2 : Generate the support L of the Goppa code ===
            while True:
                L = self.GF2m(np.random.choice(self.GF2m.elements, self.n, replace=False))
                if np.all(self.G(L) != 0):
                    break
        
            self.L = L

            # === Step 3 : Generate the parity-check matrix H_F2m in F_{2^m}^{t x n} ===
            # H_F2m[j, i] = alpha_i^j · G(alpha_i)^{-1}
            gL_inv = self.G(self.L) ** (-1)
            exponents = np.arange(self.t)
            V_t = self.L[:, None] ** exponents[None, :]           # shape (n, t) in F_{2^m}
            H_F2m = (V_t * gL_inv[:, None]).T                     # shape (t, n) in F_{2^m}
            assert isinstance(H_F2m, self.GF2m), "H_F2m should be in F_{2^m}"
            assert H_F2m.shape == (self.t, self.n)
            self.H_F2m = H_F2m

            # === Step 4 : Binary expansion of H_F2m into H in F_2^{mt x n} ===
            H_F2m_vec = H_F2m.vector()                            # shape (t, n, m)
            H_F2 = H_F2m_vec.transpose(0, 2, 1).reshape(mt, self.n)
            H_F2 = GF2(H_F2)
            assert H_F2.shape == (mt, self.n)
            self.H = H_F2

            # === Step 5 : Systematic form (I_mt | T) ===
            M = H_F2[:, :mt]   # left block of H, shape (mt, mt)
            try:
                M_inv = np.linalg.inv(M)
            except np.linalg.LinAlgError:
                # M is singular: the systematic form is not possible.
                # Following the official Classic McEliece spec, regenerate everything.
                continue
            
            # Apply M_inv on the left to put H in systematic form
            H_sys = M_inv @ H_F2                                  # shape (mt, n)
            
            # Sanity check: the left block must be the identity
            I_mt = GF2.Identity(mt)
            assert np.array_equal(H_sys[:, :mt], I_mt), "Left block is not the identity"
            
            # === Step 6 : Extract the public key T and return ===
            T = H_sys[:, mt:]                                     # shape (mt, n - mt)
            assert T.shape == (mt, self.n - mt)
            
            self.T = T
            self.S_inv = M_inv
            
            public_key = T
            secret_key = (self.G, self.L)
            return public_key, secret_key
    
    def encapsulation(self, public_key):
        """This function allows to compute the encapsulation of an element which hamming weight t

        Args:
            public_key (fieldArray): The second party of the systematic matrix 

        Returns:
            fieldArray: the encapsulation 
        """
        positions = np.random.choice(self.n, self.t, replace=False)
        e = np.zeros(self.n, dtype=int)
        e[positions] = 1
        e = self.GF2(e)
        if int(np.sum(np.array(e, dtype=int))) != self.t:
            raise ValueError(f"Error vector has weight {int(np.sum(np.array(e, dtype=int)))}, expected {self.t}")
        mt = self.m*self.t
        I_mt = self.GF2.Identity(mt)
        
        H_pub = np.hstack((I_mt, public_key))
        assert H_pub.shape == (mt, self.n)
        return e, H_pub @ e
    
    def decapsulation_modified(self, cipher):
        """Decapsulates the ciphertext following the Goppa decoding algorithm
        (Algorithm 1 from the reference paper), and returns additionally the
        Hamming weight matrix of H_{priv, γ^2} for side-channel analysis.
        
        Args:
            cipher: ciphertext z, shape (mt,) over F_2
        
        Returns:
            e: error vector of weight t recovered from the ciphertext, shape (n,)
            HW_matrix: integer matrix of shape (2t, n), where 
                    HW_matrix[j, i] = wt(α_i^j · G(α_i)^{-2})
        
        Raises:
            RuntimeError: if the decoding fails (e.g. > t errors or corruption)
        """
        # === Step 1 : Pad the ciphertext to length n ===
        v = np.concatenate([np.array(cipher, dtype=int),
                            np.zeros(self.n - self.mt, dtype=int)])
        v_F2m = self.GF2m(v)
        
        # === Step 2 : Build H_{priv, γ^2} ∈ F_{2^m}^{2t × n} ===
        # H_priv_g2[j, i] = α_i^j · G(α_i)^{-2}
        gL_inv_2 = self.G(self.L) ** (-2)
        exponents = np.arange(2 * self.t)
        V_2t = self.L[:, None] ** exponents[None, :]
        H_priv_g2 = (V_2t * gL_inv_2[:, None]).T
        
        # === Step 3 : Compute the syndrome s = H_{priv, γ^2} @ v ===
        s = H_priv_g2 @ v_F2m   # shape (2t,) in F_{2^m}
        
        # === Step 4 : Berlekamp-Massey to find the error locator polynomial σ(X) ===
        sigma = gf.berlekamp_massey(s)
        
        # === Step 5 : Evaluate σ on L to recover e ===
        # Convention : σ(α_i) = 0  ⟺  i is an error position
        sigma_vals = sigma(self.L)
        e_bits = np.array(sigma_vals == 0, dtype=int)
        weight = int(np.sum(e_bits))
        
        if weight != self.t:
            raise RuntimeError(
                f"Decoding failure: found {weight} roots of σ in L, expected {self.t}. "
                f"This indicates a corrupted ciphertext or > t errors."
            )
        
        e = self.GF2(e_bits)
        
        # === Step: Hamming weight matrix for side-channel analysis ===
        HW_matrix = np.array(H_priv_g2.vector(), dtype=int).sum(axis=2)
        
        return e, HW_matrix
