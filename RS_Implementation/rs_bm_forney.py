#!/usr/bin/env python3
# rs_bm_forney.py
# RS decoder based on Berlekamp–Massey + Forney (GF(2^8), prim 0x11d)
from typing import List, Tuple, Optional
import json
from pprint import pprint

PRIM = 0x11d

class GF:
    def __init__(self, prim=PRIM):
        self.prim = prim
        self.n = 255
        self.exp = [0] * (2 * self.n + 1)
        self.log = [-1] * (self.n + 1)
        x = 1
        for i in range(self.n):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & 0x100:
                x ^= self.prim
        for i in range(self.n, 2 * self.n + 1):
            self.exp[i] = self.exp[i - self.n]

    def add(self, a: int, b: int) -> int:
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return self.exp[(self.log[a] + self.log[b]) % self.n]

    def div(self, a: int, b: int) -> int:
        assert b != 0
        if a == 0:
            return 0
        return self.exp[(self.log[a] - self.log[b]) % self.n]

    def inv(self, a: int) -> int:
        assert a != 0
        return self.exp[(self.n - self.log[a]) % self.n]

    def alpha_pow(self, p: int) -> int:
        return self.exp[p % self.n]

def poly_trim(p: List[int]) -> List[int]:
    i = 0
    while i < len(p)-1 and p[i] == 0:
        i += 1
    return p[i:]

# We'll use highest-first representation for polynomials where convenient.
# But BM commonly uses coefficients in ascending order (constant term first).
# We'll convert syndrome list into ascending order for BM.

def compute_syndromes(received: List[int], nsym: int, gf: GF) -> List[int]:
    # Compute S1..S_nsym, S_j = r(alpha^j) for j=1..nsym
    n = len(received)
    S = []
    for j in range(1, nsym+1):
        s = 0
        for i, val in enumerate(received):
            if val == 0:
                continue
            # For systematic codeword: ordering is msg (highest coefficients) + parity
            # We evaluate r at alpha^j using index mapping: position i -> power (n-1 - i)
            power = (n - 1 - i) % gf.n
            a = gf.alpha_pow((power * j) % gf.n)
            s ^= gf.mul(val, a)
        S.append(s)
    return S

def berlekamp_massey(synd: List[int], gf: GF) -> List[int]:
    # synd list S1..S2t (highest to lowest j), produce locator polynomial coefficients (constant-first)
    # Implementation follows standard BM producing sigma(x) with sigma[0]=1
    n = len(synd)
    C = [1] + [0]*n
    B = [1] + [0]*n
    L = 0
    m = 1
    b = 1
    for i in range(n):
        # compute discrepancy
        d = synd[i]
        for j in range(1, L+1):
            if C[j] != 0 and synd[i-j] != 0:
                d ^= gf.mul(C[j], synd[i-j])
        if d == 0:
            m += 1
        else:
            T = C[:]
            coef = gf.div(d, b)
            # C = C - coef * x^m * B
            for j in range(0, n - (m) + 1):
                if B[j] != 0:
                    C[j + m] ^= gf.mul(coef, B[j])
            if 2*L <= i:
                L2 = i + 1 - L
                B = T
                b = d
                L = i + 1 - L
                m = 1
            else:
                m += 1
    # trim C to degree L
    return C[:L+1]

def chien_search(locator: List[int], gf: GF, n: int) -> List[int]:
    # locator is constant-first (sigma[0]=1)
    roots = []
    for i in range(n):
        # evaluate sigma(alpha^{-i})
        x = gf.alpha_pow((-i) % gf.n)
        val = 0
        powx = 1
        for coef in locator:
            if coef != 0:
                val ^= gf.mul(coef, powx)
            powx = gf.mul(powx, x)
        if val == 0:
            roots.append(i)
    return roots

def poly_mul(a: List[int], b: List[int], gf: GF) -> List[int]:
    # a,b constant-first
    res = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            res[i+j] ^= gf.mul(ai, bj)
    return res

def poly_scale(a: List[int], scalar: int, gf: GF) -> List[int]:
    return [gf.mul(coef, scalar) for coef in a]

def poly_divmod_constant_first(a: List[int], b: List[int], gf: GF) -> Tuple[List[int], List[int]]:
    # Convert to highest-first for division then back
    def to_highest(p):
        return list(reversed(p))
    def to_lowest(p):
        return list(reversed(p))
    ah = to_highest(a)
    bh = to_highest(b)
    # long division highest-first:
    qh = [0]*(len(ah)-len(bh)+1) if len(ah) >= len(bh) else [0]
    rh = ah[:]
    if len(ah) < len(bh):
        return [0], a[:]
    for i in range(len(qh)):
        coef = rh[i]
        if coef != 0:
            factor = gf.div(coef, bh[0])
            qh[i] = factor
            for j in range(len(bh)):
                rh[i+j] ^= gf.mul(bh[j], factor)
    # remainder:
    rem = rh[-(len(bh)-1):] if len(bh) > 1 else [0]
    return to_lowest(qh), to_lowest(rem)

def forney_evaluator(synd: List[int], locator: List[int], gf: GF, nsym: int) -> List[int]:
    # compute Omega(x) = (S(x) * Lambda(x)) mod x^{nsym}
    # represent S(x) as polynomial with coefficients S1..Snsym in ascending order for powers x^0.. (S1 is coeff of x^0)
    S_poly = synd[:]  # S1..Snsym => S[0] is coeff for x^0
    Omega = poly_mul(S_poly, locator, gf)
    # take first nsym coefficients
    return Omega[:nsym]

def forney(omega: List[int], locator: List[int], err_pos: List[int], gf: GF) -> List[int]:
    # compute error magnitudes using Forney formula:
    # E_i = - Omega(alpha^{-pos}) / Lambda'(alpha^{-pos})
    # locator and omega are constant-first
    err_mags = []
    # derivative of locator
    deriv = []
    for i in range(1, len(locator)):
        if i % 2 == 1:
            deriv.append(locator[i])
        else:
            deriv.append(0)
    # deriv is constant-first but one degree shorter
    for pos in err_pos:
        x = gf.alpha_pow((-pos) % gf.n)
        # evaluate Omega at x
        num = 0
        powx = 1
        for coef in omega:
            if coef != 0:
                num ^= gf.mul(coef, powx)
            powx = gf.mul(powx, x)
        # evaluate deriv at x
        den = 0
        powx = 1
        for coef in deriv:
            if coef != 0:
                den ^= gf.mul(coef, powx)
            powx = gf.mul(powx, x)
        if den == 0:
            return None
        val = gf.div(num, den)
        # error value is -val, but subtraction == addition in GF(2^m)
        err_mags.append(val)
    return err_mags

def rs_bm_forney_decode(received: List[int], nsym: int, prim: int = PRIM, verbose: bool = False) -> Tuple[List[int], dict]:
    gf = GF(prim)
    n = len(received)
    info = {'syndromes': None, 'locator': None, 'error_positions': [], 'error_magnitudes': [], 'corrected': False}
    if verbose:
        print(f"Starting RS decoding with {nsym} syndromes on received codeword of length {n}")
    if verbose:
        print("Computing syndromes...")
    S = compute_syndromes(received, nsym, gf)
    info['syndromes'] = S
    if verbose:
        print(f"Syndromes computed: {S}")
    if all(x == 0 for x in S):
        if verbose:
            print("All syndromes zero, no errors detected.")
        info['corrected'] = True
        return received[:], info
    if verbose:
        print("Running Berlekamp-Massey algorithm to find error locator polynomial...")
    sigma = berlekamp_massey(S, gf)  # sigma constant-first
    info['locator'] = sigma
    if verbose:
        print(f"Error locator polynomial coefficients (constant-first): {sigma}")
    if verbose:
        print("Performing Chien search to find error positions...")
    roots = chien_search(sigma, gf, n)
    if verbose:
        print(f"Chien search found error positions (roots): {roots}")
    if len(roots) == 0:
        if verbose:
            print("No error positions found, decoding failed.")
        return received[:], info
    if verbose:
        print("Computing error evaluator polynomial Omega(x)...")
    omega = forney_evaluator(S, sigma, gf, nsym)
    if verbose:
        print(f"Error evaluator polynomial Omega coefficients: {omega}")
    if verbose:
        print("Computing error magnitudes using Forney's algorithm...")
    errs = forney(omega, sigma, roots, gf)
    if errs is None:
        if verbose:
            print("Error magnitudes could not be computed (division by zero), decoding failed.")
        return received[:], info
    if verbose:
        print(f"Error magnitudes computed: {errs}")
    corrected = received[:]
    # in our syndrome eval we used position -> power = (n-1 - i)
    # chien_search returned pos such that alpha^{-pos} is root; this pos corresponds to power index
    # mapping from pos to index i: power = pos => index = n-1 - pos
    mapped_positions = [(n - 1 - p) % n for p in roots]
    if verbose:
        print(f"Mapped error positions to received indices: {mapped_positions}")
        print("Applying corrections to received codeword...")
    for idx, p in enumerate(mapped_positions):
        if verbose:
            print(f"Correcting position {p} with error magnitude {errs[idx]}")
        corrected[p] ^= errs[idx]
    if verbose:
        print("Recomputing syndromes after correction to verify...")
    S2 = compute_syndromes(corrected, nsym, gf)
    if verbose:
        print(f"Syndromes after correction: {S2}")
    if all(x == 0 for x in S2):
        if verbose:
            print("Syndromes zero after correction, decoding successful.")
        info['corrected'] = True
        info['error_positions'] = mapped_positions
        info['error_magnitudes'] = errs
        return corrected, info
    if verbose:
        print("Syndromes not zero after correction, decoding failed.")
    return received[:], info

if __name__ == "__main__":
    print("rs_bm_forney module. Import rs_bm_forney_decode(received, nsym).")