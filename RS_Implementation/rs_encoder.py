#!/usr/bin/env python3
# rs_encoder.py
# Systematic Reed–Solomon encoder (GF(2^8), primitive 0x11d)

from typing import List

PRIM = 0x11d
FIELD_SIZE = 256
GENERATOR = 2

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

    def sub(self, a: int, b: int) -> int:
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

    def pow(self, a: int, p: int) -> int:
        if a == 0:
            return 0
        return self.exp[(self.log[a] * p) % self.n]

    def alpha_pow(self, p: int) -> int:
        return self.exp[p % self.n]

# polynomial utilities (highest-degree first)
def poly_trim(p: List[int]) -> List[int]:
    i = 0
    while i < len(p)-1 and p[i] == 0:
        i += 1
    return p[i:]

def poly_add(a: List[int], b: List[int], gf: GF) -> List[int]:
    # highest-first representation
    if len(a) < len(b):
        a, b = b, a
    # pad b
    b = [0] * (len(a)-len(b)) + b
    return [gf.add(x,y) for x,y in zip(a,b)]

def poly_mul(a: List[int], b: List[int], gf: GF) -> List[int]:
    # highest-first
    res = [0]*(len(a)+len(b)-1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            res[i+j] = gf.add(res[i+j], gf.mul(ai, bj))
    return res

def poly_divmod(dividend: List[int], divisor: List[int], gf: GF):
    # both highest-first; return (quotient, remainder) (both highest-first)
    dividend = dividend[:]  # copy
    dividend = poly_trim(dividend)
    divisor = poly_trim(divisor)
    if len(divisor) == 0 or (len(divisor) == 1 and divisor[0] == 0):
        raise ZeroDivisionError("Division by zero polynomial")
    if len(dividend) < len(divisor):
        return [0], dividend
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    rem = dividend[:]
    for i in range(len(quotient)):
        coef = rem[i]
        if coef != 0:
            factor = gf.div(coef, divisor[0])
            quotient[i] = factor
            for j in range(len(divisor)):
                rem[i+j] = gf.sub(rem[i+j], gf.mul(divisor[j], factor))
    remainder = rem[-(len(divisor)-1):] if len(divisor) > 1 else [0]
    return poly_trim(quotient), poly_trim(remainder)

def rs_generator_poly(nsym: int, gf: GF) -> List[int]:
    # g(x) = (x - alpha^1)(x - alpha^2)...(x - alpha^nsym)
    g = [1]
    for i in range(nsym):
        # multiply g by [1, alpha^ (i) ]? careful: we want x - alpha^{i+1}
        # highest-first representation: [1, -alpha^k] but since in GF subtraction == addition:
        factor = [1, gf.alpha_pow(i+1)]  # corresponds to x + alpha^{i+1} but sign is XOR (same)
        g = poly_mul(g, factor, gf)
    return poly_trim(g)

def rs_encode_msg(msg: List[int], nsym: int, prim: int = PRIM, verbose: bool = False) -> List[int]:
    """
    Systematic RS encoder: returns message + parity (each symbol 0..255).
    msg: list of ints length k
    nsym: number of parity symbols
    verbose: if True, print internal steps for debugging/learning
    """
    gf = GF(prim)
    if len(msg) + nsym > 255:
        raise ValueError("Message too long for RS(255,k)")

    # 1) Build generator polynomial
    gen = rs_generator_poly(nsym, gf)  # highest-first
    if verbose:
        print("=== ENCODER STEP: Generator polynomial (highest-first coeffs) ===")
        print(gen)

    # 2) Form message polynomial shifted by x^nsym (append nsym zeros)
    msg_poly = msg[:] + [0] * nsym  # highest-first representation
    if verbose:
        print("=== ENCODER STEP: Message polynomial shifted by x^nsym ===")
        print("msg_poly (highest-first):", msg_poly)

    # 3) Polynomial division: (msg(x) * x^nsym) // g(x)
    quotient, remainder = poly_divmod(msg_poly, gen, gf)
    if verbose:
        print("=== ENCODER STEP: Polynomial division result ===")
        print("quotient (highest-first):", quotient)
        print("remainder (highest-first):", remainder)

    # 4) Normalize remainder length to nsym
    if len(remainder) < nsym:
        rem = [0] * (nsym - len(remainder)) + remainder
    else:
        rem = remainder[-nsym:]
    if verbose:
        print("=== ENCODER STEP: Parity (padded to nsym) ===")
        print("parity bytes (highest-first):", rem)

    # 5) Build systematic codeword: message + parity
    codeword = msg + rem
    if verbose:
        print("=== ENCODER STEP: Final codeword (systematic) ===")
        print("codeword:", codeword)

    # ensure byte values are 0..255
    return [x & 0xFF for x in codeword]

if __name__ == "__main__":
    msg = [32,91,11,120,209]
    nsym = 4
    cw = rs_encode_msg(msg, nsym)
    print("Encoded cw:", cw)