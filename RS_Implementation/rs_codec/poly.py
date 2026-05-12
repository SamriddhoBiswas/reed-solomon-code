# rs_codec/poly.py
# Polynomial helpers using highest-degree-first representation.
from typing import List, Tuple
from .gf import GF

def trim(p: List[int]) -> List[int]:
    """Remove leading zeros (highest-first). Keep at least one term."""
    i = 0
    while i < len(p) - 1 and p[i] == 0:
        i += 1
    return p[i:]

def add(a: List[int], b: List[int], gf: GF) -> List[int]:
    """Polynomial addition (XOR) highest-first."""
    if len(a) < len(b):
        a, b = b, a
    b = [0] * (len(a) - len(b)) + b
    return [gf.add(x, y) for x, y in zip(a, b)]

def mul(a: List[int], b: List[int], gf: GF) -> List[int]:
    """Polynomial multiplication (highest-first)."""
    if not a or not b:
        return [0]
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j, bj in enumerate(b):
            if bj == 0:
                continue
            res[i + j] = gf.add(res[i + j], gf.mul(ai, bj))
    return trim(res)

def divmod_poly(dividend: List[int], divisor: List[int], gf: GF) -> Tuple[List[int], List[int]]:
    """Polynomial long division (both highest-first). Returns (quotient, remainder)."""
    dividend = trim(dividend[:])
    divisor = trim(divisor[:])
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
                rem[i + j] = gf.sub(rem[i + j], gf.mul(divisor[j], factor))
    remainder = rem[-(len(divisor) - 1):] if len(divisor) > 1 else [0]
    return trim(quotient), trim(remainder)
