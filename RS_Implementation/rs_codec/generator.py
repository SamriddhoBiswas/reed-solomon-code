# rs_codec/generator.py
# Build the RS generator polynomial g(x) = (x - alpha^1)(x - alpha^2)...(x - alpha^nsym)
from typing import List
from .gf import GF
from .poly import mul, trim

def rs_generator_poly(nsym: int, gf: GF) -> List[int]:
    g = [1]
    for i in range(nsym):
        # factor corresponds to (x - alpha^(i+1)) ; subtraction == addition in GF(2^m)
        factor = [1, gf.alpha_pow(i + 1)]
        g = mul(g, factor, gf)
    return trim(g)
