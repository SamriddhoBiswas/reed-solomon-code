# rs_codec/gf.py
# GF(2^8) arithmetic (primitive polynomial 0x11d)
from typing import List

PRIM = 0x11d
FIELD_SIZE = 256
GENERATOR = 2

class GF:
    def __init__(self, prim: int = PRIM):
        self.prim = prim
        self.n = 255
        # exp/log tables for fast mul/div
        self.exp: List[int] = [0] * (2 * self.n + 1)
        self.log: List[int] = [-1] * (self.n + 1)
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

    def alpha_pow(self, p: int) -> int:
        return self.exp[p % self.n]
