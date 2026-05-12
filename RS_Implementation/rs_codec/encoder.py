# rs_codec/encoder.py
# High-level RS encoder (systematic). External-facing API.
from typing import List
from .gf import GF
from .generator import rs_generator_poly 
from .poly import divmod_poly, trim

class RSCodec:
    """
    RS encoder instance.

    Usage:
        codec = RSCodec(nsym=4)
        cw = codec.encode(msg)           # returns msg + parity (shortened)
        parity = codec.parity_for(msg)   # returns parity bytes only
    """

    def __init__(self, nsym: int, prim: int = None):
        self.nsym = int(nsym)
        self.gf = GF(prim) if prim is not None else GF()
        self.n = self.gf.n  # 255
        self._gen = None

    def generator(self) -> List[int]:
        if self._gen is None:
            self._gen = rs_generator_poly(self.nsym, self.gf)
        return self._gen

    def parity_for(self, msg: List[int]) -> List[int]:
        """
        Return parity bytes for the given message.
        msg: highest-degree-first list of ints (0..255).
        This uses canonical shortening: it encodes as full-255 and returns parity for that mapping.
        """
        k = len(msg)
        if k + self.nsym > self.n:
            raise ValueError("Message too long for RS(255)")
        # Left-pad so message sits at right of full 255-length polynomial
        pad_len = self.n - (k + self.nsym)
        msg_full = [0] * pad_len + msg
        # Multiply by x^nsym: append nsym zeros
        msg_poly = msg_full + [0] * self.nsym
        gen = self.generator()
        _, remainder = divmod_poly(msg_poly, gen, self.gf)
        # normalize remainder length to nsym
        if len(remainder) < self.nsym:
            rem = [0] * (self.nsym - len(remainder)) + remainder
        else:
            rem = remainder[-self.nsym:]
        return rem

    def encode(self, msg: List[int], shorten: bool = True) -> List[int]:
        """
        Return systematic codeword (message + parity).
        If shorten=True (default), returns right-most k+nsym bytes (shortened).
        If shorten=False, returns full 255-length codeword (msg left-padded).
        """
        if any((not isinstance(x, int)) or x < 0 or x > 255 for x in msg):
            raise ValueError("msg values must be ints in range 0..255")
        k = len(msg)
        if k + self.nsym > self.n:
            raise ValueError("Message too long for RS(255)")
        pad_len = self.n - (k + self.nsym)
        # compute parity for canonical full-length mapping
        parity = self.parity_for(msg)
        # build full codeword
        msg_full = [0] * pad_len + msg
        cw_full = msg_full + parity
        if shorten:
            return cw_full[-(k + self.nsym):]
        return cw_full

# convenience functional API
def rs_encode_msg(msg: List[int], nsym: int) -> List[int]:
    codec = RSCodec(nsym)
    return codec.encode(msg)
