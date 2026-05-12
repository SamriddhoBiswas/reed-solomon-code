#!/usr/bin/env python3
# test_full_cycle_bm.py
from rs_encoder import rs_encode_msg
from rs_bm_forney import rs_bm_forney_decode
import random, json
from pprint import pprint

# configure
#message = [32, 91, 11, 120, 209]  # message (highest-degree first)
user_input = input("Enter a message string: ")
message = [ord(c) for c in user_input]  # Convert string to array of byte values
nsym = 4  # parity symbols
k = len(message)
n = k + nsym

print("Original message:", message)

# 1) Encode
codeword = rs_encode_msg(message, nsym)
print("Encoded message:", codeword)

# 2) Demo: add 2 symbol errors
rx = codeword[:]
rx[2] ^= 5
rx[5] ^= 3
print("Received:", rx)
corrected, info = rs_bm_forney_decode(rx, nsym)
print("\nDecoder info:")
pprint(info)
print("Corrected:", corrected)
print("Recovered message:", corrected[:k] if info.get('corrected') else None)

recovered_string = ''.join(chr(byte) for byte in corrected[:k])
print("Recovered string:", recovered_string)

# 3) Full sweep: test 100 trials per error count
trials = 100
results = {}
for ecount in range(0, nsym+2):  # 0..nsym+1
    succ = 0
    for _ in range(trials):
        rx = codeword[:]
        if ecount > 0:
            pos = random.sample(range(n), ecount)
            for p in pos:
                rx[p] ^= random.randint(1,255)
        corr, inf = rs_bm_forney_decode(rx, nsym)
        if inf.get('corrected') and corr[:k] == message:
            succ += 1
    results[ecount] = {'successes': succ, 'trials': trials, 'rate': succ / trials}
print("\nResults (formatted):")
print(json.dumps(results, indent=2))