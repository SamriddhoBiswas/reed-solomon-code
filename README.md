# Reed-Solomon Code: Implementation and Decomposition Analysis

**Institute:** Indian Institute of Information Technology, Kalyani  
**Project Title:** Decomposition Analysis of Reed-Solomon Codes  
**Supervisor:** Dr. Bhaskar Biswas

**Team Members:**

| Name | Roll Number |
|------|-------------|
| Raman Sah | CSE/23075/1123 |
| Ritaban Chaudhuri | CSE/23079/1127 |
| Samriddho Biswas | CSE/23085/1133 |
| Swastik Gupta | CSE/23104/1152 |

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Repository Structure](#2-repository-structure)
3. [Mathematical Background](#3-mathematical-background)
4. [Module: RS_Implementation](#4-module-rs_implementation)
   - 4.1 [rs_codec Package](#41-rs_codec-package)
   - 4.2 [rs_encoder.py (Standalone Encoder)](#42-rs_encoderpy-standalone-encoder)
   - 4.3 [rs_bm_forney.py (Decoder)](#43-rs_bm_forneypy-decoder)
   - 4.4 [Test Scripts](#44-test-scripts)
5. [Module: RS_Decomposition](#5-module-rs_decomposition)
   - 5.1 [Decomposition Theory](#51-decomposition-theory)
   - 5.2 [rs_splitter.py](#52-rs_splitterpy)
6. [Running the Code](#6-running-the-code)
7. [Dependencies](#7-dependencies)
8. [Key Results](#8-key-results)
9. [References](#9-references)

---

## 1. Project Overview

This repository contains a complete Python implementation of Reed-Solomon (RS) error-correcting codes, together with a systematic analysis of their decomposability into smaller valid sub-codes.

A Reed-Solomon code is denoted RS(n, k, t) where:

- `n` = total codeword length (number of symbols)
- `k` = number of data symbols
- `t` = number of correctable symbol errors

and the fundamental identity `n = k + 2t` holds throughout.

The project is divided into two independent components:

- **RS_Implementation**: A ground-up implementation of RS encoding over GF(2^8) and decoding via the Berlekamp-Massey and Forney algorithms.
- **RS_Decomposition**: A framework that determines whether any RS(n, k, t) code can be split into two smaller valid RS codes, derives a closed-form theorem for this condition, and enumerates all valid splits.

The central theorem derived and proved in this project is:

> RS(n, k, t) is decomposable into two valid RS sub-codes if and only if `k >= 2` and `t >= 2`.

---

## 2. Repository Structure

```
reed-solomon-code-main/
|
|-- RS_Implementation/                  # Encoder and decoder implementation
|   |-- rs_codec/                       # Reusable codec package
|   |   |-- __init__.py                 # Exports RSCodec, rs_encode_msg
|   |   |-- gf.py                       # GF(2^8) finite field arithmetic
|   |   |-- poly.py                     # Polynomial operations over GF
|   |   |-- generator.py                # RS generator polynomial construction
|   |   |-- encoder.py                  # High-level systematic RS encoder (class + function)
|   |
|   |-- rs_encoder.py                   # Standalone encoder with verbose debugging
|   |-- rs_bm_forney.py                 # Berlekamp-Massey + Forney decoder
|   |-- test_encode.py                  # Basic encoding smoke test
|   |-- test_full_cycle_bm.py           # End-to-end encode-corrupt-decode test suite
|
|-- RS_Decomposition/                   # Decomposition analysis framework
|   |-- rs_splitter.py                  # Full decomposition engine (brute-force + O(1) theorem)
|
|-- README.md                           
|-- .gitignore
|-- .gitattributes
```

---

## 3. Mathematical Background

### Finite Field GF(2^8)

All arithmetic is performed over GF(2^8), the Galois field of 256 elements. The primitive polynomial used is:

```
p(x) = x^8 + x^4 + x^3 + x^2 + 1   (hexadecimal: 0x11d)
```

Multiplication is implemented using precomputed exponent and logarithm tables (anti-log tables), making field multiplication O(1). Since GF(2^m) has characteristic 2, addition and subtraction are both XOR.

### RS Code Structure

An RS code over GF(2^8) with `t` parity symbols has a generator polynomial:

```
g(x) = (x - alpha^1)(x - alpha^2) ... (x - alpha^(2t))
```

where `alpha` is a primitive element of GF(2^8). Encoding is systematic: the codeword takes the form `[message | parity]`, where the parity symbols are the remainder of `(message(x) * x^(2t)) / g(x)`.

### Decoding (Berlekamp-Massey + Forney)

Error correction follows the standard algebraic decoding pipeline:

1. **Syndrome computation**: Evaluate the received polynomial at `alpha^1` through `alpha^(2t)`.
2. **Berlekamp-Massey algorithm**: Use the syndrome sequence to find the error-locator polynomial `sigma(x)`.
3. **Chien search**: Find the roots of `sigma(x)` over GF(2^8) to identify error positions.
4. **Forney algorithm**: Compute the error magnitudes from the error-evaluator polynomial `Omega(x)`.
5. **Correction**: XOR the received codeword at the identified positions with the computed magnitudes.

---

## 4. Module: RS_Implementation

### 4.1 rs_codec Package

The `rs_codec` package provides a clean, importable encoder. It is structured as four submodules:

#### `rs_codec/gf.py` — Galois Field Arithmetic

Implements the `GF` class, which builds exponent and logarithm tables for GF(2^8) at construction time.

| Method | Description |
|--------|-------------|
| `add(a, b)` | XOR (addition in GF(2^8)) |
| `sub(a, b)` | XOR (subtraction = addition in characteristic 2) |
| `mul(a, b)` | Table-lookup multiplication |
| `div(a, b)` | Table-lookup division |
| `inv(a)` | Multiplicative inverse |
| `alpha_pow(p)` | Returns `alpha^p`, the p-th power of the primitive element |

The table construction runs once in `__init__` and generates `exp[0..510]` and `log[0..255]`.

#### `rs_codec/poly.py` — Polynomial Operations

All polynomials are represented with the highest-degree coefficient first (big-endian).

| Function | Description |
|----------|-------------|
| `trim(p)` | Strip leading zero coefficients |
| `add(a, b, gf)` | Coefficient-wise XOR after zero-padding the shorter polynomial |
| `mul(a, b, gf)` | Convolution using GF multiplication |
| `divmod_poly(dividend, divisor, gf)` | Polynomial long division; returns `(quotient, remainder)` |

#### `rs_codec/generator.py` — Generator Polynomial

`rs_generator_poly(nsym, gf)` builds `g(x)` by successively multiplying linear factors `(x - alpha^i)` for `i = 1` to `nsym`. The result is a polynomial of degree `nsym` with all coefficients in GF(2^8).

#### `rs_codec/encoder.py` — Systematic Encoder

Provides the `RSCodec` class and a convenience function.

```python
from rs_codec import RSCodec, rs_encode_msg

# Class-based
codec = RSCodec(nsym=4)
codeword = codec.encode([32, 91, 11, 120, 209])

# Functional
codeword = rs_encode_msg([32, 91, 11, 120, 209], nsym=4)
```

Encoding steps performed internally:
1. Validate message length (`k + nsym <= 255`).
2. Build or retrieve cached generator polynomial.
3. Compute parity: divide `message * x^nsym` by `g(x)`, take the remainder.
4. Concatenate message and parity bytes to form the systematic codeword.

Shortened codes are supported by default; pass `shorten=False` to obtain the full 255-symbol codeword with zero padding.

---

### 4.2 rs_encoder.py (Standalone Encoder)

`rs_encoder.py` is a self-contained single-file encoder that mirrors the `rs_codec` package but adds a `verbose` flag. When `verbose=True`, it prints each intermediate step of the encoding process — generator polynomial, message polynomial, division quotient and remainder, and final codeword — making it useful for learning and debugging.

```python
from rs_encoder import rs_encode_msg

codeword = rs_encode_msg([32, 91, 11, 120, 209], nsym=4, verbose=True)
```

It can also be run directly:

```bash
python3 RS_Implementation/rs_encoder.py
```

---

### 4.3 rs_bm_forney.py (Decoder)

`rs_bm_forney.py` implements the full algebraic decoding pipeline. All major steps are individually implemented.

| Function | Description |
|----------|-------------|
| `compute_syndromes(received, nsym, gf)` | Returns the syndrome vector `[S1, ..., S_nsym]` |
| `berlekamp_massey(synd, gf)` | Returns error-locator polynomial `sigma(x)` (constant-first) |
| `chien_search(locator, gf, n)` | Finds roots of `sigma(x)` by evaluating at all field elements |
| `forney_evaluator(synd, locator, gf, nsym)` | Computes `Omega(x) = (S(x) * sigma(x)) mod x^nsym` |
| `forney(omega, locator, err_pos, gf)` | Applies Forney formula to compute error magnitudes |
| `rs_bm_forney_decode(received, nsym, ...)` | Top-level decode function; returns `(corrected, info_dict)` |

The `info_dict` returned by `rs_bm_forney_decode` contains:

```python
{
  'syndromes': [...],
  'locator': [...],
  'error_positions': [...],
  'error_magnitudes': [...],
  'corrected': True or False
}
```

The decoder supports a `verbose=True` flag that prints a step-by-step trace of the decoding process.

---

### 4.4 Test Scripts

#### `test_encode.py`

A minimal smoke test that encodes a fixed message using both the `RSCodec` class and the `rs_encode_msg` function, printing the resulting codewords for visual verification.

```bash
python3 RS_Implementation/test_encode.py
```

#### `test_full_cycle_bm.py`

A comprehensive end-to-end test that:

1. Accepts a string from stdin and converts it to a byte array.
2. Encodes the message using the RS encoder (`nsym = 4`, capable of correcting up to 2 symbol errors).
3. Introduces 2 artificial symbol errors and attempts decoding.
4. Runs a statistical sweep: for each error count from 0 to `nsym+1`, performs 100 random trials and records the success rate.

```bash
python3 RS_Implementation/test_full_cycle_bm.py
# Enter a message string: Hello
```

Expected behavior: decoding succeeds with probability 1.0 for up to 2 errors (`t = nsym/2 = 2`) and drops to 0 for 3 or more errors.

---

## 5. Module: RS_Decomposition

### 5.1 Decomposition Theory

A valid 2-way decomposition of RS(n, k, t) requires finding RS(n1, k1, t1) and RS(n2, k2, t2) such that:

```
n = n1 + n2
t = t1 + t2
k = k1 + k2   (automatically satisfied given the above and the RS identity)
```

and each sub-code independently satisfies `ni = ki + 2*ti`, `ki >= 1`, `ti >= 1`, `ni >= 3`.

The key insight is that once the split of `t` and `n` is fixed, the split of `k` is fully determined: `ki = ni - 2*ti`. This reduces the search space to pairs `(t1, n1)` only.

The derived theorem, validated by exhaustive brute-force over all codes with `n <= 25`:

```
RS(n, k, t) is decomposable  <=>  k >= 2  AND  t >= 2
```

**Proof sketch (sufficiency):** Set `t1 = 1`, `t2 = t - 1`. A valid `n1` must satisfy `3 <= n1 <= k + 1`. This interval is non-empty if and only if `k + 1 >= 3`, i.e., `k >= 2`.

**Proof sketch (necessity):** If `t < 2`, no valid split of `t` into two positive integers exists. If `k < 2`, the interval for `n1` is empty and no valid split exists.

---

### 5.2 rs_splitter.py

`rs_splitter.py` is a self-contained module that provides three layers of functionality.

#### Core Data Type

```python
@dataclass(frozen=True)
class Split:
    n1: int; k1: int; t1: int
    n2: int; k2: int; t2: int
```

#### Key Functions

| Function | Complexity | Description |
|----------|------------|-------------|
| `is_valid_rs(n, k, t)` | O(1) | Returns True if RS(n, k, t) is a valid code |
| `validate_parent(n, k, t)` | O(1) | Returns `(ok, error_message)` for the parent code |
| `can_split(n, k, t)` | O(1) | Applies the theorem; returns `(bool, explanation)` |
| `find_all_splits(n, k, t)` | O(n * t) | Enumerates all valid 2-way splits |
| `print_splits(n, k, t, splits)` | — | Formatted table output with rates and field sizes |
| `run_tests()` | — | 51-case automated test suite |
| `demo()` | — | Demonstration on standard and real-world codes |

#### Algorithm (find_all_splits)

```
for t1 in {1, ..., t-1}:
    t2 = t - t1
    for n1 in {2*t1 + 1, ..., n - 2*t2 - 1}:
        n2 = n - n1
        k1 = n1 - 2*t1
        k2 = n2 - 2*t2
        yield Split(n1, k1, t1, n2, k2, t2)
```

Every `(t1, n1)` pair within the stated bounds yields a valid split with no additional validity checks required, because the bounds are derived precisely from `k1 >= 1` and `k2 >= 1`.

#### Command-Line Usage

```bash
# Run built-in tests and full demonstration
python3 RS_Decomposition/rs_splitter.py

# Analyze a specific code from the command line
python3 RS_Decomposition/rs_splitter.py 15 9 3

# Analyze the CD/DVD code
python3 RS_Decomposition/rs_splitter.py 255 223 16
```

#### Import Usage

```python
from RS_Decomposition.rs_splitter import find_all_splits, can_split, print_splits

splittable, reason = can_split(255, 223, 16)
splits = find_all_splits(255, 223, 16)
print_splits(255, 223, 16, splits)
```

---

## 6. Running the Code

Ensure Python 3.8 or later is installed. No third-party packages are required.

### Encoding only

```bash
cd RS_Implementation
python3 test_encode.py
```

### Full encode-decode cycle

```bash
cd RS_Implementation
python3 test_full_cycle_bm.py
```

### Standalone encoder with verbose trace

```bash
cd RS_Implementation
python3 rs_encoder.py
```

### Decomposition analysis (tests + demo)

```bash
cd RS_Decomposition
python3 rs_splitter.py
```

### Decomposition for a custom code

```bash
python3 RS_Decomposition/rs_splitter.py <n> <k> <t>
# Example:
python3 RS_Decomposition/rs_splitter.py 31 15 8
```

---

## 7. Dependencies

| Requirement | Version | Purpose |
|-------------|---------|---------|
| Python | >= 3.8 | All modules |
| Standard library only | — | `sys`, `dataclasses`, `typing`, `random`, `json`, `pprint` |

No external packages (NumPy, SciPy, etc.) are required. All finite field and polynomial arithmetic is implemented from scratch.

---

## 8. Key Results

### Theorem

RS(n, k, t) is decomposable if and only if `k >= 2` and `t >= 2`.

### Test Results (rs_splitter.py)

| Metric | Value |
|--------|-------|
| Total test cases | 51 |
| Passed | 51 |
| Failed | 0 |
| Theorem-brute-force agreement (all n <= 25) | 100% |

### Real-World Code Analysis

| Application | RS Code | Decomposable |
|-------------|---------|--------------|
| CD / DVD storage | RS(255, 223, 16) | Yes |
| DVB-S satellite | RS(255, 239, 8) | Yes |
| DVB-T terrestrial | RS(204, 188, 8) | Yes |
| NASA (TENEX) | RS(31, 15, 8) | Yes |

Example decomposition for the CD/DVD code:

```
RS(255, 223, 16) = RS(3, 1, 1) + RS(252, 222, 15)
```

### Algorithm Complexity

| Component | Time Complexity |
|-----------|----------------|
| GF table construction | O(n) where n = 255 |
| Encoding | O(k * nsym) |
| Syndrome computation | O(n * nsym) |
| Berlekamp-Massey | O(nsym^2) |
| Chien search | O(n * t) |
| Forney evaluation | O(t^2) |
| Decomposition check (can_split) | O(1) |
| Full split enumeration (find_all_splits) | O(n * t) |

---

## 9. References

1. I. S. Reed and G. Solomon, "Polynomial Codes Over Certain Finite Fields," *Journal of the Society for Industrial and Applied Mathematics*, vol. 8, no. 2, pp. 300-304, 1960.
2. R. E. Blahut, *Theory and Practice of Error Control Codes*. Addison-Wesley, 1983.
3. S. Lin and D. J. Costello, *Error Control Coding: Fundamentals and Applications*, 2nd ed. Prentice Hall, 2004.
4. S. B. Wicker and V. K. Bhargava, *Reed-Solomon Codes and Their Applications*. IEEE Press, 1994.
5. W. W. Peterson and E. J. Weldon, *Error-Correcting Codes*, 2nd ed. MIT Press, 1972.
