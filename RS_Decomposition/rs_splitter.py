import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
"""
Reed-Solomon Code Splitter
==========================
Given RS(n, k, t) with n = k + 2t, finds every valid 2-way split into
RS(n1, k1, t1)  ⊕  RS(n2, k2, t2)  such that:

    n  = n1 + n2          (block lengths add up)
    t  = t1 + t2          (error-correction capacities add up)
    k  = k1 + k2          (data symbols add up — derived automatically)

Each sub-code must independently satisfy the RS validity conditions:
    ni = ki + 2*ti,   ki >= 1,   ti >= 1,   ni >= 3

Key insight: once n and t are split, k splits are fully determined:
    ki = ni - 2*ti
so the search space is over (n1, t1) pairs only.

Usage
-----
    python3 rs_splitter.py                   # runs built-in demo
    python3 rs_splitter.py 15 9 3            # RS(15,9,3) from command line
    python3 rs_splitter.py 255 223 16        # CD/DVD code

Or import and call directly:
    from rs_splitter import find_all_splits, print_splits
    splits = find_all_splits(15, 9, 3)
    print_splits(15, 9, 3, splits)
"""

import sys
from dataclasses import dataclass
from typing import List, Tuple, Optional


# ─────────────────────────────────────────────────────────────────────────────
# Core data type
# ─────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Split:
    """One valid 2-way split of a parent RS code."""
    # Sub-code 1
    n1: int; k1: int; t1: int
    # Sub-code 2
    n2: int; k2: int; t2: int

    def as_tuple(self):
        return (self.n1, self.k1, self.t1, self.n2, self.k2, self.t2)


# ─────────────────────────────────────────────────────────────────────────────
# Validation helpers
# ─────────────────────────────────────────────────────────────────────────────

def is_valid_rs(n: int, k: int, t: int) -> bool:
    """
    True iff RS(n, k, t) satisfies all Reed-Solomon validity conditions:
      1.  n = k + 2t           (fundamental RS identity)
      2.  k >= 1               (at least one data symbol)
      3.  t >= 1               (corrects at least one error)
      4.  n >= 3               (minimum possible RS code: k=1, t=1 → n=3)
    """
    return (
        n == k + 2 * t  and
        k >= 1          and
        t >= 1          and
        n >= 3
    )


def validate_parent(n: int, k: int, t: int) -> Tuple[bool, str]:
    """
    Validate the parent code parameters and return (ok, error_message).
    Returns (True, '') if valid, (False, reason) if not.
    """
    if n != k + 2 * t:
        return False, (
            f"n = k + 2t violated: {n} ≠ {k} + 2×{t} = {k + 2*t}. "
            f"This is not a valid RS code."
        )
    if k < 1:
        return False, f"k = {k} is invalid: must be ≥ 1."
    if t < 1:
        return False, f"t = {t} is invalid: must be ≥ 1."
    if n < 3:
        return False, f"n = {n} is invalid: must be ≥ 3."
    return True, ""


# ─────────────────────────────────────────────────────────────────────────────
# Main algorithm: find all valid splits
# ─────────────────────────────────────────────────────────────────────────────

def find_all_splits(n: int, k: int, t: int) -> List[Split]:
    """
    Find every valid 2-way split of RS(n, k, t) into two RS sub-codes.

    Algorithm
    ---------
    Iterate over all (t1, n1) pairs:

        t1 ∈ {1, 2, …, t-1}          → t2 = t - t1 ≥ 1
        n1 ∈ {2t1+1, …, n-2t2-1}     → k1 = n1-2t1 ≥ 1  AND  k2 = n2-2t2 ≥ 1

    The inner loop bounds are tight: every (t1, n1) in this range gives a
    valid split with no further checks needed, because:
        k1 = n1 - 2t1 ≥ (2t1+1) - 2t1 = 1  ✓
        k2 = (n-n1) - 2t2 ≥ (2t2+1) - 2t2 = 1  ✓ (when n1 ≤ n-2t2-1)

    Complexity: O(n × t) — at most n×t iterations regardless of result count.

    Parameters
    ----------
    n, k, t : RS code parameters with n = k + 2t

    Returns
    -------
    List of Split objects, in lexicographic order of (t1, n1).
    Raises ValueError if (n, k, t) is not a valid RS code.
    """
    ok, reason = validate_parent(n, k, t)
    if not ok:
        raise ValueError(f"RS({n},{k},{t}) is not a valid RS code: {reason}")

    splits: List[Split] = []

    for t1 in range(1, t):           # t1 in {1, …, t-1}
        t2 = t - t1
        n1_min = 2 * t1 + 1          # smallest n1 that gives k1 = n1-2t1 ≥ 1
        n1_max = n - 2 * t2 - 1      # largest n1 that gives k2 = n2-2t2 ≥ 1
        for n1 in range(n1_min, n1_max + 1):
            n2 = n - n1
            k1 = n1 - 2 * t1
            k2 = n2 - 2 * t2
            splits.append(Split(n1, k1, t1, n2, k2, t2))

    return splits


# ─────────────────────────────────────────────────────────────────────────────
# O(1) theorem: can this code be split at all?
# ─────────────────────────────────────────────────────────────────────────────

def can_split(n: int, k: int, t: int) -> Tuple[bool, str]:
    """
    Decide in O(1) whether any valid split exists, and explain why.

    THEOREM: RS(n,k,t) is splittable  ⟺  k ≥ 2  AND  t ≥ 2.

    Proof sketch:
      Choose t1=1, t2=t-1. Need n1 in [2t1+1, n-2t2-1] = [3, k+1].
      Non-empty iff k+1 ≥ 3 iff k ≥ 2.  ∴ sufficient.
      k<2 or t<2 → no (t1,n1) can satisfy both constraints.  ∴ necessary.
    """
    ok, reason = validate_parent(n, k, t)
    if not ok:
        return False, f"Invalid RS code: {reason}"
    if t < 2:
        return False, (
            f"t = {t} < 2.  Need t ≥ 2 so that t can be split into "
            f"t1 ≥ 1 and t2 ≥ 1."
        )
    if k < 2:
        return False, (
            f"k = {k} < 2.  Need k ≥ 2 so that k = k1+k2 with both k1 ≥ 1 "
            f"and k2 ≥ 1 is possible."
        )
    # Canonical minimum example
    ex_n1, ex_t1 = 3, 1
    ex_n2, ex_t2 = n - 3, t - 1
    ex_k1, ex_k2 = ex_n1 - 2 * ex_t1, ex_n2 - 2 * ex_t2
    return True, (
        f"k = {k} ≥ 2 and t = {t} ≥ 2.  "
        f"Canonical minimum split: RS({ex_n1},{ex_k1},{ex_t1}) ⊕ RS({ex_n2},{ex_k2},{ex_t2})."
    )


# ─────────────────────────────────────────────────────────────────────────────
# Display
# ─────────────────────────────────────────────────────────────────────────────

def _field_size(n: int) -> int:
    """Smallest m such that 2^m - 1 >= n.  GF(2^m) is the natural field."""
    m = 2
    while (1 << m) - 1 < n:
        m += 1
    return m


def print_splits(n: int, k: int, t: int, splits: List[Split],
                 verbose: bool = True) -> None:
    """
    Print a formatted table of all valid splits for RS(n, k, t).
    """
    W = 74
    print()
    print("=" * W)
    print(f"  RS({n}, {k}, {t})  |  d_min = {2*t+1}  |  rate = {k/n:.4f}"
          f"  |  parity symbols = {2*t}")
    splittable, reason = can_split(n, k, t)
    print(f"  Splittable: {'YES' if splittable else 'NO'}  —  {reason}")
    print(f"  Total valid splits: {len(splits)}")
    print("=" * W)

    if not splits:
        print()
        return

    if verbose:
        # Header
        print(f"\n  {'#':>4}  {'Sub-code 1':^24}  {'Sub-code 2':^24}  {'Fields':^14}")
        print(f"  {'':>4}  {'RS(n1,k1,t1)':^24}  {'RS(n2,k2,t2)':^24}  {'GF(2^m1|m2)':^14}")
        print("  " + "-" * (W - 2))

        for i, s in enumerate(splits, 1):
            m1 = _field_size(s.n1)
            m2 = _field_size(s.n2)
            c1 = f"RS({s.n1},{s.k1},{s.t1})"
            c2 = f"RS({s.n2},{s.k2},{s.t2})"
            r1 = f"{s.k1/s.n1:.3f}"
            r2 = f"{s.k2/s.n2:.3f}"
            fields = f"GF(2^{m1}|2^{m2})"
            print(f"  {i:>4}.  {c1:<12}  rate={r1}  {c2:<12}  rate={r2}  {fields}")
    else:
        # Compact
        for i, s in enumerate(splits, 1):
            print(f"  {i:>4}.  RS({s.n1},{s.k1},{s.t1}) ⊕ RS({s.n2},{s.k2},{s.t2})")

    print()


# ─────────────────────────────────────────────────────────────────────────────
# Tests
# ─────────────────────────────────────────────────────────────────────────────

def run_tests() -> bool:
    passed = failed = 0

    def check(cond, name, detail=""):
        nonlocal passed, failed
        if cond:
            passed += 1
            print(f"  [PASS] ✓  {name}")
        else:
            failed += 1
            print(f"  [FAIL] ✗  {name}" + (f"\n         → {detail}" if detail else ""))

    print("\n" + "=" * 60)
    print("  TESTS")
    print("=" * 60)

    # ── 1. Parent validation ─────────────────────────────────────────────────
    print("\n── Validation ──")

    # Valid parents
    for (n, k, t) in [(3,1,1), (7,3,2), (15,9,3), (255,223,16)]:
        ok, _ = validate_parent(n, k, t)
        check(ok, f"RS({n},{k},{t}) accepted as valid")

    # Invalid parents
    for (n, k, t, why) in [
        (7, 4, 2, "n≠k+2t"),
        (4, 0, 2, "k<1"),
        (2, 0, 1, "n<3"),
        (6, 2, 1, "valid: 6=2+4? No, 2+2=4≠6 → n≠k+2t"),
    ]:
        ok, _ = validate_parent(n, k, t)
        check(not ok, f"RS({n},{k},{t}) rejected ({why})")

    # ── 2. can_split O(1) theorem ────────────────────────────────────────────
    print("\n── can_split theorem ──")

    yes_cases = [(6,2,2),(7,3,2),(8,4,2),(9,3,3),(10,4,3),(15,5,5),(255,223,16)]
    for (n, k, t) in yes_cases:
        ok, _ = can_split(n, k, t)
        check(ok, f"RS({n},{k},{t}) → splittable")

    no_cases = [(3,1,1,"t=1"),(4,2,1,"t=1"),(5,3,1,"t=1"),
                (5,1,2,"k=1"),(7,1,3,"k=1"),(9,1,4,"k=1")]
    for (n, k, t, why) in no_cases:
        ok, _ = can_split(n, k, t)
        check(not ok, f"RS({n},{k},{t}) → not splittable ({why})")

    # ── 3. Split counts ──────────────────────────────────────────────────────
    print("\n── Split counts ──")

    # Minimum decomposable: RS(6,2,2) has exactly 1 split
    splits = find_all_splits(6, 2, 2)
    check(len(splits) == 1, "RS(6,2,2): exactly 1 split",
          f"got {len(splits)}")
    if splits:
        s = splits[0]
        check(s.n1==3 and s.k1==1 and s.t1==1 and
              s.n2==3 and s.k2==1 and s.t2==1,
              "RS(6,2,2) → RS(3,1,1) ⊕ RS(3,1,1)")

    # Non-splittable codes give 0 splits
    for (n, k, t) in [(3,1,1),(4,2,1),(5,1,2),(7,1,3)]:
        splits = find_all_splits(n, k, t)
        check(len(splits) == 0, f"RS({n},{k},{t}): 0 splits",
              f"got {len(splits)}")

    # Known counts
    known = [
        (7, 3, 2, 2),
        (8, 4, 2, 3),
        (9, 3, 3, 4),
        (10, 4, 3, 6),
        (15, 5, 5, 16),
        (15, 9, 3, 16),
    ]
    for (n, k, t, expected) in known:
        splits = find_all_splits(n, k, t)
        check(len(splits) == expected,
              f"RS({n},{k},{t}): {expected} splits",
              f"got {len(splits)}")

    # ── 4. Split correctness ─────────────────────────────────────────────────
    print("\n── Split correctness ──")

    for (n, k, t) in [(7,3,2),(10,4,3),(15,9,3),(20,8,6)]:
        splits = find_all_splits(n, k, t)
        all_ok = True
        for s in splits:
            # n1 + n2 = n
            ok_n = s.n1 + s.n2 == n
            # t1 + t2 = t
            ok_t = s.t1 + s.t2 == t
            # k = k1 + k2
            ok_k = s.k1 + s.k2 == k
            # Both sub-codes valid
            ok_c1 = is_valid_rs(s.n1, s.k1, s.t1)
            ok_c2 = is_valid_rs(s.n2, s.k2, s.t2)
            if not all([ok_n, ok_t, ok_k, ok_c1, ok_c2]):
                all_ok = False
                break
        check(all_ok, f"RS({n},{k},{t}): all {len(splits)} splits are valid")

    # ── 5. Theorem vs brute-force consistency (all n ≤ 25) ──────────────────
    print("\n── Theorem ↔ brute-force consistency (n ≤ 25) ──")

    mismatches = 0
    total = 0
    for n in range(3, 26):
        for t in range(1, n // 2 + 1):
            k = n - 2 * t
            if k < 1:
                continue
            total += 1
            theorem_says, _ = can_split(n, k, t)
            brute_says = len(find_all_splits(n, k, t)) > 0
            if theorem_says != brute_says:
                mismatches += 1
    check(mismatches == 0,
          f"Theorem and brute-force agree on all {total} codes with n ≤ 25")

    # ── 6. Symmetric splits ──────────────────────────────────────────────────
    print("\n── Symmetric splits ──")

    # RS(8,4,2): should contain RS(4,2,1)⊕RS(4,2,1)
    splits = find_all_splits(8, 4, 2)
    has_sym = any(s.n1==4 and s.k1==2 and s.t1==1 and
                  s.n2==4 and s.k2==2 and s.t2==1 for s in splits)
    check(has_sym, "RS(8,4,2) contains symmetric split RS(4,2,1)⊕RS(4,2,1)")

    # RS(6,2,2): only one split = symmetric
    splits = find_all_splits(6, 2, 2)
    check(len(splits) == 1 and splits[0].n1 == splits[0].n2,
          "RS(6,2,2) only split is symmetric")

    # ── 7. ValueError on invalid parent ──────────────────────────────────────
    print("\n── Error handling ──")

    for (n, k, t) in [(7,4,2),(4,0,2)]:
        try:
            find_all_splits(n, k, t)
            check(False, f"find_all_splits({n},{k},{t}) should raise ValueError")
        except ValueError:
            check(True, f"find_all_splits({n},{k},{t}) raises ValueError correctly")

    # ── Summary ──────────────────────────────────────────────────────────────
    total_tests = passed + failed
    print(f"\n{'='*60}")
    print(f"  Total: {total_tests}  |  Passed: {passed}  |  Failed: {failed}")
    print(f"  Pass rate: {100*passed/total_tests:.1f}%")
    print("=" * 60)
    return failed == 0


# ─────────────────────────────────────────────────────────────────────────────
# Demo: show splits for a range of codes
# ─────────────────────────────────────────────────────────────────────────────

def demo() -> None:
    print("\n" + "=" * 74)
    print("  REED-SOLOMON CODE SPLITTER — DEMONSTRATION")
    print("=" * 74)

    examples = [
        # (n,  k,  t,  label)
        (3,  1,  1,  "minimum RS code — not splittable"),
        (5,  1,  2,  "k=1 — not splittable"),
        (6,  2,  2,  "minimum splittable RS code"),
        (7,  3,  2,  "GF(8) classic example"),
        (8,  4,  2,  "includes symmetric split"),
        (9,  3,  3,  "t=3 with small n"),
        (10, 4,  3,  "6 splits"),
        (15, 9,  3,  "GF(16), 16 splits"),
        (15, 5,  5,  "high t, 16 splits"),
        (31, 15, 8,  "NASA IAC14-TENEX code"),
        (255, 223, 16, "CD/DVD RS code"),
    ]

    for n, k, t, label in examples:
        print(f"\n  ── {label} ──")
        splits = find_all_splits(n, k, t)
        # Limit display to first 8 splits for large codes
        display = splits[:8]
        print_splits(n, k, t, display, verbose=True)
        if len(splits) > 8:
            print(f"  ... and {len(splits) - 8} more (showing first 8 only)\n")

    # Summary table
    print("\n" + "=" * 74)
    print("  SUMMARY TABLE — Split counts for all valid RS codes (n ≤ 20)")
    print("=" * 74)
    print(f"\n  {'Code':>14}  {'k≥2':>4}  {'t≥2':>4}  {'Splits':>7}  Smallest split")
    print("  " + "-" * 65)
    for n in range(3, 21):
        for t in range(1, n // 2 + 1):
            k = n - 2 * t
            if k < 1:
                continue
            splits = find_all_splits(n, k, t)
            k2 = "✓" if k >= 2 else "✗"
            t2 = "✓" if t >= 2 else "✗"
            if splits:
                s = splits[0]
                smallest = f"RS({s.n1},{s.k1},{s.t1}) ⊕ RS({s.n2},{s.k2},{s.t2})"
            else:
                smallest = "—"
            print(f"  RS({n:>2},{k:>2},{t:>2})  {k2:>4}  {t2:>4}  {len(splits):>7}  {smallest}")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main():
    args = sys.argv[1:]

    if len(args) == 3:
        # Called from command line: python3 rs_splitter.py n k t
        try:
            n, k, t = int(args[0]), int(args[1]), int(args[2])
        except ValueError:
            print("Usage: python3 rs_splitter.py <n> <k> <t>")
            sys.exit(1)

        ok, reason = validate_parent(n, k, t)
        if not ok:
            print(f"\nERROR: {reason}")
            sys.exit(1)

        splits = find_all_splits(n, k, t)
        print_splits(n, k, t, splits, verbose=True)

    else:
        # No args: run tests then demo
        all_passed = run_tests()
        demo()
        if not all_passed:
            sys.exit(1)


if __name__ == "__main__":
    main()