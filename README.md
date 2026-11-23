# tonelli-shanks-algorithm-utils

A small, dependency-free Python library that computes square roots modulo an odd prime using the Tonelli-Shanks algorithm.

```python
from tonelli_shanks_algorithm_utils import tonelli_shanks, is_quadratic_residue

p = 13
n = 10
root = tonelli_shanks(n, p)
print(root)          # 6 or 7
print((root * root) % p == n)  # True

print(is_quadratic_residue(3, 13))  # False
print(tonelli_shanks(3, 13))        # None
```

## Why this library exists

Computing square roots modulo a prime is a common building block in number theory, cryptography, and competitive programming.  Python's standard library provides `pow` for modular exponentiation but no modular square root.  This library fills that gap with a single, focused implementation of the classic Tonelli-Shanks algorithm.

The trade-off is explicit: the library does **not** verify that `p` is prime.  A deterministic primality check would add significant code and runtime cost, and the algorithm is only valid for prime moduli anyway.  The caller is expected to supply a prime.  For the same reason, the function returns `None` rather than raising an exception when `n` is not a quadratic residue; this keeps the control flow simple for callers who are already testing residues with `is_quadratic_residue`.

## Edge cases

- `n` is reduced modulo `p` first, so negative inputs are accepted.
- `n == 0` returns `0`.
- `p == 2` is rejected with `ValueError` because the algorithm as implemented requires an odd prime.  (In GF(2) the only square root of 1 is 1, and of 0 is 0, but handling that separately would complicate the API for no practical gain.)
- The function returns one of the two square roots; which one is unspecified.

## Exports

- `tonelli_shanks(n: int, p: int) -> Optional[int]`
- `mod_sqrt(n: int, p: int) -> Optional[int]` (alias for `tonelli_shanks`)
- `is_quadratic_residue(n: int, p: int) -> bool`
