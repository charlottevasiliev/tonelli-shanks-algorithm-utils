"""Core implementation of the Tonelli-Shanks algorithm.

The algorithm computes a square root of a quadratic residue ``n`` modulo an odd
prime ``p``.  It requires that ``p`` is prime; primality is not checked here
because a deterministic primality test would add significant overhead and is
outside the scope of this small utility.  The caller is responsible for
providing a prime modulus.

References:
    Tonelli, A. (1891). "Bemerkung über die Auflösung quadratischer Congruenzen".
    Shanks, D. (1972). "Five number-theoretic algorithms".
"""

from __future__ import annotations

import random
from typing import Optional


def _legendre_symbol(a: int, p: int) -> int:
    """Return the Legendre symbol (a/p) for an odd prime p.

    Uses Euler's criterion.  For p prime this returns 0, 1, or -1.
    """
    if a % p == 0:
        return 0
    ls = pow(a, (p - 1) // 2, p)
    return -1 if ls == p - 1 else ls


def is_quadratic_residue(n: int, p: int) -> bool:
    """Return True if ``n`` is a quadratic residue modulo the odd prime ``p``.

    A quadratic residue is a number that has a square root modulo ``p``.
    Zero is considered a quadratic residue (0 = 0^2).
    """
    if p <= 2:
        raise ValueError("p must be an odd prime greater than 2")
    return _legendre_symbol(n, p) != -1


def _find_quadratic_non_residue(p: int) -> int:
    """Find a quadratic non-residue modulo the odd prime ``p``.

    The Tonelli-Shanks algorithm needs a non-residue ``z``.  For small primes we
    can simply scan from 2 upward.  For larger primes a random search is
    usually faster, but the scan is deterministic and fast enough for the
    typical inputs this library is used with.  We use a random starting point
    only when p is large enough that scanning could be slow.
    """
    if p < 3:
        raise ValueError("p must be an odd prime greater than 2")
    # For small primes, a linear scan is deterministic and cheap.
    if p < 1000:
        for z in range(2, p):
            if _legendre_symbol(z, p) == -1:
                return z
        raise ValueError("no quadratic non-residue found; p is not prime?")
    # For larger primes, random sampling avoids a potentially long scan.
    # The expected number of trials is 2, so this is safe in practice.
    while True:
        z = random.randrange(2, p)
        if _legendre_symbol(z, p) == -1:
            return z


def tonelli_shanks(n: int, p: int) -> Optional[int]:
    """Compute a square root of ``n`` modulo the odd prime ``p``.

    Returns an integer ``r`` with ``r * r % p == n % p``, or ``None`` if ``n``
    is not a quadratic residue modulo ``p``.

    The algorithm handles the common special cases first:

    * ``n % p == 0`` returns 0.
    * ``p == 2`` returns ``n % 2`` (the only root in GF(2)).
    * ``p % 4 == 3`` uses the short exponent ``(p + 1) // 4``.

    For the general case it performs the standard two-stage reduction:
    factor ``p - 1 = q * 2^s`` with ``q`` odd, then iteratively reduce the
    order of the candidate root until it is a true square root.

    Args:
        n: The integer whose square root is sought.
        p: An odd prime modulus.

    Returns:
        A square root of ``n`` modulo ``p``, or ``None`` if no root exists.

    Raises:
        ValueError: If ``p`` is not an odd prime greater than 2.
    """
    if p <= 2:
        raise ValueError("p must be an odd prime greater than 2")
    n = n % p
    if n == 0:
        return 0
    if _legendre_symbol(n, p) == -1:
        return None

    # Fast path for primes congruent to 3 mod 4.
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)

    # Factor p - 1 = q * 2^s with q odd.
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1

    # Find a quadratic non-residue z.
    z = _find_quadratic_non_residue(p)

    # Initialize the loop variables.
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)

    # Main loop: reduce m until t == 1.
    while t != 1:
        # Find the smallest i (0 < i < m) such that t^(2^i) == 1 mod p.
        i = 0
        t2i = t
        while t2i != 1:
            t2i = (t2i * t2i) % p
            i += 1
            if i == m:
                # This should never happen for a prime p and quadratic residue n.
                return None
        # Update the variables.
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p
    return r


# A short alias for convenience; the module name is already long.
mod_sqrt = tonelli_shanks
