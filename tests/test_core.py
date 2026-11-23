"""Tests for tonelli_shanks_algorithm_utils.core."""

import unittest

from tonelli_shanks_algorithm_utils import (
    is_quadratic_residue,
    mod_sqrt,
    tonelli_shanks,
)


class TestTonelliShanks(unittest.TestCase):
    def test_small_primes(self):
        # p = 3, 5, 7, 11, 13, 17, 19
        cases = [
            (1, 3, 1),
            (2, 7, 4),  # 4^2 = 16 = 2 mod 7
            (4, 5, 2),  # 2^2 = 4 mod 5
            (9, 13, 3),  # 3^2 = 9 mod 13
            (10, 13, 6),  # 6^2 = 36 = 10 mod 13
            (2, 17, 6),  # 6^2 = 36 = 2 mod 17
            (4, 17, 2),  # 2^2 = 4 mod 17
            (9, 19, 3),  # 3^2 = 9 mod 19
            (0, 7, 0),
        ]
        for n, p, expected in cases:
            with self.subTest(n=n, p=p):
                root = tonelli_shanks(n, p)
                self.assertIsNotNone(root)
                self.assertEqual((root * root) % p, n % p)
                # The returned root must be one of the two roots.  We accept
                # either; the expected value is just a sanity check.
                if expected is not None:
                    self.assertIn(root, (expected, (-expected) % p))

    def test_non_residues_return_none(self):
        # 2 is a non-residue mod 3, 3 is a non-residue mod 5, etc.
        # Note: 3 is a quadratic residue mod 13 (4^2 = 16 = 3), so it is
        # not included here.
        for n, p in [(2, 3), (3, 5), (5, 7), (7, 11), (2, 11), (5, 13)]:
            with self.subTest(n=n, p=p):
                self.assertIsNone(tonelli_shanks(n, p))

    def test_p_mod_4_eq_3_shortcut(self):
        # p = 7, 11, 19 are all 3 mod 4.
        for p in (7, 11, 19):
            for n in range(1, p):
                root = tonelli_shanks(n, p)
                if is_quadratic_residue(n, p):
                    self.assertIsNotNone(root)
                    self.assertEqual((root * root) % p, n)
                else:
                    self.assertIsNone(root)

    def test_p_mod_4_eq_1_general_case(self):
        # p = 5, 13, 17, 29 are all 1 mod 4.
        for p in (5, 13, 17, 29):
            for n in range(1, p):
                root = tonelli_shanks(n, p)
                if is_quadratic_residue(n, p):
                    self.assertIsNotNone(root)
                    self.assertEqual((root * root) % p, n)
                else:
                    self.assertIsNone(root)

    def test_large_prime(self):
        p = 1000003  # prime
        n = 123456
        root = tonelli_shanks(n, p)
        if root is None:
            self.assertFalse(is_quadratic_residue(n, p))
        else:
            self.assertEqual((root * root) % p, n)

    def test_negative_n(self):
        # Negative n should be reduced modulo p first.
        p = 13
        n = -1  # -1 mod 13 = 12, which is a quadratic residue (5^2 = 25 = 12)
        root = tonelli_shanks(n, p)
        self.assertIsNotNone(root)
        self.assertEqual((root * root) % p, n % p)  # 12
        # Also test a negative non-residue.
        n = -2  # -2 mod 13 = 11, which is a non-residue
        self.assertIsNone(tonelli_shanks(n, p))

    def test_zero(self):
        for p in (3, 5, 7, 11, 13):
            self.assertEqual(tonelli_shanks(0, p), 0)

    def test_invalid_p(self):
        with self.assertRaises(ValueError):
            tonelli_shanks(1, 2)
        with self.assertRaises(ValueError):
            tonelli_shanks(1, 1)
        with self.assertRaises(ValueError):
            tonelli_shanks(1, 0)
        with self.assertRaises(ValueError):
            tonelli_shanks(1, -3)

    def test_is_quadratic_residue(self):
        self.assertTrue(is_quadratic_residue(1, 7))
        self.assertTrue(is_quadratic_residue(2, 7))  # 3^2=2 mod 7
        self.assertTrue(is_quadratic_residue(4, 7))
        self.assertFalse(is_quadratic_residue(3, 7))
        self.assertFalse(is_quadratic_residue(5, 7))
        self.assertFalse(is_quadratic_residue(6, 7))
        self.assertTrue(is_quadratic_residue(0, 7))

    def test_mod_sqrt_alias(self):
        self.assertIs(mod_sqrt, tonelli_shanks)
        self.assertEqual(mod_sqrt(4, 7), tonelli_shanks(4, 7))


if __name__ == "__main__":
    unittest.main()
