from src.polynomial import Polynomial
from src.gf_polynomial import GF_Polynomial
from src.hyperelliptic import HC
from src.integer import ZP
from src.utils import isPrime
from random import randint


class GaloisField:
    def __init__(self, p, m=1, polynomial_coeff=None):
        # TODO:
        # Check i polynomial is set and if it's irreducible
        if m < 1:
            raise ValueError(f"Invalid exponent parameter for GF: {m}")
        if m > 1 and polynomial_coeff == None:
            raise ValueError(
                "Field with q != p must have defined irreducible polynomial"
            )
        if not isPrime(p):
            raise ValueError(f"{p} is not prime")

        self.p = p
        self.m = m
        self.q = p**m
        self._poly = (
            Polynomial(self._parse_coeff(polynomial_coeff), self)
            if polynomial_coeff
            else None
        )

    @property
    def zero(self):
        if self.m == 1:
            return self.int_zero
        return GF_Polynomial(self, [self.int_zero])

    @property
    def one(self):
        if self.m == 1:
            return self.int_one
        return GF_Polynomial(self, [self.int_one])

    @property
    def int_zero(self):
        return ZP(self, 0)

    @property
    def int_one(self):
        return ZP(self, 1)

    @property
    def poly_zero(self):
        return Polynomial([self.zero], self)

    @property
    def poly_one(self):
        return Polynomial([self.one], self)

    def int(self, value):
        return ZP(self, value)

    def element(self, value):
        if self.m > 1:
            if not isinstance(value, list):
                raise ValueError(
                    "Element of field is not defined as array of coefficients"
                )
            return GF_Polynomial(self, self._parse_coeff(value))
        return ZP(self, value)

    def rand_int(self):
        return self.int(randint(0, self.p - 1))

    def rand_poly(self, deg):
        return self.poly([self.rand_int() for _ in range(deg + 1)])

    def poly(self, coeff, symbol="x"):
        parsed_coeff = self._parse_coeff(coeff)
        return Polynomial(parsed_coeff, self, symbol)

    def hyperelliptic(self, h, f):
        return HC(self, h, f)

    # NOTE: Solution can be lifted vai hensell lemma from p to p^m
    def sqrt(self, a):
        if not self.is_quadratic_residue(a):
            raise ValueError(f"Argument {a} has no square root")

        if a == 0:
            return self.int(0)

        a = a.value if isinstance(a, ZP) else a
        p = self.p
        q = p - 1
        s = 0
        while q & 1 == 0:
            q //= 2
            s += 1
        if s == 1:
            return pow(a, (p + 1) // 4, p)
        for test in range(2, p):
            if pow(test, (p - 1) // 2, p) == p - 1:
                break

        c = pow(test, q, p)
        r = pow(a, (q + 1) // 2, p)
        t = pow(a, q, p)
        m = s
        t2 = 0
        while (t - 1) % p != 0:
            t2 = t**2 % p
            for i in range(1, m):
                if (t2 - 1) % p == 0:
                    break
                t2 = t2**2 % p
            b = pow(c, 1 << (m - i - 1), p)
            r = (r * b) % p
            c = (b * b) % p
            t = (t * c) % p
            m = i
        return self.int(r)

    def is_quadratic_residue(self, a):
        if self.legendre(a) == 1 or a == 0:
            return True
        return False

    def legendre(self, a):
        p = self.p
        return a ** ((p - 1) // 2)

    def mod_inv(self, a):
        if a == 0:
            return 0

        p = self.p
        t0, t1 = 0, 1
        r0, r1 = p, a

        while r1 != 0:
            q = r0 // r1
            t0, t1 = t1, t0 - q * t1
            r0, r1 = r1, r0 - q * r1

        if r0 > 1:
            raise Exception(f"Element {a} is not inversable mod {p}")

        if t0 < 0:
            t0 += p

        return t0

    def get_elements(self, limit=0):
        if self.m == 1:
            for _ in range(self.p if limit == 0 else limit):
                yield self.int(_)
            return

        if limit <= 0:
            limit = self.q

        a = self.element([1, 0])
        result = a**0

        for _ in range(min(limit, self.q)):
            yield result
            result = result * a

    def __call__(self, arg):
        if isinstance(arg, int):
            return self.int(arg)
        elif isinstance(arg, list):
            return self.poly(arg)
        return

    def __eq__(self, other):
        if not isinstance(other, GaloisField):
            return False
        return self.q == other.q and self._poly == other._poly

    def _parse_coeff(self, coeff):
        return list(map(lambda x: ZP(self, x) if isinstance(x, int) else x, coeff))
