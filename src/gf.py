from src.polynomial import Polynomial
from src.gf_polynomial import GF_Polynomial
from src.hyperelliptic import HC
from src.integer import ZP
from src.utils import is_prime, factors
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
        if not is_prime(p):
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

    def sqrt(self, a):
        if isinstance(a, int):
            a = self.int(a)
        if isinstance(a, ZP):
            return a.sqrt()
        raise ValueError(f"{a} does not have sqrt method")

    def factors(self, n):
        return factors(n)

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
