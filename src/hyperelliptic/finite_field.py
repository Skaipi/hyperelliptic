from random import randint

from .utils import is_prime
from .galois_field import GaloisField
from .ring_polynomial import RingPolynomial
from .hyperelliptic import HC
from .integer import ZP


class FiniteField:
    def __init__(self, p):
        if not is_prime(p):
            raise ValueError(f"{p} is not prime")

        self.p = p

    def zero(self):
        return ZP(self, 0)

    def one(self):
        return ZP(self, 1)

    def element(self, value):
        return ZP(self, value)

    def get_elements(self):
        return range(0, self.p)

    def poly(self, coeff, symbol="x"):
        parsed_coeff = self._parse_coeff(coeff)
        return RingPolynomial(self, parsed_coeff, symbol)

    def extension(self, polynomial):
        return GaloisField(self, polynomial)

    def rand_element(self):
        return self.element(randint(0, self.p - 1))

    def rand_poly(self, deg):
        return self.poly([self.rand_element() for _ in range(deg + 1)])

    def rand_irreducible_poly(self, deg):
        leading_coeff = self.one()
        poly = self.poly([leading_coeff] + [self.rand_element() for _ in range(deg)])
        while not poly.is_irreducible():
            poly = self.poly(
                [leading_coeff] + [self.rand_element() for _ in range(deg)]
            )
        return poly

    def hyperelliptic(self, h, f):
        return HC(self, h, f)

    def _is_field_element(self, value):
        return isinstance(value, ZP) and value.gf == self

    def _parse_coeff(self, coeff):
        return list(
            map(lambda x: x if self._is_field_element(x) else self.element(x), coeff)
        )

    def __call__(self, value):
        if isinstance(value, int):
            return self.element(value)
        raise ValueError(f"Invalid value for {self} element")

    def __eq__(self, other):
        return self.p == other.p

    def __str__(self):
        return f"Finite Field mod {self.p}"
