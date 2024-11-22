from random import randint

from .utils import is_prime
from .galois_field import GaloisField
from .ring_polynomial import RingPolynomial
from .hyperelliptic import HC
from .integer import ZP


class FiniteField:
    """Provides set of tools which allow operations over Finite Field"""

    def __init__(self, p: int) -> None:
        if not is_prime(p):
            raise ValueError(f"{p} is not prime")

        self.p: int = p

    def zero(self) -> ZP:
        """Return addition neutral element of Finite Field"""
        return ZP(self, 0)

    def one(self) -> ZP:
        """Return multiplication neutral element of Finite Field"""
        return ZP(self, 1)

    def element(self, value: int | ZP) -> ZP:
        """Return element of Finite Field (integer) from scalar value"""
        return ZP(self, value)

    def get_elements(self) -> range:
        """Returns range of all elements in Finite Field"""
        return range(0, self.p)

    def poly(self, coeff: list[ZP] | list[int], symbol: str = "x") -> RingPolynomial:
        """Returns polynomial over Finite Field"""
        parsed_coeff = self._parse_coeff(coeff)
        return RingPolynomial(self, parsed_coeff, symbol)

    def extension(self, polynomial: RingPolynomial) -> GaloisField:
        """Return linear space (Galois Field) over Finite Field using irreducible polynomial provided in argument"""
        return GaloisField(self, polynomial)

    def rand_element(self):
        """Returns random element from Finite Field"""
        return self.element(randint(0, self.p - 1))

    def rand_poly(self, deg: int):
        """Returns random vector (polynomial) defined over Finite Field"""
        return self.poly([self.rand_element() for _ in range(deg + 1)])

    def rand_irreducible_poly(self, deg: int) -> RingPolynomial:
        """Returns random irreducible polynomial defined over Finite Field"""
        leading_coeff = self.one()
        poly = self.poly([leading_coeff] + [self.rand_element() for _ in range(deg)])
        while not poly.is_irreducible():
            poly = self.poly(
                [leading_coeff] + [self.rand_element() for _ in range(deg)]
            )
        return poly

    def hyperelliptic(self, h: RingPolynomial, f: RingPolynomial) -> HC:
        """Returns hyperelliptic curve defined over Finite Field"""
        return HC(self, h, f)

    def _is_field_element(self, value) -> bool:
        return isinstance(value, ZP) and value.gf == self

    def _parse_coeff(self, coeff) -> list[ZP]:
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
