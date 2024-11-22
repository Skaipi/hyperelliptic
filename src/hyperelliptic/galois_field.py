from .ring_polynomial import RingPolynomial
from .gf_polynomial import GF_Polynomial
from .hyperelliptic import HC
from .integer import ZP


class GaloisField:
    """Provides set of tools which allow operations over Galois Field"""

    def __init__(self, base: "FiniteField", polynomial: RingPolynomial):
        if base != polynomial.gf:
            raise ValueError(f"{polynomial} must be defined over base field {base}")
        if not polynomial.is_irreducible():
            raise ValueError(f"{polynomial} is not irreducible")

        self._poly: RingPolynomial = polynomial
        self.base: "FiniteField" = base
        self.p: int = base.p
        self.m: int = polynomial.deg
        self.q: int = self.p**self.m

    def zero(self) -> GF_Polynomial:
        """Return addition neutral element of Galois Field"""
        return GF_Polynomial(self, [self.base.zero()])

    def one(self) -> GF_Polynomial:
        """Return multiplication neutral element of Galois Field"""
        return GF_Polynomial(self, [self.base.one()])

    def element(self, value: list[ZP | int] | ZP | int) -> GF_Polynomial:
        """Return element of Galois Field (polynomial) from list of coefficients or scalar value"""
        if isinstance(value, int) or isinstance(value, ZP):
            value = [value]
        if not isinstance(value, list):
            raise ValueError(f"{self} element must be defined by list object")
        parsed_coeff = self.base._parse_coeff(value)  # pylint: disable=protected-access
        return GF_Polynomial(self, parsed_coeff)

    def get_elements(self):
        """Generator function for iterating over all elements in Galois Field"""
        a = self.element([1, 0])
        result = a**0

        for _ in range(self.q - 1):
            yield result
            result = result * a

    def poly(self, coeff: list[ZP] | list[int]) -> RingPolynomial:
        """Returns polynomial over Galois Field. A polynomial whose coefficients are GF elements."""
        parsed_coeff = self._parse_coeff(coeff)
        return RingPolynomial(self, parsed_coeff)

    def rand_element(self):
        """Returns random element from Galois Field"""
        return self.element([self.base.rand_element() for _ in range(self.m)])

    def hyperelliptic(self, h: RingPolynomial, f: RingPolynomial) -> HC:
        """Returns hyperelliptic curve defined over Galois Field"""
        return HC(self, h, f)

    def _is_field_element(self, value: object) -> bool:
        return isinstance(value, GF_Polynomial) and value.gf == self

    def _parse_coeff(self, coeff: list[ZP] | list[int]) -> list[GF_Polynomial]:
        return list(map(self.element, coeff))

    def __call__(self, value) -> GF_Polynomial:
        if isinstance(value, int) or self.base._is_field_element(value):
            return self.element([value])
        if isinstance(value, list):
            return self.element(value)
        raise ValueError(f"Invalid value for {self} element")

    def __eq__(self, other: object):
        return (
            isinstance(other, GaloisField)
            and self.base == other.base
            and self._poly == other._poly
        )

    def __str__(self):
        return f"Galois Field mod {self.p} mod {self._poly}"
