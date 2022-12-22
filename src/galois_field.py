from src.ring_polynomial import RingPolynomial
from src.gf_polynomial import GF_Polynomial
from src.hyperelliptic import HC


class GaloisField:
    def __init__(self, base, polynomial):
        if base != polynomial.gf:
            raise ValueError(f"{polynomial} must be defined over base field {base}")
        if not polynomial.is_irreducible():
            raise ValueError(f"{polynomial} is not irreducible")

        self._poly = polynomial
        self.base = base
        self.p = base.p
        self.m = polynomial.deg
        self.q = self.p**self.m

    def zero(self):
        return GF_Polynomial(self, [self.base.zero()])

    def one(self):
        return GF_Polynomial(self, [self.base.one()])

    def element(self, value):
        if not isinstance(value, list):
            raise ValueError(f"{self} element must be defined by list object")
        parsed_coeff = self.base._parse_coeff(value)
        return GF_Polynomial(self, parsed_coeff)

    def get_elements(self):
        a = self.element([1, 0])
        result = a**0

        for _ in range(self.q - 1):
            yield result
            result = result * a

    def poly(self, coeff):
        parsed_coeff = self._parse_coeff(coeff)
        return RingPolynomial(self, parsed_coeff)

    def rand_element(self):
        return self.element([self.base.rand_element() for _ in range(self.m)])

    def hyperelliptic(self, h, f):
        return HC(self, h, f)

    def _is_field_element(self, value):
        return isinstance(value, GF_Polynomial) and value.gf == self

    def _parse_coeff(self, coeff):
        return list(
            map(lambda x: x if self._is_field_element(x) else self.element(x), coeff)
        )

    def __call__(self, value):
        if isinstance(value, int) or self.base._is_field_element(value):
            return self.element([value])
        if isinstance(value, list):
            return self.element(value)
        raise ValueError(f"Invalid value for {self} element")

    def __eq__(self, other):
        return self.base == other.base and self._poly == other._poly

    def __str__(self):
        return f"Galois Field mod {self.p} mod {self._poly}"
