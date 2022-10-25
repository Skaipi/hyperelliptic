from src.polynomial import Polynomial
from src.integer import ZP
from src.utils import xgcd


class GF_Polynomial:
    def __init__(self, gf, coeff, symbol="a") -> None:
        self._poly = Polynomial(coeff, gf, symbol) % gf._poly
        self.coeff = self._poly.coeff
        self.symbol = symbol
        self.gf = gf

    def _from_zp_polynomial(self, poly):
        return GF_Polynomial(self.gf, poly.coeff, poly.symbol)

    def zero(self):
        coeff = [self.coeff[0].zero()]
        return GF_Polynomial(self.gf, coeff, self.symbol)

    def one(self):
        coeff = [self.coeff[0].one()]
        return GF_Polynomial(self.gf, coeff, self.symbol)

    @property
    def deg(self) -> int:
        return len(self.coeff) - 1

    def inverse(self):
        if self == self.zero():
            raise ZeroDivisionError("Element 0 has no inverse")

        d, a, b = xgcd(self._poly, self.gf._poly)
        if not d.isConst():
            raise Exception("Element has no inverse!")
        d_inv = d.toInt().inverse()
        result = a * d_inv % self.gf._poly
        return self._from_zp_polynomial(result)

    def __add__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._from_zp_polynomial(self._poly + other)
        result = self._poly + other._poly
        return self._from_zp_polynomial(result)

    def __sub__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._from_zp_polynomial(self._poly - other)
        result = self._poly - other._poly
        return self._from_zp_polynomial(result)

    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, ZP):
            return self._poly * other
        result = self._poly * other._poly
        return self._from_zp_polynomial(result)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, other):
        result = self._poly**other
        return self._from_zp_polynomial(result)

    def __divmod__(self, other):
        r1, r2 = divmod(self._poly, other._poly)
        return self._from_zp_polynomial(r1), self._from_zp_polynomial(r2)

    def __div__(self, other):
        return self._from_zp_polynomial(self._poly / other)

    def __truediv__(self, other):
        return self._from_zp_polynomial(self._poly / other)

    def __neg__(self):
        return self._from_zp_polynomial(-self._poly)

    def __call__(self, x):
        result = self.gf.int_zero
        for i, c in enumerate(self.coeff[::-1]):
            result += x**i * c
        return result

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __str__(self):
        return str(self._poly)

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._poly == other
        return self._poly == other._poly
