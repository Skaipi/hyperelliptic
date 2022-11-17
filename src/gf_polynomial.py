from src.polynomial import Polynomial
from src.utils import gf_operation
from src.integer import ZP


class GF_Polynomial:
    def __init__(self, gf, coeff, symbol="a") -> None:
        self._poly = Polynomial(coeff, gf, symbol) % gf._poly
        self.coeff = self._poly.coeff
        self.symbol = symbol
        self.gf = gf

    @property
    def deg(self) -> int:
        return len(self.coeff) - 1

    def zero(self):
        coeff = [self.coeff[0].zero()]
        return GF_Polynomial(self.gf, coeff, self.symbol)

    def one(self):
        coeff = [self.coeff[0].one()]
        return GF_Polynomial(self.gf, coeff, self.symbol)

    def xgcd(self, other):
        return self._poly.xgcd(other._poly)

    def to_monic(self):
        return self._from_zp_polynomial(self._poly.to_monic())

    def inverse(self):
        if self == self.zero():
            raise ZeroDivisionError("Element 0 has no inverse")

        d, a, b = self.xgcd(self.gf)
        if not d.isConst():
            raise Exception("Element has no inverse!")
        return self._from_zp_polynomial(a.to_monic())

    def _from_zp_polynomial(self, poly):
        return GF_Polynomial(self.gf, poly.coeff, poly.symbol)

    @gf_operation
    def __add__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._from_zp_polynomial(self._poly + other)
        result = self._poly + other._poly
        return self._from_zp_polynomial(result)

    @gf_operation
    def __sub__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._from_zp_polynomial(self._poly - other)
        result = self._poly - other._poly
        return self._from_zp_polynomial(result)

    @gf_operation
    def __mul__(self, other):
        if isinstance(other, int) or isinstance(other, ZP):
            return self._from_zp_polynomial(self._poly * other)
        result = self._poly * other._poly
        return self._from_zp_polynomial(result)

    @gf_operation
    def __pow__(self, other):
        result = self._poly**other
        return self._from_zp_polynomial(result)

    @gf_operation
    def __divmod__(self, other):
        r1, r2 = divmod(self._poly, other._poly)
        return self._from_zp_polynomial(r1), self._from_zp_polynomial(r2)

    @gf_operation
    def __truediv__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._from_zp_polynomial(self._poly / other)
        return self._from_zp_polynomial(self._poly * other.inverse())

    def __neg__(self):
        return self._from_zp_polynomial(-self._poly)

    def __call__(self, x):
        result = self.gf.int_zero
        for i, c in enumerate(self.coeff[::-1]):
            result += x**i * c
        return result

    @gf_operation
    def __mod__(self, other):
        return divmod(self, other)[1]

    @gf_operation
    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __str__(self):
        return str(self._poly)

    def __repr__(self):
        return str(self)

    @gf_operation
    def __eq__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._poly == other
        if isinstance(other, GF_Polynomial):
            return self._poly == other._poly
        return False
