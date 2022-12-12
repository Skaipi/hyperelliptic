from src.field_polynomial import FieldPolynomial
from src.utils import gf_operation
from src.integer import ZP


class GF_Polynomial:
    def __init__(self, gf, coeff, symbol="a") -> None:
        self._poly = FieldPolynomial(coeff, gf, symbol) % gf._poly
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

    def to_int(self):
        result = 0
        for i, c in enumerate(self.coeff[::-1]):
            result += c * self.gf.p**i
        return result

    def inverse(self):
        if self == self.zero():
            raise ZeroDivisionError("Element 0 has no inverse")

        d, a, b = self.xgcd(self.gf)
        if not d.isConst():
            raise Exception("Element has no inverse!")

        return self._from_zp_polynomial(a / d.toInt())

    def sqrt(self):
        # Tonelli-Shanks algorithm
        if not self.is_quadratic_residue():
            raise ValueError(f"Argument {self} has no square root")

        if self == self.zero():
            return self.zero()

        p = self.gf.q
        q = p - 1
        s = 0  # S must be an integer, not ZP

        # find Q and S such that p - 1 = Q * 2^S
        while q % 2 == 0:
            q //= 2
            s += 1

        if s == 1:
            return pow(self, (p + 1) // 4)

        # Find z that is not quadratic residue
        z = self.gf.rand_element()
        while z.is_quadratic_residue():
            z = self.gf.rand_element()

        c = pow(z, q)
        t = pow(self, q)
        r = pow(self, (q + 1) // 2)
        m = s

        while True:
            if t == 1:
                return r

            # Find the least i such that t^2^i = 1
            i = 1
            while pow(t, pow(2, i)) != 1:
                i += 1

            b = pow(c, pow(2, m - i - 1))
            m = i
            c = pow(b, 2)
            t = t * c
            r = r * b

    def is_quadratic_residue(self):
        if self.legendre() == self.one() or self == self.zero():
            return True
        return False

    def legendre(self):
        return pow(self, ((self.gf.q - 1) // 2))

    def _from_zp_polynomial(self, poly):
        return GF_Polynomial(self.gf, poly.coeff, poly.symbol)

    @gf_operation
    def __add__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._from_zp_polynomial(self._poly + other)
        result = self._poly + other._poly
        return self._from_zp_polynomial(result)

    @gf_operation
    def __radd__(self, other):
        return self.__add__(other)

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
    def __rmul__(self, other):
        return self.__mul__(other)

    @gf_operation
    def __pow__(self, other):
        result = pow(self._poly, other, self.gf._poly)
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

    @gf_operation
    def __mod__(self, other):
        return divmod(self, other)[1]

    @gf_operation
    def __floordiv__(self, other):
        return divmod(self, other)[0]

    @gf_operation
    def __eq__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            return self._poly == other
        if isinstance(other, GF_Polynomial):
            return self._poly == other._poly
        return False

    def __neg__(self):
        return self._from_zp_polynomial(-self._poly)

    def __call__(self, x):
        result = self.gf.int_zero
        for i, c in enumerate(self.coeff[::-1]):
            result += x**i * c
        return result

    def __str__(self):
        return str(self._poly)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return self.to_int().value
