from .polynomial import Polynomial


class GF_Polynomial(Polynomial):
    def __init__(self, field, coeff, symbol="a"):
        self.gf = field
        super().__init__(coeff, symbol)

        if self.deg >= field._poly.deg:
            self.coeff = (self % field._poly).coeff

    def coeff_zero(self):
        return self.gf.base.zero()

    def coeff_one(self):
        return self.gf.base.one()

    def inverse(self):
        if self == self.zero():
            raise ZeroDivisionError("Element 0 has no inverse")

        d, a, b = self.xgcd(self.gf._poly)
        if not d.isConst():
            raise Exception("Element has no inverse!")
        return a

    def sqrt(self):
        # Handle special case of char(2) fields
        if self.gf.p == 2:
            return pow(self, self.gf.q // 2)

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

        while t != 1:
            # Find the least i such that t^2^i = 1
            i = 1
            while pow(t, pow(2, i)) != 1:
                i += 1

            b = pow(c, pow(2, m - i - 1))
            m = i
            c = pow(b, 2)
            t = t * c
            r = r * b
        return r

    def is_quadratic_residue(self):
        return self.legendre() == self.one() or self == self.zero()

    def legendre(self):
        return pow(self, ((self.gf.q - 1) // 2))

    def _from_coeff(self, coeff):
        return GF_Polynomial(self.gf, coeff, self.symbol)

    def __truediv__(self, other):
        if isinstance(other, GF_Polynomial):
            return self * other.inverse()
        return super().__truediv__(other)

    def __hash__(self):
        result = 0
        for i, c in enumerate(self.coeff[::-1]):
            result += c.value * self.gf.p**i
        return result
