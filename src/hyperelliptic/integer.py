"""(module) containing ZP class"""

import math

from .utils import gf_operation


class ZP:
    """Class representing integers over finite field"""

    def __init__(self, gf, value):
        self.gf = gf
        self.p = gf.p
        self.value = value % self.p

    def zero(self):
        """Neutral element of addition"""
        return self._from_value(0)

    def one(self):
        """Neutral element of multiplication"""
        return self._from_value(1)

    @gf_operation
    def gcd(self, other):
        """Euclidean algorithm"""
        return math.gcd(self.value, other.value)

    @gf_operation
    def xgcd(self, other):
        """Extended Euclidean algorithm"""
        r0 = other.value if isinstance(other, ZP) else other
        r1 = self.value
        s1, s0 = 1, 0
        t1, t0 = 0, 1

        while r0 != 0:
            q = r1 // r0
            r1, r0 = r0, r1 - q * r0
            s1, s0 = s0, s1 - q * s0
            t1, t0 = t0, t1 - q * t0

        return self._from_value(r1), self._from_value(s1), self._from_value(t1)

    def sqrt(self):
        """Tonelli-Shanks algorithm"""
        if not self.is_quadratic_residue():
            raise ValueError(f"Argument {self} has no square root")

        if self == self.zero():
            return self.zero()

        p = self.p
        q = p - 1
        s = 0  # S must be an integer, not ZP

        # find Q and S such that p - 1 = Q * 2^S
        while q % 2 == 0:
            q //= 2
            s += 1

        if s == 1:
            return pow(self, (p + 1) // 4)

        # Find z that is not quadratic residue
        z = self.one()
        while z.is_quadratic_residue():
            z += 1

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
        """Returns true if square root of element exists"""
        if self.legendre() == self.one() or self == self.zero():
            return True
        return False

    def legendre(self):
        """Legendre symbol of an element"""
        return pow(self, ((self.p - 1) // 2))

    def inverse(self):
        """Find inverse y of an element x such that x^{-1} = y and xy = 1"""
        if self == 0:
            raise ZeroDivisionError("Element 0 is not inversable")

        r, s, t = self.xgcd(self.p)
        if r > 1:
            raise TypeError(f"Element {self.value} is not inversable mod {self.p}")
        return s

    def _from_value(self, value):
        return ZP(self.gf, value)

    @gf_operation
    def __add__(self, other):
        if isinstance(other, int):
            return self + self._from_value(other)
        if not isinstance(other, ZP):  # Avoid circular dependency with Polynomials
            return other + self
        return self._from_value(self.value + other.value)

    @gf_operation
    def __radd__(self, other):
        return self.__add__(other)

    @gf_operation
    def __sub__(self, other):
        if isinstance(other, int):
            return self - self._from_value(other)
        if not isinstance(other, ZP):  # Avoid circular dependency with Polynomials
            return other - self
        return self._from_value(self.value - other.value)

    @gf_operation
    def __rsub__(self, other):
        return self.__sub__(other)

    @gf_operation
    def __mul__(self, other):
        if isinstance(other, int):
            return self * self._from_value(other)
        if not isinstance(other, ZP):  # Avoid circular dependency with Polynomials
            return other * self
        return self._from_value(self.value * other.value)

    @gf_operation
    def __rmul__(self, other):
        return self.__mul__(other)

    @gf_operation
    def __truediv__(self, other):
        if isinstance(other, int):
            return self._from_value(
                self.value * self._from_value(other).inverse().value
            )
        return self._from_value(self.value * other.inverse().value)

    @gf_operation
    def __rtruediv__(self, other):
        return self.inverse() * other

    @gf_operation
    def __divmod__(self, other):
        if isinstance(other, int):
            q, r = divmod(self.value, other)
            return self._from_value(q), self._from_value(r)
        q, r = divmod(self.value, other.value)
        return self._from_value(q), self._from_value(r)

    @gf_operation
    def __mod__(self, other):
        return divmod(self, other)[1]

    @gf_operation
    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __pow__(self, other, mod=None):
        if isinstance(other, ZP):
            return self._pow(other.value, mod)
        if isinstance(other, int):
            return self._pow(other, mod)
        raise ValueError("Exponent must be an integer")

    def __neg__(self):
        return self._from_value(-self.value)

    def _pow(self, other, mod=None):
        if mod is not None:
            return self._from_value(pow(self.value, other, mod))
        return self._from_value(pow(self.value, other, self.p))

    @gf_operation
    def __gt__(self, other):
        if isinstance(other, int):
            return self.value > (other % self.p)
        return self.value > other.value

    @gf_operation
    def __lt__(self, other):
        if isinstance(other, int):
            return self.value < (other % self.p)
        return self.value < other.value

    def __eq__(self, other):
        if isinstance(other, int):
            return self.value == (other % self.p)
        if isinstance(other, ZP):
            return self.value == other.value
        return False

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return self.value
