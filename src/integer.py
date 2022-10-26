import math

from src.utils import gf_operation


class ZP:
    def __init__(self, gf, value):
        self.gf = gf
        self.p = gf.p
        self.value = value % self.p

    def zero(self):
        return self._from_value(0)

    def one(self):
        return self._from_value(1)

    @gf_operation
    def gcd(self, other):
        return math.gcd(self.value, other.value)

    @gf_operation
    def xgcd(self, other):
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

    def inverse(self):
        if self == 0:
            raise ZeroDivisionError(f"Element 0 is not inversable")

        r, s, t = self.xgcd(self.p)
        if r > 1:
            raise Exception(f"Element {self.value} is not inversable mod {self.p}")
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
    def __div__(self, other):
        if isinstance(other, int):
            return self._from_value(
                self.value * self._from_value(other).inverse().value
            )
        return self._from_value(self.value * other.inverse().value)

    @gf_operation
    def __truediv__(self, other):
        return self.__div__(other)

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

    def __pow__(self, other):
        if isinstance(other, int):
            return self._from_value(pow(self.value, other, self.p))
        raise NotImplementedError("")

    def __neg__(self):
        return self._from_value(-self.value)

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
