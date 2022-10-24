import math


class ZP:
    def __init__(self, gf, value):
        self.gf = gf
        self.p = gf.p
        self.value = value % self.p

    def _from_value(self, value):
        return ZP(self.gf, value)

    def zero(self):
        return self._from_value(0)

    def one(self):
        return self._from_value(1)

    def gcd(self, other):
        return math.gcd(self.value, other.value)

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

    def __add__(self, other):
        if isinstance(other, int):
            return self + self._from_value(other)
        if not isinstance(other, ZP):  # Avoid circular dependency with Polynomials
            return other + self
        return self._from_value(self.value + other.value)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, int):
            return self - self._from_value(other)
        if not isinstance(other, ZP):  # Avoid circular dependency with Polynomials
            return other - self
        return self._from_value(self.value - other.value)

    def __rsub__(self, other):
        return self.__sub__(other)

    def __mul__(self, other):
        if isinstance(other, int):
            return self * self._from_value(other)
        if not isinstance(other, ZP):  # Avoid circular dependency with Polynomials
            return other * self
        return self._from_value(self.value * other.value)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __div__(self, other):
        if isinstance(other, int):
            return self._from_value(
                self.value * self._from_value(other).inverse().value
            )
        return self._from_value(self.value * other.inverse().value)

    def __truediv__(self, other):
        return self.__div__(other)

    def __divmod__(self, other):
        if isinstance(other, int):
            q, r = divmod(self.value, other)
            return self._from_value(q), self._from_value(r)
        q, r = divmod(self.value, other.value)
        return self._from_value(q), self._from_value(r)

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __pow__(self, other):
        if isinstance(other, int):
            return self._from_value(pow(self.value, other, self.p))
        raise NotImplementedError("")

    def __neg__(self):
        return self._from_value(-self.value)

    def __gt__(self, other):
        if isinstance(other, int):
            return self.value > (other % self.p)
        elif isinstance(other, ZP):
            return self.value > other.value
        raise Exception("")

    def __lt__(self, other):
        if isinstance(other, int):
            return self.value < (other % self.p)
        elif isinstance(other, ZP):
            return self.value < other.value
        raise Exception("")

    def __eq__(self, other):
        if isinstance(other, int):
            return self.value == (other % self.p)
        elif isinstance(other, ZP):
            return self.value == other.value
        return False

    def __str__(self):
        return str(self.value)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return self.value
