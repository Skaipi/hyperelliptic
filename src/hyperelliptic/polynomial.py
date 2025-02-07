"""(module) containing general implementation of polynomial with coeffiecients from finite field"""

from copy import copy
from .integer import ZP


def is_int_like(obj):
    """Check if argument is integer like"""
    return isinstance(obj, int) or isinstance(obj, ZP)


def is_number_like(obj):
    """Check if argument is numeric type"""
    return is_int_like(obj) or isinstance(obj, float)


def same_type_coeff(function):
    """Decorator that ensures arguments for polynomial method have coefficients with the same type"""

    def wrapper(poly_a, poly_b):
        # pylint: disable=C0123
        if is_number_like(poly_a) or is_number_like(poly_b):
            return function(poly_a, poly_b)

        ca, cb = poly_a.leading_coeff, poly_b.leading_coeff
        is_int_type = is_int_like(ca) or is_int_like(cb)
        if not is_int_type and type(ca) != type(cb):
            raise ValueError(
                f"Types of polynomial {poly_a} and {poly_b} must not differ"
            )
        return function(poly_a, poly_b)

    return wrapper


class Polynomial:
    """General use class implementing basic operations on arbitrary polynomials"""

    def __init__(self, coeff, symbol="x"):
        self.coeff = self._strip(coeff)
        self.symbol = symbol
        self._has_int_coeff = is_int_like(self.leading_coeff)

    @property
    def deg(self):
        """Polynomial degree"""
        return len(self.coeff) - 1

    @property
    def leading_coeff(self):
        """Most signifficant coefficient of a polynomial"""
        return self.coeff[0]

    def zero(self):
        """Polynomial neutral for addition"""
        return self._from_coeff([self.coeff_zero()])

    def one(self):
        """Polynomial neutral to multiplication"""
        return self._from_coeff([self.coeff_one()])

    def coeff_zero(self):
        """Coefficient neutral to addition"""
        return 0

    def coeff_one(self):
        """Coefficient neutral to multiplication"""
        return 1

    def to_monic(self):
        """Reduce polynomial to form where leading coefficient is equal to one"""
        c = self.leading_coeff
        if c != 1:
            return self / c
        return self

    def is_const(self):
        """Check if polynial is degree 0"""
        return self.deg <= 0

    def to_int(self):
        """Cast degree 0 polynomial to base element"""
        if not self.is_const():
            raise TypeError(f"Can not cast non-0 degree poly {self} to an integer")
        return self.leading_coeff

    def derivative(self):
        """Derivative of polynomial"""
        if self.deg == 0:
            return self.zero()
        coeff = [c * (self.deg - i) for i, c in enumerate(self.coeff[:-1])]
        return self._from_coeff(coeff)

    @same_type_coeff
    def gcd(self, other) -> "Polynomial":
        """Euclidean algorithm for polynomials"""
        r1, r0 = self, other
        while r0 != self.zero():
            r1, r0 = r0, r1 % r0

        return r1.to_monic()

    @same_type_coeff
    def xgcd(self, other: "Polynomial"):
        """Extended Euclidean algorithm for polynomials"""
        # pylint: disable=W0212
        r1, r0 = self._copy(), other._copy()
        s1, s0 = self.one(), self.zero()
        t1, t0 = self.zero(), self.one()

        while r0 != self.zero():
            q = r1 // r0
            r1, r0 = r0, r1 - q * r0
            s1, s0 = s0, s1 - q * s0
            t1, t0 = t0, t1 - q * t0

        leading_coeff = r1.leading_coeff
        if leading_coeff != self.one():
            r1 = r1 / leading_coeff
            s1 = s1 / leading_coeff
            t1 = t1 / leading_coeff

        return r1, s1, t1

    def _copy(self):
        return Polynomial(self.coeff)

    def _from_coeff(self, coeff):
        return Polynomial(coeff, self.symbol)

    def _strip(self, coeff) -> list[int | ZP]:
        element = next(filter(lambda x: x != 0, coeff), None)
        return [self.coeff_zero()] if element is None else coeff[coeff.index(element) :]

    @same_type_coeff
    def __add__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            coeff = copy(self.coeff)  # Preventing mutation on self.coeff
            coeff[-1] = coeff[-1] + other
            return self._from_coeff(coeff)

        zero_coeff = self.coeff_zero()
        size = max(self.deg, other.deg)
        s_coeff = [zero_coeff] * (size - self.deg) + self.coeff
        o_coeff = [zero_coeff] * (size - other.deg) + other.coeff

        result = list(map(lambda s, o: s + o, s_coeff, o_coeff))
        return self._from_coeff(result)

    @same_type_coeff
    def __radd__(self, other):
        return self.__add__(other)

    @same_type_coeff
    def __sub__(self, other):
        return self + (-other)

    @same_type_coeff
    def __rsub__(self, other):
        return -self + other

    @same_type_coeff
    def __mul__(self, other):
        if is_number_like(other):
            return self._from_coeff(list(map(lambda x: x * other, self.coeff)))

        result = [self.coeff_zero()] * (self.deg + other.deg + 1)
        for e1, c1 in enumerate(self.coeff):
            for e2, c2 in enumerate(other.coeff):
                result[e1 + e2] += c1 * c2
        return self._from_coeff(result)

    @same_type_coeff
    def __rmul__(self, other):
        return self.__mul__(other)

    @same_type_coeff
    def __truediv__(self, other):
        if isinstance(other, self.leading_coeff.__class__) or is_int_like(other):
            return self._from_coeff(list(map(lambda x: x / other, self.coeff)))
        return NotImplemented

    @same_type_coeff
    def __divmod__(self, other):
        if self.deg < other.deg:
            return self.zero(), self
        if other.is_const():
            return self / other.to_int(), self.zero()

        zero = self.coeff_zero()
        one = self.coeff_one()
        remainder = self._copy()
        quotient = Polynomial([zero])

        while remainder != self.zero() and other.deg <= remainder.deg:
            t = remainder.leading_coeff / other.leading_coeff
            m = Polynomial([one] + [zero] * (remainder.deg - other.deg))
            quotient = quotient + t * m
            remainder = remainder - t * other * m

        return self._from_coeff(quotient.coeff), self._from_coeff(remainder.coeff)

    @same_type_coeff
    def __mod__(self, other):
        return divmod(self, other)[1]

    @same_type_coeff
    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __pow__(self, other, mod=None):
        if not is_int_like(other):
            raise NotImplementedError(f"Polynomial exp error (exp = {other})")

        sq = self
        exp = other if isinstance(other, int) else other.value
        result = self.one()

        while exp > 0:
            if exp % 2 == 1:
                result = result * sq
                if mod is not None:
                    result = result % mod
            sq *= sq
            if mod is not None:
                sq = sq % mod
            exp = exp // 2
        return result

    def __neg__(self):
        return self._from_coeff(list(map(lambda x: -x, self.coeff)))

    def __eq__(self, other):
        if isinstance(other, int):
            return len(self.coeff) == 1 and self.leading_coeff == other
        if isinstance(other, Polynomial):
            if self.deg != other.deg:
                return False
            return all(a == b for a, b in zip(self.coeff, other.coeff))
        return False

    def __repr__(self):
        return str(self)

    def __call__(self, x):
        return sum(
            [x**i * c for i, c in enumerate(self.coeff[::-1])], self.coeff_zero()
        )

    def __str__(self):
        if all(x == 0 for x in self.coeff):
            return "0"

        def wrap_expr(expr):
            if self._has_int_coeff:
                return expr
            return f"({expr})"

        def get_nth_expr(n, c):
            if c != 1:
                return f"{wrap_expr(c)}{self.symbol}{append_pow(n)}"
            return f"{self.symbol}{append_pow(n)}"

        # pylint: disable=C3001
        result = ""
        first = 0
        last = self.deg
        append_pow = lambda exp: "" if exp == last - 1 else f"^{last-i}"

        for i, c in enumerate(self.coeff):
            if c == 0:
                continue
            elif i == last:
                result += f" + {c}" if first != last else f"{c}"
            elif i == first:
                result += get_nth_expr(i, c)
            else:
                result += f" + {get_nth_expr(i, c)}"
        return result
