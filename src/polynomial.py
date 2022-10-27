from src.integer import ZP
from copy import copy


def is_int_like(obj):
    return isinstance(obj, int) or isinstance(obj, ZP)


def same_type_coeff(function):
    def wrapper(poly_a, poly_b):
        if is_int_like(poly_a) or is_int_like(poly_b):
            return function(poly_a, poly_b)

        ca, cb = poly_a.coeff[0], poly_b.coeff[0]
        is_int_type = is_int_like(ca) or is_int_like(cb)
        if not is_int_type and type(ca) != type(cb):
            raise ValueError(
                f"Types of polynomial {poly_a} and {poly_b} must not differ"
            )
        return function(poly_a, poly_b)

    return wrapper


class Polynomial:
    def __init__(self, coeff, field=None, symbol="x") -> None:
        self.coeff = self._strip(coeff)
        self.symbol = symbol
        self.gf = field

    @property
    def deg(self) -> int:
        return len(self.coeff) - 1

    def zero(self):
        return self._from_coeff([self.coeff_zero()])

    def one(self):
        return self._from_coeff([self.coeff_one()])

    def to_monic(self):
        c = self.coeff[0]
        if c > 1:
            return self / c
        return self

    def coeff_zero(self):
        if self.gf != None:
            return self.gf.int_zero
        return 0

    def coeff_one(self):
        if self.gf != None:
            return self.gf.int_one
        return 1

    def isConst(self):
        return self.deg <= 0

    def toInt(self):
        if not self.isConst():
            raise Exception(f"Can not cast {self} to integer")
        return self.coeff[0]

    def derivative(self):
        l = len(self.coeff[:-1])
        if l == 0:
            return self.zero()
        coeff = [c * (l - i) for i, c in enumerate(self.coeff[:-1])]
        return self._from_coeff(coeff)

    @same_type_coeff
    def gcd(self, other):
        r1, r0 = self, other
        while r0 != self.zero():
            r1, r0 = r0, r1 % r0

        return r1.to_monic()

    @same_type_coeff
    def xgcd(self, other):
        r1, r0 = self, other
        s1, s0 = self.one(), self.zero()
        t1, t0 = self.zero(), self.one()

        while r0 != self.zero():
            q = r1 // r0
            r1, r0 = r0, r1 - q * r0
            s1, s0 = s0, s1 - q * s0
            t1, t0 = t0, t1 - q * t0

        return r1, s1, t1

    # https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
    def factor(self):
        if self.gf == None:
            raise ValueError(
                "Polynomial must be defined over finite field to be factored"
            )

        factors = []
        sf_factors = self.square_free_factors()

        for sf_factor in sf_factors:
            df_factors = sf_factor.distinct_degree_factors()

            for df_factor, deg in df_factors:
                ed_factors = df_factor.equal_degree_factors(deg)
                factors.extend(ed_factors)

        return factors

    def square_free_factors(self):
        gf = self.gf
        factors = []

        c = self.gcd(self.derivative())
        w = self // c

        # Step one
        i = 1
        while w != self.one():
            y = w.gcd(c)
            fac = w // y
            if fac != self.one() and i % gf.p != 0:
                [factors.append(fac) for _ in range(i)]
            w = y
            c = c // y
            i = i + 1

        # Step two
        if c != self.one():
            coeff = [x ** (gf.p ** (gf.m - 1)) for x in c.coeff]
            factors.extend(gf.poly(coeff).factor())

        return factors

    def distinct_degree_factors(self):
        gf = self.gf
        factors = []
        poly = self._from_coeff(self.coeff)
        x = self._from_coeff([gf(1), gf(0)])

        i = 1
        while i <= self.deg // 2:
            h = pow(x, gf.q, poly)
            g = poly.gcd(h - x)
            if g != self.one():
                factors.append((g, i))
                poly = poly // g
                h = h % poly
            i = i + 1

        if poly != self.one():
            factors.append((poly, poly.deg))
        return factors

    def equal_degree_factors(self, deg):
        gf = self.gf
        factors = [self._from_coeff(self.coeff)]

        while len(factors) < self.deg // deg:
            rand = gf.rand_poly(deg)
            g = self.gcd(rand)
            if g == self.one():
                g = pow(rand, (gf.q**deg - 1) // 2, self) - self.one()
            i = 0
            for fac in factors:
                if fac.deg <= deg:
                    continue
                d = g.gcd(fac)
                if d not in [self.one(), fac]:
                    factors.remove(fac)
                    factors.append(d)
                    factors.append(fac // d)
                i = i + 1

        return factors

    def _from_coeff(self, coeff):
        return Polynomial(coeff, self.gf, self.symbol)

    def _strip(self, arr):
        index = -1
        for i in range(len(arr)):
            if arr[i] != 0:
                index = i
                break
        return arr[index:]

    @same_type_coeff
    def __add__(self, other):
        if isinstance(other, ZP) or isinstance(other, int):
            coeff = copy(self.coeff)  # Preventing mutation on self.coeff
            coeff[-1] = coeff[-1] + other
            return self._from_coeff(coeff)

        result = []
        zero_coeff = self.coeff_zero()
        size = max(self.deg, other.deg)
        s_coeff = [zero_coeff] * (size - self.deg) + self.coeff
        o_coeff = [zero_coeff] * (size - other.deg) + other.coeff

        for i in range(size + 1):
            coeff = s_coeff[i] + o_coeff[i]
            result.append(coeff)
        return self._from_coeff(result)

    @same_type_coeff
    def __sub__(self, other):
        return self + (-other)

    @same_type_coeff
    def __mul__(self, other):
        if is_int_like(other):
            return self._from_coeff(list(map(lambda x: x * other, self.coeff)))

        result = [self.coeff_zero()] * (self.deg + other.deg + 1)
        for e1, c1 in enumerate(self.coeff):
            for e2, c2 in enumerate(other.coeff):
                result[e1 + e2] += c1 * c2
        return self._from_coeff(result)

    @same_type_coeff
    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, exp, mod=None):
        if not isinstance(exp, int):
            raise NotImplementedError(f"Polynomial exp error (exp = {exp})")
        sq = self
        result = self.one()
        while exp > 0:
            if exp & 1:
                result *= sq
            sq *= sq
            exp >>= 1
        return result if mod is None else result % mod

    def __div__(self, other):
        if isinstance(other, self.coeff[0].__class__):
            return self._from_coeff(list(map(lambda x: x / other, self.coeff)))
        raise ValueError("Cannot divide polynomial by other type than int")

    def __truediv__(self, other):
        return self.__div__(other)

    @same_type_coeff
    def __divmod__(self, other):
        if self.deg < other.deg:
            return self.zero(), self
        if other.isConst():
            return self / other.toInt(), self.zero()

        coeff_zero = self.coeff_zero()
        remainder = self._from_coeff(self.coeff)
        quotient = self.zero()

        while remainder != self.zero() and other.deg <= remainder.deg:
            t = remainder.coeff[0] / other.coeff[0]
            m = self._from_coeff([1] + [coeff_zero] * (remainder.deg - other.deg))
            quotient = quotient + t * m
            remainder = remainder - t * other * m

        return quotient, remainder

    @same_type_coeff
    def __mod__(self, other):
        return divmod(self, other)[1]

    @same_type_coeff
    def __floordiv__(self, other):
        return divmod(self, other)[0]

    def __neg__(self):
        return self._from_coeff(list(map(lambda x: -x, self.coeff)))

    def __eq__(self, other):
        if isinstance(other, int):
            return len(self.coeff) == 1 and self.coeff[0] == other
        if isinstance(other, Polynomial):
            if self.deg != other.deg:
                return False
            return all(a == b for a, b in zip(self.coeff, other.coeff))
        return False

    def __repr__(self) -> str:
        return str(self)

    def __call__(self, x):
        return sum([x**i * c for i, c in enumerate(self.coeff[::-1])])

    def __str__(self) -> str:
        if all(x == 0 for x in self.coeff):
            return "0"

        result = ""
        first = 0
        last = self.deg
        append_pow = lambda exp: "" if exp == last - 1 else f"^{last-i}"

        for i, c in enumerate(self.coeff):
            if c == 0:
                continue
            elif i == first and first != last:
                result += f"{c if c != 1 else str()}{self.symbol}{append_pow(i)}"
            elif i == last:
                result += f" + {c}" if first != last else f"{c}"
            else:
                result += f" + {c if c != 1 else str()}{self.symbol}{append_pow(i)}"
        return result
