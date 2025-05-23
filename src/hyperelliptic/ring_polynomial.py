"""(module) containing implementation of polynomial ring over arbitrary finite field"""

from .utils import factors as int_factors
from .polynomial import Polynomial


class RingPolynomial(Polynomial):
    """RingPolynomial implements polynomials with coefficients from finite field"""

    def __init__(self, field, coeff, symbol="x"):
        self.gf = field
        super().__init__(coeff, symbol)

    def coeff_zero(self):
        return self.gf.zero()

    def coeff_one(self):
        return self.gf.one()

    def factors(self):
        """
        Returns factors of polynomial. Returns list of poolynomials that divide this element.
        https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
        """
        factors = []
        sf_factors = self.square_free_factors()

        for sf_factor in sf_factors:
            df_factors = sf_factor.distinct_degree_factors()

            for df_factor, deg in df_factors:
                ed_factors = df_factor.equal_degree_factors(deg)
                factors.extend(ed_factors)

        return factors

    def gcd(self, other: "RingPolynomial") -> "RingPolynomial":
        return self._from_coeff(super().gcd(other).coeff)

    def square_free_factors(self) -> list["RingPolynomial"]:
        """Yun's algorithm"""
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
                for _ in range(i):
                    factors.append(fac)
            w = y
            c = c // y
            i = i + 1

        # Step two
        if c != self.one():
            coeff = [x ** (gf.p ** (gf.m - 1)) for x in c.coeff]
            factors.extend(gf.poly(coeff).square_free_factors())

        return factors

    def distinct_degree_factors(self) -> list[tuple["RingPolynomial", int]]:
        """Split a square-free polynomial into a product of polynomials
        whose irreducible factors all have the same degree."""
        gf = self.gf
        factors = []
        poly = self._from_coeff(self.coeff)
        x = self._from_coeff([self.coeff_one(), self.coeff_zero()])

        i = 1
        while i <= self.deg // 2:
            h = pow(x, gf.p, poly)
            g = poly.gcd(h - x)
            if g != self.one():
                factors.append((g, i))
                poly = poly // g
                h = h % poly
            i = i + 1

        if poly != self.one():
            factors.append((poly, poly.deg))
        return factors

    def equal_degree_factors(self, deg) -> list["RingPolynomial"]:
        """Cantor Zassenhaus algorithm"""
        gf = self.gf
        factors = [self._from_coeff(self.coeff)]

        while len(factors) < self.deg // deg:
            rand = gf.rand_poly(deg)
            g = self.gcd(rand)
            if g == self.one():
                g = (rand ** ((gf.p**deg - 1) // 2) - self.one()) % self

            for fac in factors:
                if fac.deg <= deg:
                    continue
                d = g.gcd(fac)
                if d not in [self.one(), fac]:
                    # pylint: disable=W4701
                    factors.remove(fac)
                    factors.append(d)
                    factors.append(fac // d)

        return factors

    def is_irreducible(self):
        """Rabin test of irreducibility"""
        gf = self.gf
        deg_factors = int_factors(self.deg)
        quotients = set([self.deg // factor for factor in deg_factors])
        x = self._from_coeff([self.coeff_one(), self.coeff_zero()])

        prev_h = x
        prev_q = 0
        for quotient in quotients:
            h = pow(prev_h, pow(gf.p, quotient - prev_q), self)
            g = self.gcd((h - x) % self)
            if g != self.one():
                return False
            prev_h, prev_q = h, quotient
        g = pow(prev_h, pow(gf.p, self.deg - prev_q), self)
        g = (g - x) % self

        if g == self.zero():
            return True
        return False

    def _from_coeff(self, coeff):
        return RingPolynomial(self.gf, coeff, self.symbol)

    def __floordiv__(self, other) -> "RingPolynomial":
        return self._from_coeff(super().__floordiv__(other).coeff)
