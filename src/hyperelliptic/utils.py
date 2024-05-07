from .primes import PRIMES
from random import randrange
from math import gcd, log, floor
from collections import Counter
from itertools import combinations
from functools import reduce


def gf_operation(function):
    def function_wrapper(a, b):
        is_field_operation = type(a) == type(b)
        if is_field_operation and a.gf != b.gf:
            raise ValueError(f"{a} field does not match {b} field")
        return function(a, b)

    return function_wrapper


def is_prime(p, k=32):
    # Miller-Rabin primality test
    if p < 2:
        return False
    if p == 2:
        return True
    if p == 3:
        return True
    if p % 2 == 0:
        return False

    r, s = 0, p - 1
    while s % 2 == 0:
        r += 1
        s //= 2
    for _ in range(k):
        a = randrange(2, p - 1)
        x = pow(a, s, p)
        if x == 1 or x == p - 1:
            continue
        for _ in range(r - 1):
            x = pow(x, 2, p)
            if x == p - 1:
                break
        else:
            return False
    return True


def all_factors(N):
    result = set([1, N])
    divisors = factors(N)

    for i in range(1, len(divisors) + 1):
        for c in combinations(divisors, i):
            result.add(reduce((lambda x, y: x * y), c))

    return result


def factors(n):
    factors = []

    while n > 1:
        x = factor(n)
        factors.append(x)
        n = n // x

    return factors


def factor(n, smoothness_bound=1):
    # Factor small numbers with pollard's p - 1
    if is_prime(n):
        return n

    max_bound = 1000
    if smoothness_bound > max_bound:
        raise ValueError(f"Can not factor {n}")

    primes = list(filter(lambda x: x <= smoothness_bound, PRIMES))

    m = 1
    for p in primes:
        m *= p ** floor(log(smoothness_bound, p))

    a = 3 if n % 2 == 0 else 2
    g = gcd(pow(a, m, n) - 1, n)

    if 1 < g < n:
        return g

    if g == 1:
        return factor(n, smoothness_bound * 2)

    if g == n:
        raise ValueError(f"Can not factor {n}")
