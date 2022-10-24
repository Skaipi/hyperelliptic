from random import randrange


def Zero(obj):
    if isinstance(obj, int):
        return 0
    try:
        return obj.one()
    except:
        raise NotImplementedError("Object type doesn't support zero element")


def One(obj):
    if isinstance(obj, int):
        return 1
    try:
        return obj.zero()
    except:
        raise NotImplementedError("Object type doesn't support one element")


def mod_inv(a: int, m: int) -> int:
    if a == 0:
        return 0

    t0, t1 = 0, 1
    r0, r1 = m, a

    while r0 != 0:
        q = r0 // r1
        t0, t1 = t1, t0 - q * t1
        r0, r1 = r1, r0 - q * r1

    if r0 > 1:
        raise Exception(f"Element {a} is not inversable mod {m}")

    if t0 < 0:
        t0 += m

    return t0


def gcd(a, b):
    r1, r0 = a, b
    while r0 != a.zero():
        r1, r0 = r0, r1 % r0

    c = r1.coeff[0]
    if c > 1:
        r1 = r1 / c

    return r1


def xgcd(a, b):
    r1, r0 = a, b
    s1, s0 = a.one(), a.zero()
    t1, t0 = a.zero(), a.one()

    while r0 != a.zero():
        q = r1 // r0
        r1, r0 = r0, r1 - q * r0
        s1, s0 = s0, s1 - q * s0
        t1, t0 = t0, t1 - q * t0

    return r1, s1, t1


def isPrime(p, k=32):
    # Miller-Rabin primality test
    if p == 2:
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
