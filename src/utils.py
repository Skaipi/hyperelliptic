from random import randrange


def gf_operation(function):
    def function_wrapper(a, b):
        is_int_operation = isinstance(a, int) or isinstance(b, int)
        if not is_int_operation and a.gf != b.gf:
            raise ValueError(f"{a} field does not match {b} field")
        return function(a, b)

    return function_wrapper


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
