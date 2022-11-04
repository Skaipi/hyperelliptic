import pytest
from src.gf import GaloisField


def test_constructor():
    prime = 2**16 + 1
    trivial_gf = GaloisField(2)
    valid_gf = GaloisField(prime)

    with pytest.raises(ValueError):
        GaloisField(prime - 1)


def test_different_gf():
    g1 = GaloisField(11)
    g2 = GaloisField(23)
    g3 = GaloisField(11)
    g4 = GaloisField(23)

    a = g1(5)
    b = g2(5)
    c = g3(7)
    d = g4(7)

    with pytest.raises(ValueError):
        a + b

    assert a + c == GaloisField(11)(1)
    assert b + d == GaloisField(23)(12)


def test_get_elements():
    gf = GaloisField(7, 4, [1, 0, 5, 4, 3])

    assert len([_ for _ in gf.get_elements()]) == gf.q

    gf = GaloisField(11)
    assert len([_ for _ in gf.get_elements()]) == gf.q


def test_factorization():
    gf = GaloisField(2)

    expected_209_factors = [11, 19]
    factors = gf.factors(209)
    assert all(f in expected_209_factors for f in factors)
    assert len(factors) == len(expected_209_factors)

    expected_299_factors = [13, 23]
    factors = gf.factors(299)
    assert all(f in expected_299_factors for f in factors)
    assert len(factors) == len(expected_299_factors)

    expected_439_factors = [439]
    factors = gf.factors(439)
    assert all(f in expected_439_factors for f in factors)
    assert len(factors) == len(expected_439_factors)

    expected_64_factors = [2, 2, 2, 2, 2, 2]
    factors = gf.factors(64)
    assert all(f in expected_64_factors for f in factors)
    assert len(factors) == len(expected_64_factors)

    expected_2401_factors = [7, 7, 7, 7]
    factors = gf.factors(2401)
    assert all(f in expected_2401_factors for f in factors)
    assert len(factors) == len(expected_2401_factors)

    too_big_to_factor = 728332861387732709516448268243094614312200863702341084222464
    with pytest.raises(ValueError):
        gf.factors(too_big_to_factor)
