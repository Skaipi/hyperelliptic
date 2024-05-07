import pytest
from hyperelliptic import factors, all_factors


def test_factorization():
    expected_209_factors = [11, 19]
    factors_209 = factors(209)
    assert all(f in expected_209_factors for f in factors_209)
    assert len(factors_209) == len(expected_209_factors)

    expected_299_factors = [13, 23]
    factors_299 = factors(299)
    assert all(f in expected_299_factors for f in factors_299)
    assert len(factors_299) == len(expected_299_factors)

    expected_439_factors = [439]
    factors_439 = factors(439)
    assert all(f in expected_439_factors for f in factors_439)
    assert len(factors_439) == len(expected_439_factors)

    expected_64_factors = [2, 2, 2, 2, 2, 2]
    factors_64 = factors(64)
    assert all(f in expected_64_factors for f in factors_64)
    assert len(factors_64) == len(expected_64_factors)

    expected_2401_factors = [7, 7, 7, 7]
    factors_2401 = factors(2401)
    assert all(f in expected_2401_factors for f in factors_2401)
    assert len(factors_2401) == len(expected_2401_factors)

    too_big_to_factor = 728332861387732709516448268243094614312200863702341084222464
    with pytest.raises(ValueError):
        factors(too_big_to_factor)


def test_all_factors():
    expected_209_factors = [1, 11, 19, 209]
    factors_209 = all_factors(209)
    assert all(f in expected_209_factors for f in factors_209)
    assert len(factors_209) == len(expected_209_factors)

    expected_439_factors = [1, 439]
    factors_439 = all_factors(439)
    assert all(f in expected_439_factors for f in factors_439)
    assert len(factors_439) == len(expected_439_factors)

    expected_2401_factors = [1, 7, 49, 343, 2401]
    factors_2401 = all_factors(2401)
    assert all(f in expected_2401_factors for f in factors_2401)
    assert len(factors_2401) == len(expected_2401_factors)
