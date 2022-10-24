import pytest
from src.gf import GaloisField


def test_comparison():
    gf = GaloisField(7)
    x = gf(5)

    assert x == 5
    assert x < 6
    assert x > 4


def test_overflow_contructor_comparison():
    gf = GaloisField(7)
    x = gf(9)

    assert x == 2
    assert x < 3
    assert x > 1


def test_addition():
    gf = GaloisField(11)

    assert gf(4) + gf(5) == 9
    assert gf(7) + gf(7) == 3
    assert gf(2) + gf(-2) == 0
    assert gf(13) + gf(21) == 1
    assert gf(4) + 5 == 9


def test_substraction():
    gf = GaloisField(11)

    assert gf(4) - gf(5) == -1
    assert gf(4) - gf(5) == 10
    assert gf(2) - gf(2) == 0
    assert gf(10) - gf(122) == 9
    assert gf(4) - 5 == -1


def test_multiplication():
    gf = GaloisField(11)

    assert gf(4) * gf(4) == 5
    assert gf(7) * gf(7) == 5
    assert gf(8) * gf(1) == 8
    assert gf(125) * gf(0) == 0
    assert gf(4) * 4 == 5


def test_division():
    gf = GaloisField(11)

    assert gf(2) / gf(4) == 6
    assert gf(7) / gf(6) == 3
    assert gf(2) / 4 == 6

    with pytest.raises(ZeroDivisionError):
        gf(7) / gf(11)


def test_exponentiation():
    gf = GaloisField(11)

    assert gf(2) ** 2 == 4
    assert gf(7) ** 11 == 7


def test_sqrt():
    gf = GaloisField(11)

    assert gf.sqrt(gf(0)) == 0
    assert gf.sqrt(gf(1)) == 1
    assert gf.sqrt(gf(3)) == 5

    with pytest.raises(ValueError):
        gf.sqrt(gf(2))


def test_negation():
    gf = GaloisField(11)

    assert -gf(0) == 0
    assert -gf(1) == 10
    assert -gf(10) == 1


def test_hash():
    gf = GaloisField(11)

    random_integers = [gf.rand_int() for i in range(20)]
    list(set(random_integers))
