import pytest
from src.gf import GaloisField


def test_addition():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p = gf.poly([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p + o == gf.poly([1, 0, 1, 1])
    assert o + p == gf.poly([1, 0, 1, 1])
    assert p + z == p
    assert z + p == p


def test_substraction():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p = gf.poly([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p - o == p + o
    assert o - p == -p + o
    assert p - z == p + z
    assert z - p == -p + z


def test_multiplication():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p = gf.poly([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p * o == p
    assert p * z == gf.poly_zero
    assert o * p == p
    assert z * p == gf.poly_zero

    gf = GaloisField(11)
    p = gf.poly([4, 0, 4, 10, 3, 6])
    x = gf(3)
    assert p * x == gf.poly([1, 0, 1, 8, 9, 7])


def test_division():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p = gf.poly([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p / o == p

    with pytest.raises(ZeroDivisionError):
        p / z

    gf = GaloisField(11)
    p = gf.poly([1, 0, 1, 8, 9, 7])
    x = gf(3)
    assert p / x == gf.poly([4, 0, 4, 10, 3, 6])
