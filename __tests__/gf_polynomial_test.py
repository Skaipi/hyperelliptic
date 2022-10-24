import pytest
from src.gf import GaloisField


def test_constructors():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    assert gf.poly([1, 0, 1, 0, 0, 0]) == gf.poly([1, 1, 0, 1])
    assert gf.poly([1, 0, 0, 1, 0, 0]) == gf.poly_one
    assert gf.poly([1, 0, 0, 1, 0, 1]) == gf.poly_zero


def test_addition():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p1 = gf.poly([1, 0, 1, 0, 1, 1])
    p2 = gf.poly([1, 0, 1, 0])
    assert p1 + p2 == gf.poly([1, 0, 0, 0, 0, 1])
    assert p1 + p2 == gf.poly([1, 0, 0])
    assert p1 + gf.poly_zero == p1
    assert p1 + gf.poly_one == gf.poly([1, 0, 1, 0, 1, 0])


def test_substraction():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p1 = gf.poly([1, 1, 1, 0])
    p2 = gf.poly([1, 0, 1, 0])
    assert p1 + p2 == p1 - p2
    assert p1 - gf.poly_one == gf.poly([1, 1, 1, 1])
    assert p1 - gf.poly_zero == p1


def test_multiplication():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p1 = gf.poly([1, 0, 1, 0, 0])
    p2 = gf.poly([1, 0])
    assert p1 * p2 == gf.poly([1, 1, 0, 1])
    assert p1 * gf.poly_one == p1
    assert p1 * gf.poly_zero == gf.poly_zero


def test_division():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p1 = gf.poly([1, 0, 1, 0, 0])
    p2 = gf.poly([1, 0, 1])
    p3 = gf.poly([1, 0])
    assert divmod(p1, p2) == (gf.poly([1, 0, 0]), gf.poly_zero)
    assert divmod(p1, p3) == (gf.poly([1, 0, 1, 0]), gf.poly_zero)
    assert divmod(p2, p3) == (gf.poly([1, 0]), gf.poly_one)


def test_inversion():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p1 = gf.poly([1, 0, 1, 0, 0])
    assert p1.inverse() == gf.poly([1, 1, 1, 1, 0])
    assert p1.inverse() * p1 == gf.poly_one

    with pytest.raises(ZeroDivisionError):
        gf.poly_zero.inverse()


def test_exponentiation():
    gf = GaloisField(2, 5, [1, 0, 0, 1, 0, 1])

    p1 = gf.poly([1, 0])
    assert p1**0 == gf.poly_one
    assert p1**1 == p1
    assert p1**6 == gf.poly([1, 0, 1, 0])
    assert p1**15 == gf.poly([1, 1, 1, 1, 1])
