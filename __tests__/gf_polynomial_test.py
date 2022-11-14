import pytest
from src.gf import GaloisField


def test_constructors():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    assert gf.element([1, 0, 1, 0, 0, 0]) == gf.element([1, 1, 0, 1])
    assert gf.element([1, 0, 0, 1, 0, 0]) == gf.one
    assert gf.element([1, 0, 0, 1, 0, 1]) == gf.zero


def test_addition():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 0, 1, 0, 1, 1])
    p2 = gf.element([1, 0, 1, 0])
    assert p1 + p2 == gf.element([1, 0, 0, 0, 0, 1])
    assert p1 + p2 == gf.element([1, 0, 0])
    assert p1 + gf.zero == p1
    assert p1 + gf.one == gf.element([1, 0, 1, 0, 1, 0])


def test_substraction():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 1, 1, 0])
    p2 = gf.element([1, 0, 1, 0])
    assert p1 + p2 == p1 - p2
    assert p1 - gf.one == gf.element([1, 1, 1, 1])
    assert p1 - gf.zero == p1


def test_multiplication():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 0, 1, 0, 0])
    p2 = gf.element([1, 0])
    assert p1 * p2 == gf.element([1, 1, 0, 1])
    assert p1 * gf.one == p1
    assert p1 * gf.zero == gf.zero


def test_division():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 0, 1, 0, 0])
    p2 = gf.element([1, 0, 1])
    p3 = gf.element([1, 0])

    assert divmod(p1, p2) == (gf.element([1, 0, 0]), gf.zero)
    assert divmod(p1, p3) == (gf.element([1, 0, 1, 0]), gf.zero)
    assert divmod(p2, p3) == (gf.element([1, 0]), gf.one)

    assert p1 % p2 == gf.zero
    assert p1 % p3 == gf.zero
    assert p2 % p3 == gf.one

    assert p1 // p2 == gf.element([1, 0, 0])
    assert p1 // p3 == gf.element([1, 0, 1, 0])
    assert p2 // p3 == gf.element([1, 0])


def test_inversion():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 0, 1, 0, 0])
    assert p1.inverse() == gf.element([1, 1, 1, 1, 0])
    assert p1.inverse() * p1 == gf.one

    with pytest.raises(ZeroDivisionError):
        gf.zero.inverse()


def test_exponentiation():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 0])
    assert p1**0 == gf.one
    assert p1**1 == p1
    assert p1**6 == gf.element([1, 0, 1, 0])
    assert p1**15 == gf.element([1, 1, 1, 1, 1])


def test_to_monic():
    gf = GaloisField(7)
    poly = gf.poly([1, 6, 0, 4])
    gf = gf.extension(poly)

    p1 = gf.element([3, 4, 2])
    assert p1.to_monic() == gf.element([1, 6, 3])


def test_to_string():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 0, 1, 0, 0])
    assert str(p1) == "a^4 + a^2"
    assert str(gf.element([1])) == "1"
    assert str(gf.element([0])) == "0"


def test_repr():
    gf = GaloisField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p1 = gf.element([1, 0, 1, 0, 0])
    assert repr(p1) == "a^4 + a^2"
    assert repr(gf.element([1])) == "1"
    assert repr(gf.element([0])) == "0"
