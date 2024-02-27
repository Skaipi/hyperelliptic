import pytest
from hyperelliptic import FiniteField


def test_addition():
    gf = FiniteField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p = gf.element([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p + o == gf.element([1, 0, 1, 1])
    assert o + p == gf.element([1, 0, 1, 1])
    assert p + z == p
    assert z + p == p


def test_substraction():
    gf = FiniteField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p = gf.element([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p - o == p + o
    assert o - p == -p + o
    assert p - z == p + z
    assert z - p == -p + z


def test_multiplication():
    gf = FiniteField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p = gf.element([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p * o == p
    assert p * z == gf.zero()
    assert o * p == p
    assert z * p == gf.zero()


def test_division():
    gf = FiniteField(2)
    poly = gf.poly([1, 0, 0, 1, 0, 1])
    gf = gf.extension(poly)

    p = gf.element([1, 0, 1, 0])
    o = gf(1)
    z = gf(0)
    assert p / o == p

    with pytest.raises(ZeroDivisionError):
        p / z
