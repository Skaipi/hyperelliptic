import pytest
from hyperelliptic import FiniteField


def test_constructor():
    gf = FiniteField(11)
    poly = gf.poly([1, 3, 3])
    gf = gf.extension(poly)

    gf = FiniteField(11)
    poly = gf.poly([1, 0, 10])

    with pytest.raises(ValueError):
        gf.extension(poly)


def test_different_gf():
    g1 = FiniteField(11)
    p1 = g1.poly([1, 3, 3])
    g1 = g1.extension(p1)

    g2 = FiniteField(23)
    p2 = g2.poly([1, 1, 12])
    g2 = g2.extension(p2)

    g3 = FiniteField(11)
    p3 = g3.poly([1, 3, 3])
    g3 = g3.extension(p3)

    g4 = FiniteField(23)
    p4 = g4.poly([1, 1, 12])
    g4 = g4.extension(p4)

    a = g1.element([1, 5])
    b = g2.element([1, 5])
    c = g3.element([1, 7])
    d = g4.element([1, 7])

    with pytest.raises(ValueError):
        a + b

    assert a + c == g3([2, 1])
    assert b + d == g4([2, 12])


def test_to_string():
    gf = FiniteField(11)
    poly = gf.poly([1, 3, 3])
    gf = gf.extension(poly)

    assert str(gf) == "Galois Field mod 11 mod x^2 + 3x + 3"
