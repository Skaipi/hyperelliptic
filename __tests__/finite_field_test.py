import pytest
from hyperelliptic import FiniteField


def test_constructor():
    prime = 2**16 + 1
    trivial_gf = FiniteField(2)
    valid_gf = FiniteField(prime)

    with pytest.raises(ValueError):
        FiniteField(prime - 1)


def test_different_gf():
    g1 = FiniteField(11)
    g2 = FiniteField(23)
    g3 = FiniteField(11)
    g4 = FiniteField(23)

    a = g1(5)
    b = g2(5)
    c = g3(7)
    d = g4(7)

    with pytest.raises(ValueError):
        a + b

    assert a + c == FiniteField(11).element(1)
    assert b + d == FiniteField(23).element(12)


def test_random_irreducible_poly():
    gf = FiniteField(11)

    for i in range(2, 5):
        p = gf.rand_irreducible_poly(i)
        assert p.is_irreducible()


def test_to_string():
    gf = FiniteField(11)

    assert str(gf) == "Finite Field mod 11"
