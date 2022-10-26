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
