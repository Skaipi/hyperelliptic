import pytest
from src.gf import GaloisField


def test_division():
    gf = GaloisField(11)
    p = gf.poly([4, 0, 4, 10, 3, 6])
    x = gf(3)
    assert p * x == gf.poly([1, 0, 1, 8, 9, 7])


def test_multiplication():
    gf = GaloisField(11)
    p = gf.poly([1, 0, 1, 8, 9, 7])
    x = gf(3)

    assert p / x == gf.poly([4, 0, 4, 10, 3, 6])

    with pytest.raises(ZeroDivisionError):
        p / gf.int_zero
