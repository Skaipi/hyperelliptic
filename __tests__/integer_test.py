import pytest
from hyperelliptic import FiniteField


def test_comparison():
    gf = FiniteField(7)
    x = gf(5)

    assert x == 5 == gf(5)
    assert x < 6
    assert x < gf(6)
    assert x > 4
    assert x > gf(4)


def test_overflow_contructor_comparison():
    gf = FiniteField(7)
    x = gf(9)

    assert x == 2
    assert x < 3
    assert x > 1


def test_addition():
    gf = FiniteField(11)

    assert gf(4) + gf(5) == 9
    assert gf(7) + gf(7) == 3
    assert gf(2) + gf(-2) == 0
    assert gf(13) + gf(21) == 1
    assert gf(4) + 5 == 9

    gf_2 = FiniteField(7)
    with pytest.raises(ValueError):
        gf(10) + gf_2(7)


def test_substraction():
    gf = FiniteField(11)

    assert gf(4) - gf(5) == -1
    assert gf(4) - gf(5) == 10
    assert gf(2) - gf(2) == 0
    assert gf(10) - gf(122) == 9
    assert gf(4) - 5 == -1

    gf_2 = FiniteField(7)
    with pytest.raises(ValueError):
        gf(10) - gf_2(7)


def test_multiplication():
    gf = FiniteField(11)

    assert gf(4) * gf(4) == 5
    assert gf(7) * gf(7) == 5
    assert gf(8) * gf(1) == 8
    assert gf(125) * gf(0) == 0
    assert gf(4) * 4 == 5

    gf_2 = FiniteField(7)
    with pytest.raises(ValueError):
        gf(10) * gf_2(7)


def test_division():
    gf = FiniteField(11)

    assert gf(2) / gf(4) == 6
    assert gf(7) / gf(6) == 3
    assert gf(2) / 4 == 6

    with pytest.raises(ZeroDivisionError):
        gf(7) / gf(11)

    gf_2 = FiniteField(7)
    with pytest.raises(ValueError):
        gf(10) / gf_2(7)


def test_exponentiation():
    gf = FiniteField(11)

    assert gf(2) ** 2 == 4
    assert gf(7) ** 11 == 7


def test_exponentiation_with_other_element():
    gf = FiniteField(11)

    assert gf(2) ** gf(2) == 4
    assert gf(7) ** gf(11) == 1


def test_sqrt():
    gf = FiniteField(11)

    assert gf(0).sqrt() == 0
    assert gf(1).sqrt() == 1
    assert gf(3).sqrt() == 5

    with pytest.raises(ValueError):
        gf(2).sqrt()


def test_negation():
    gf = FiniteField(11)

    assert -gf(0) == 0
    assert -gf(1) == 10
    assert -gf(10) == 1


def test_hash():
    gf = FiniteField(11)

    random_integers = [gf.rand_element() for _ in range(20)]
    list(set(random_integers))


def test_to_string():
    gf = FiniteField(11)

    assert str(gf(10)) == "10"
    assert str(gf(7)) == "7"
    assert str(gf(5)) == "5"


def test_repr():
    gf = FiniteField(11)

    assert repr(gf(10)) == "10"
    assert repr(gf(7)) == "7"
    assert repr(gf(5)) == "5"
