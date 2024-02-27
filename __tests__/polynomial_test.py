from hyperelliptic import Polynomial


def test_constructors():
    assert Polynomial([1, 0, 0, 4]) == Polynomial([0, 0, 1, 0, 0, 4])
    assert Polynomial([3, 0, 1, 0, 2, 2]).deg == 5


def test_addition():
    p1 = Polynomial([1, 0, 1, 8, 9, 7])
    p2 = Polynomial([1, 0, 3, 4])
    assert p1 + p2 == Polynomial([1, 0, 2, 8, 12, 11])


def test_substraction():
    p1 = Polynomial([1, 0, 1, 8, 9, 7])
    p2 = Polynomial([1, 0, 3, 4])
    assert p2 - p1 == Polynomial([-1, 0, 0, -8, -6, -3])
    assert p1 - p2 == Polynomial([1, 0, 0, 8, 6, 3])


def test_multiplication():
    p1 = Polynomial([1, 0, 1, 8, 9, 7])
    p2 = Polynomial([1, 0, 3, 4])
    p3 = Polynomial([1, 0, 10])
    p4 = Polynomial([7, 9])
    assert p1 * p2 == Polynomial([1, 0, 4, 12, 12, 35, 59, 57, 28])
    assert p2 * p1 == Polynomial([1, 0, 4, 12, 12, 35, 59, 57, 28])
    assert p3 * p4 == Polynomial([7, 9, 70, 90])


def test_division():
    p1 = Polynomial([1, 4, 2, 8, 1, 4])
    p2 = Polynomial([1, 4, 1, 4])
    p3 = Polynomial([1, 0, 1, 8, 9, 7])
    p4 = Polynomial([1, 0, 3, 4])

    assert divmod(p1, p2) == (Polynomial([1, 0, 1]), Polynomial([0]))
    assert divmod(p3, p4) == (Polynomial([1, 0, -2]), Polynomial([4, 15, 15]))

    assert p3 // p4 == Polynomial([1, 0, -2])
    assert p3 % p4 == Polynomial([4, 15, 15])


def test_comparison():
    p1 = Polynomial([1, 0, 10])
    assert p1 != 10
    assert p1 == Polynomial([0, 1, 0, 10])


def test_exponentiation():
    p1 = Polynomial([1, 0, 1])
    assert p1**2 == Polynomial([1, 0, 2, 0, 1])
    assert p1**3 == Polynomial([1, 0, 3, 0, 3, 0, 1])
    assert p1 * p1 * p1 * p1 == p1**4
    assert p1**0 == 1


def test_to_string():
    assert str(Polynomial([1, 0, 3, 1, 1])) == "x^4 + 3x^2 + x + 1"
    assert str(Polynomial([1, 0, 3, 1, 0])) == "x^4 + 3x^2 + x"
    assert (
        str(
            Polynomial(
                [
                    Polynomial([1, 0, 3, 1, 1]),
                    Polynomial([1, 0, 3, 1, 0]),
                    Polynomial([0]),
                ],
                "y",
            )
        )
        == "(x^4 + 3x^2 + x + 1)y^2 + (x^4 + 3x^2 + x)y"
    )
    assert (
        str(
            Polynomial(
                [
                    Polynomial([1]),
                    Polynomial([1, 0, 3, 1, 1]),
                    Polynomial([1, 0, 3, 1, 0]),
                    Polynomial([0]),
                ],
                "y",
            )
        )
        == "y^3 + (x^4 + 3x^2 + x + 1)y^2 + (x^4 + 3x^2 + x)y"
    )


def test_repr():
    assert repr(Polynomial([1, 0, 3, 1, 1])) == "x^4 + 3x^2 + x + 1"
    assert repr(Polynomial([1, 0, 3, 1, 0])) == "x^4 + 3x^2 + x"


def test_gcd():
    p1 = Polynomial([1, 7, 6])
    p2 = Polynomial([1, 6, 5])

    p3 = Polynomial([1, 4, 1, 4])
    p4 = Polynomial([1, 0, 1])

    assert p1.gcd(p2) == Polynomial([1, 1])
    assert p3.gcd(p4) == p4
