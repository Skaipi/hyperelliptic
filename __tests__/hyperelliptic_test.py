import pytest
from src.finite_field import FiniteField


def test_constructor():
    gf = FiniteField(11)

    f = gf.poly([2, 0, 3, 7, 1, 2])
    h = gf.poly([0])

    with pytest.raises(ValueError) as e:
        gf.hyperelliptic(h, f)
    assert str(e.value) == "Function f must be monic"

    f = gf.poly([1, 3, 7, 1, 2])

    with pytest.raises(ValueError) as e:
        gf.hyperelliptic(h, f)
    assert str(e.value) == "Degree of function f must be an odd number"

    f = gf.poly([1, 0, 3, 7, 1, 2])
    h = gf.poly([1, 0, 0])

    with pytest.raises(ValueError) as e:
        gf.hyperelliptic(h, f)
    assert str(e.value) == "h(x) must be 0 for char(F) != 2"


def test_points_on_curve():
    gf = FiniteField(11)

    f = gf.poly([1, 0, 3, 7, 1, 2])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    expected_points = [
        ("Inf", "Inf"),
        (1, 6),
        (1, 5),
        (2, 0),
        (4, 6),
        (4, 5),
        (6, 7),
        (6, 4),
        (7, 7),
        (7, 4),
        (9, 7),
        (9, 4),
        (10, 2),
        (10, 9),
    ]
    points = c.get_all_points()

    assert all([p in expected_points for p in points]) == True
    assert len(points) == len(expected_points)

    gf = FiniteField(5)

    f = gf.poly([1, 0, 3, 7, 1, 2])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    expected_points = [("Inf", "Inf"), (1, 3), (1, 2), (4, 2), (4, 3)]
    points = c.get_all_points()

    assert all([p in expected_points for p in points]) == True
    assert len(points) == len(expected_points)


def test_points_on_curve_over_gf():
    gf = FiniteField(11)
    poly = gf.poly([1, 3, 3])
    gf = gf.extension(poly)

    f = gf.poly([1, 0, 3, 7, 1, 2])
    h = gf.poly([0])

    c = gf.hyperelliptic(h, f)
    points = c.get_all_points()

    assert len(points) == 147


def test_non_zero_divisors_on_curve():
    gf = FiniteField(5)

    f = gf.poly([1, 0, 3, 7, 1, 2])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    expected_divisors = [  # Note those are not unique divisors
        c.divisor(gf.poly([1, 3, 1]), gf.poly([2])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([2, 0])),
        c.divisor(gf.poly([1, 4]), gf.poly([2])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([2])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([2, 0])),
        c.divisor(gf.poly([1, 2, 1]), gf.poly([3])),
        c.divisor(gf.poly([1, 1]), gf.poly([3])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([3])),
        c.divisor(gf.poly([1, 4]), gf.poly([2])),
        c.divisor(gf.poly([1, 1]), gf.poly([3])),
        c.divisor(gf.poly([1, 1]), gf.poly([2])),
        c.divisor(gf.poly([1, 4]), gf.poly([3])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([2])),
        c.divisor(gf.poly([1, 1]), gf.poly([2])),
        c.divisor(gf.poly([1, 2, 1]), gf.poly([2])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([3])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([3, 0])),
        c.divisor(gf.poly([1, 4]), gf.poly([3])),
        c.divisor(gf.poly([1, 0, 4]), gf.poly([3, 0])),
        c.divisor(gf.poly([1, 3, 1]), gf.poly([3])),
    ]
    divisors = []
    points = c.get_all_points()
    unique_points = list(set(points))  # remove special points (p == -p)

    for p1 in unique_points:
        for p2 in unique_points:
            if p1 == c.point_inverse(p2):
                continue
            divisors.append(c.divisor_from_points([p1, p2]))

    assert all([d in expected_divisors for d in divisors])
    assert len(divisors) == len(expected_divisors)
