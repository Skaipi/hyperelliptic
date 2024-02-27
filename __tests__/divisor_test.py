from hyperelliptic import FiniteField


def test_constructor():
    gf = FiniteField(5)

    f = gf.poly([1, 0, 3, 2, 0, 3])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    p1 = (gf(3), gf(0))
    p2 = (gf(1), gf(2))
    q1 = (gf(4), gf(1))
    q2 = (gf(3), gf(0))

    dz = c.divisor_from_points([("Inf", "Inf"), ("Inf", "Inf")])
    d1 = c.divisor_from_points([p1, p2])
    d2 = c.divisor_from_points([q1, q2])
    d3 = c.divisor_from_points([q1, q1])

    expected_dz = c.zero_divisor()
    expected_d1 = c.divisor(gf.poly([1, 1, 3]), gf.poly([4, 3]))
    expected_d2 = c.divisor(gf.poly([1, 3, 2]), gf.poly([1, 2]))
    expected_d3 = c.divisor(gf.poly([1, 2, 1]), gf.poly([1]))

    assert dz == expected_dz
    assert d1 == expected_d1
    assert d2 == expected_d2
    assert d3 == expected_d3

    assert (d1.v**2 + d1.v * d1.c.h - d1.c.f) % d1.u == 0
    assert (d2.v**2 + d2.v * d2.c.h - d2.c.f) % d2.u == 0
    assert (d3.v**2 + d3.v * d3.c.h - d3.c.f) % d3.u == 0


def test_points_from_divisor():
    gf = FiniteField(5)

    f = gf.poly([1, 0, 3, 2, 0, 3])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    dz = c.zero_divisor()
    d1 = c.divisor(gf.poly([1, 1, 3]), gf.poly([4, 3]))
    d2 = c.divisor(gf.poly([1, 3, 2]), gf.poly([1, 2]))
    d3 = c.divisor(gf.poly([1, 2, 1]), gf.poly([1]))

    assert dz.points == [("Inf", "Inf"), ("Inf", "Inf")]
    assert all(p in [(3, 0), (1, 2)] for p in d1.points) and len(d1.points) == 2
    assert all(p in [(4, 1), (3, 0)] for p in d2.points) and len(d2.points) == 2
    assert all(p in [(4, 1), (4, 1)] for p in d3.points) and len(d3.points) == 2


def test_random_divisors():
    gf = FiniteField(5)

    f = gf.poly([1, 0, 3, 2, 0, 3])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    for _ in range(10):
        c.get_random_divisor()


def test_addition():
    gf = FiniteField(11)

    f = gf.poly([1, 0, 3, 7, 1, 2])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    d1 = c.divisor(gf.poly([1, 7, 10]), gf.poly([1, 9]))
    d2 = c.divisor(gf.poly([1, 0, 10]), gf.poly([7, 9]))
    d3 = c.divisor(gf.poly([1, 10]), gf.poly([6]))
    dz = c.zero_divisor()

    assert d1 + d2 == d3
    assert d1 + dz == d1
    assert d2 + dz == d2
    assert d3 + dz == d3
    assert dz + dz == dz
