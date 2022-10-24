import sys

sys.path.insert(0, "..")

from src.gf import GaloisField
from ElGamal import encode_divisor


if __name__ == "__main__":
    SMALL_PRIME = 884666614024826252892955729547
    MEDIUM_PRIME = 728332861387732709516448268243094614312200863702341084222463
    BIG_PRIME = 963438379803025380029290083247131353689695039066556841283773521973356839991301911935703108462287808104585159620942384189

    gf = GaloisField(SMALL_PRIME)
    f = gf.poly([1, 0, 3, 2, 0, 3])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)

    d = c.get_random_divisor()
    e = gf.rand_int()
    print(e)
    print(d)
    print(d * e)
