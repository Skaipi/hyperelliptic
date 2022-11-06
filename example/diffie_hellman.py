import sys

sys.path.insert(0, "..")

from src.gf import GaloisField

if __name__ == "__main__":
    SMALL_PRIME = 884666614024826252892955729547
    MEDIUM_PRIME = 728332861387732709516448268243094614312200863702341084222463
    BIG_PRIME = 963438379803025380029290083247131353689695039066556841283773521973356839991301911935703108462287808104585159620942384189

    gf = GaloisField(BIG_PRIME)
    f = gf.poly([1, 0, 0, 0, 0, 17])
    h = gf.poly([0])
    c = gf.hyperelliptic(h, f)
    g = c.get_random_divisor()

    alice_number = gf.rand_int()
    alice_element = g * alice_number

    bob_number = gf.rand_int()
    bob_element = g * bob_number

    alice_common_element = bob_element * alice_number
    bob_common_element = alice_element * bob_number

    print(alice_common_element == bob_common_element)
    print(bob_common_element)
