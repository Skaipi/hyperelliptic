# hyperalliptic

<div align=center>
  <a href="https://github.com/Skaipi/hyperalliptic/actions/workflows/test.yml"><img src="https://github.com/Skaipi/hyperalliptic/actions/workflows/test.yml/badge.svg"></a>
</div>

Set of tools for operations on hyperelliptic curves and finite fields.

## Installation
To install as python module:
 - download this repo
 - navigate to repo's directory
 - run the following command to install from source
```
python -m pip install -e ./
```

## Example
```python
# Diffie-Hellman algorithm
gf = GaloisField(PRIME_NUMBER)
h = gf.poly([0])
f = gf.poly(COEFFICIENTS)
c = gf.hyperelliptic(h, f)
g = c.get_random_divisor()

alice_secret = gf.rand_int()
alice_element = g * alice_secret

bob_secret = gf.rand_int()
bob_element = g * bob_secret

common_element = alice_element * bob_sercret
```
