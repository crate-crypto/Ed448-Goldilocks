# Field

This library will make use of two Fields:

- The base prime field `q` which both the Montgomery and Edwards Curve are defined over
- The scalar field which is instantiated from the curve being defined over `q`


- For now we will use 28-bit limbs to define an element in `F_q`. In the future, it would be nice to have an implementation for 32-bit limbs and 64-bit limbs, however this is not a priority.

- For now we have karatsuba under field elements, however in the future we should have a separate module for karatsuba so that we can include other ways to do multiplication.