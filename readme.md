# Goldilocks

VERY BIG WIP

## Field Choice

-  The field size is a Solinas trinomial prime 2^448 - 2^224 -1. This prime is called the Goldilocks prime.
- We will denote the Goldilocks prime with the variable `q`.

## Goldilocks Curve

The goldilocks curve is an untwisted edwards curve. x^2 + y^2 = 1 -39081x^2y^2.

## Twisted Curve

This library will also implement the Twisted variation by using a dual isogeny which is a variation of the affine doubling formula. This dual isogeny will effectively clear the cofactor for the computation on the Twisted curve, and it will allow us to use the faster twisted curve arithmetics. However, this strategy does not clear the cofactor entirely. If the point is of low order, then the Twisted curve arithmetic will produce the identity point, while the remaining untwisted curve arithmetic will produce a point in the small order subgroup.  

Computing the dual isogeny to clear the cofactor has the same cost as clearing the cofactor through doubling. [Link paper]

## Decaf/Ristretto

Due to the rationale in the Twisted Curve section, this library will also implement Decaf and Ristretto over the Twisted Edwards Curve.

## Completed Point vs Extensible Point

Deviating from Curve25519-Dalek, this library will implement Extensible points instead of Completed Points. Due to the following observation:

- There is a cost of 3/4 Field multiplications to switch from the CompletedPoint. SO if we were to perform repeated doubling, this would add an extra cost for each doubling.


## Interesting Notes


Goldilocks is defined using an untwisted edwards curve, Curve25519 is defined using a Montgomery Curve. Both use an isogeny to compute the corresponding Twisted Edwards  

Perhaps for interopability in naming:

- CurveXXX = MontgomeryCurve