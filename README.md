# Marlin-Math

Accurate and fastest Math functions in java, like the Marlin renderer !

Rationale:
Java supports Quadratic & Cubic curves in Java2D & JavaFX Graphics, that require good essential math functions:
- Polynoms of 1 to 3 degrees (line, quad, cubic)
- 2nd & 3rd order root finder
- extrema, inflection points
- intersections

Double-precision does not mean perfect accuracy, but BigDecimal impl is too costly.

This project aims testing accuracy to get the upper bound of numerical errors, possibly use solver refinement to increase accuracy.
However these functions must be the fastest possible ones according to the accuracy requirements.
