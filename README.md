Marlin-Math
===========

Accurate and fastest Math functions in java, like the Marlin renderer !


Rationale
=========
Java supports Quadratic & Cubic curves in Java2D & JavaFX Graphics, that require good essential math functions:
- Polynoms of 1 to 3 degrees (line, quad, cubic)
- 2nd & 3rd order root finder
- extrema, inflection points
- intersections

Double-precision does not mean perfect accuracy, but BigDecimal impl is too costly.
See [wikipedia 'Loss of significance'](https://en.wikipedia.org/wiki/Loss_of_significance) that perfectly describe the problem:
- floating-point numbers are inaccurate: ulp(x) is the number representing the smallest difference ~ 1E-16 * x
- math operations (+ - * /) can produce high losses, ie the number is really far from the exact solution >> ulp

This project aims testing accuracy to get the upper bound of numerical errors, possibly use solver refinement to increase accuracy.
However these functions must be the fastest possible ones according to the accuracy requirements.


License
=======

As some code comes from OpenJDK, its license is the OpenJDK's license = GPL2 + ClassPath exception:

GNU General Public License, version 2,
with the Classpath Exception

The GNU General Public License (GPL)

Version 2, June 1991

See License.md


Contributing
============

Contributions are welcomed, please report bugs, fork and share your improvements via pull requests.

Since we contribute parts of this library into OpenJDK, we accept contributions from people who have signed the [Oracle Contribution Agreeement](http://www.oracle.com/technetwork/community/oca-486395.html).


Related projects
================

- [Marlin-Renderer](https://github.com/bourgesl/marlin-renderer) provides the Marlin renderer for Java2D (OpenJDK)
- [Marlin-FX](https://github.com/bourgesl/marlin-fx) provides the Marlin renderer for JavaFX (OpenJFX)
