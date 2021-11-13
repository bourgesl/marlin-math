/*
 * Copyright (c) 2021, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.  Oracle designates this
 * particular file as subject to the "Classpath" exception as provided
 * by Oracle in the LICENSE file that accompanied this code.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 */
package test;

import java.awt.Shape;
import java.awt.geom.CubicCurve2D;
import java.awt.geom.PathIterator;
import java.awt.geom.QuadCurve2D;
import java.awt.geom.Rectangle2D;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import org.marlin.math.util.WelfordVariance;

/**
 * @summary This is not a test. This is an exploratory task to empirically
 *          identify what is the accuracy of the findExtrema() function.
 *          (cubic and quad curves)
 *
 * @author @mickleness @bourgesl
 */
@SuppressWarnings("UseOfSystemOutOrSystemErr")
public class FindExtremaAccuracyTest {

    private final static boolean DEBUG = false;

    private final static int N = 1000000;
    private final static int M = N / 10;

    /** statistics on condition-number */
    private static final WelfordVariance condStats = new WelfordVariance();
    /** statistics on ratio(delta / ucond) */
    private static final WelfordVariance ratioCondStats = new WelfordVariance();

    private static final EnumMap<Result, AtomicInteger> results = new EnumMap<Result, AtomicInteger>(Result.class);

    enum Result {
        PASSING, FAIL_LEFT, FAIL_RIGHT, FAIL_BOTH;
    }

    /**
     * This iterates through million random curves to determine the
     * condition-number needed to estimate an upper-limit of the 
     * numerical accuracy.
     * @param args unused
     */
    public static void main(String[] args) {
        System.out.println(FindExtremaAccuracyTest.class.getSimpleName() + ": start");
        System.out.println("N: " + N);

        if (DEBUG) {
            // m 0.0 0.0 c 2.63209572625E8 0.3 -4.700001119893978E8 0.6 0.0 1.0 
            test(0, 0.0, 0.0, 2.63209572625E8, 0.3, -4.700001119893978E8, 0.6, 0.0, 1.0);

            test(0, 1e9);
            return;
        }

        // original:
        test(0, 10);
        test(0, 1e-6);
        test(0, 1e6);
        test(0, 1e9);

        test(1e6, 10);
        test(1e6, 1e-3);
        test(1e6, 1e6);
        test(1e6, 1e9);

        test(1e9, 10);
        test(1e9, 1e6);
        test(1e9, 1e9);

        for (double scale = 1e6; scale <= 1e16; scale *= 10.0) {
            test(scale, scale);
        }
        System.out.println(FindExtremaAccuracyTest.class.getSimpleName() + ": done");
    }

    private static void test(final double off, final double rng) {
        System.out.println("-------------------");
        System.out.println("offset value: " + off);
        System.out.println("random scale: " + rng);

        resetStats();

        // Initialize random once (seed fixed):
        final Random random = new Random(0L);

        for (int n = 0; n < N; n++) {
            test(n, random, off, rng);

            if (n % M == 0) {
                System.out.println("Iteration " + n + " ---");
                // dumpStats(false);
            }
        }
        System.out.println("Test done ---");
        dumpStats(true);
    }

    private static void dumpStats(boolean verbose) {
        System.out.println("maxCond:   " + condStats.max());
        System.out.println("ratioCond: " + ratioCondStats.max());

        if (verbose) {
            System.out.println("stats(cond) : " + condStats);
            System.out.println("stats(ratio): " + ratioCondStats);

            System.out.println("stats(results):");
            for (Result r : Result.values()) {
                System.out.println(r.toString() + ":\t" + results.get(r));
            }
        }
    }

    private static void resetStats() {
        condStats.reset();
        ratioCondStats.reset();

        results.clear();
        for (Result r : Result.values()) {
            results.put(r, new AtomicInteger());
        }
    }

    private static void test(int trial, Random random, final double off, final double rng) {
        final double half = rng / 2.0;
        double x0 = off;
        double y0 = 0.0;
        double cx1 = off + (random.nextDouble() * rng - half);
        double cy1 = 0.3;
        double cx2 = off + (random.nextDouble() * rng - half);
        double cy2 = 0.6;
        double x3 = off;
        double y3 = 1.0;

        test(trial, x0, y0, cx1, cy1, cx2, cy2, x3, y3);
    }

    private static void test(int trial,
                             double x0, double y0,
                             double cx1, double cy1,
                             double cx2, double cy2,
                             double x3, double y3) {

        final CubicCurve2D curve = new CubicCurve2D.Double(x0, y0, cx1, cy1, cx2, cy2, x3, y3);

        // The incoming data from a PathIterator is always represented by doubles, so that needs
        // to be where we start. (That is: if there's machine error already baked into those
        // doubles, then that's not something we can control for or accommodate.)
        // ... but everything that follows can, technically be calculated in really high precision:
        final BigDecimal[] coeff = new BigDecimal[4];
        final BigDecimal[] deriv_coeff = new BigDecimal[3];
        final BigDecimal[] tExtrema = new BigDecimal[2];
        final BigDecimal[] xExtrema = new BigDecimal[2];

        findExtrema(x0, cx1, cx2, x3, coeff, deriv_coeff, tExtrema, xExtrema);

        final BigDecimal leftX = xExtrema[0];
        final BigDecimal rightX = xExtrema[1];

        // TODO: extract findExtrema(double[]) to be more direct & efficient in long tests:
        final Rectangle2D r = getBounds2D(curve.getPathIterator(null));

        final BigDecimal observedLeftX = new BigDecimal(r.getMinX());
        final BigDecimal observedRightX = new BigDecimal(r.getMaxX());

        final boolean badMin = observedLeftX.compareTo(leftX) > 0;
        final boolean badMax = observedRightX.compareTo(rightX) < 0;

        final Result result;
        if (badMin && badMax) {
            result = Result.FAIL_BOTH;
        } else if (badMin) {
            result = Result.FAIL_LEFT;
        } else if (badMax) {
            result = Result.FAIL_RIGHT;
        } else {
            result = Result.PASSING;
        }
        // update result stats:
        results.get(result).incrementAndGet();

        /*
        Ideally condition number: cond(t) = sum (|ai|.|t^i|) / sum (ai.t^i)
        but as t in [0-1]: |t^i| = t^i = 1
         */
        final double cond = coeff[3].abs().add(coeff[2].abs()).add(coeff[1].abs()).add(coeff[0].abs()).doubleValue();
        final double ucond = Math.ulp(cond);

        final double deltaMin = leftX.subtract(observedLeftX).abs().doubleValue();
        final double deltaMax = rightX.subtract(observedRightX).abs().doubleValue();

        final double delta = Math.max((badMin ? deltaMin : 0.0), (badMax ? deltaMax : 0.0));
        final double uratio = delta / ucond;

        // update stats:
        condStats.add(cond);
        if (delta > 0.0) {
            ratioCondStats.add(uratio);
        }

        if (delta <= ucond) {
            return;
        }

        System.out.println("Cond number:\t" + cond + "\tulp: " + ucond);
        System.out.println("delta:\t" + delta);
        System.out.println("cond KO:\t" + (delta / ucond));

        System.out.println("Examining (trial #" + trial + "), " + result + ", " + toString(curve));

        System.out.println("Cond number?: " + cond + " ulp: " + ucond);
        System.out.println("deltaMin:\t" + deltaMin);
        System.out.println("deltaMax:\t" + deltaMax);

        System.out.println("cond KO: " + (delta / ucond));

        System.out.println("coeff[3]: " + coeff[3]);
        System.out.println("coeff[2]: " + coeff[2]);
        System.out.println("coeff[1]: " + coeff[1]);
        System.out.println("coeff[0]: " + coeff[0]);

        System.out.println("tExtrema: " + Arrays.toString(tExtrema));

        final String leftStr = toUniformString(leftX);
        final String rightStr = toUniformString(rightX);
        switch (result) {
            case FAIL_BOTH:
                System.out.println("Exp:\t" + leftStr + "\t" + rightStr);
                break;
            case FAIL_LEFT:
                System.out.println("Exp:\t" + leftStr);
                break;
            case FAIL_RIGHT:
                System.out.println("Exp:\t" + rightStr);
                break;
            default:
                break;
        }

        final String leftStr2 = toComparisonString(new BigDecimal(r.getMinX()), leftStr);
        final String rightStr2 = toComparisonString(new BigDecimal(r.getMaxX()), rightStr);

        switch (result) {
            case FAIL_BOTH:
                System.out.println("Orig:\t" + leftStr2 + "\t" + rightStr2);
                break;
            case FAIL_LEFT:
                System.out.println("Orig:\t" + leftStr2);
                break;
            case FAIL_RIGHT:
                System.out.println("Orig:\t" + rightStr2);
                break;
            default:
                break;
        }
    }

    /**
     * Convert a shape into SVG-ish notation for debugging/readability.
     * @param shape shape to convert
     * @return String (SVG-ish notation)
     */
    private static String toString(Shape shape) {
        final StringBuilder sb = new StringBuilder(256);
        final double[] coords = new double[6];

        for (final PathIterator pi = shape.getPathIterator(null); pi.isDone(); pi.next()) {
            final int type = pi.currentSegment(coords);

            switch (type) {
                case PathIterator.SEG_MOVETO:
                    sb.append("m ").append(coords[0]).append(" ").append(coords[1]).append(" ");
                    break;
                case PathIterator.SEG_LINETO:
                    sb.append("l ").append(coords[0]).append(" ").append(coords[1]).append(" ");
                    break;
                case PathIterator.SEG_QUADTO:
                    sb.append("q ").append(coords[0]).append(" ").append(coords[1]).append(" ")
                            .append(coords[2]).append(" ").append(coords[3]).append(" ");
                    break;
                case PathIterator.SEG_CUBICTO:
                    sb.append("c ").append(coords[0]).append(" ").append(coords[1]).append(" ")
                            .append(coords[2]).append(" ").append(coords[3])
                            .append(" ").append(coords[4]).append(" ").append(coords[5]).append(" ");
                    break;
                case PathIterator.SEG_CLOSE:
                    sb.append("z");
                    break;
                default:
                    break;
            }
        }
        return sb.toString();
    }

    private static String toUniformString(BigDecimal decimal) {
        int DIGIT_COUNT = 40;
        String str = decimal.toPlainString();
        if (str.length() >= DIGIT_COUNT) {
            str = str.substring(0, DIGIT_COUNT - 1) + "…";
        }
        while (str.length() < DIGIT_COUNT) {
            str += " ";
        }
        return str;
    }

    private static String toComparisonString(BigDecimal target, String compareAgainst) {
        final String str = toUniformString(target);

        for (int i = 0; i < str.length(); i++) {
            char ch1 = str.charAt(i);
            char ch2 = compareAgainst.charAt(i);
            if (ch1 != ch2) {
                return str.substring(0, i) + createCircleDigit(ch1) + str.substring(i + 1);
            }
        }
        return str;
    }

    /**
     * Convert a digit 0-9 into a "circle digit". Really we just want any unobtrusive way to
     * highlight a character.
     * @param ch char to convert
     * @return unicode character "circle digit"
     */
    private static char createCircleDigit(char ch) {
        if (ch >= '1' && ch <= '9') {
            return (char) (ch - '1' + '\u2460');
        }
        if (ch == '0') {
            return '\u24ea';
        }
        return ch;
    }

    private final static BigDecimal TWO = new BigDecimal(2.0);
    private final static BigDecimal THREE = new BigDecimal(3.0);

    private static void findExtrema(final double x0, final double x1, final double x2, final double x3,
                                    final BigDecimal[] coeff, final BigDecimal[] deriv_coeff,
                                    final BigDecimal[] tExtrema, final BigDecimal[] range) {

        /*
        final double dx32 = 3.0 * (coords[2] - coords[0]);
        final double dx21 = 3.0 * (coords[0] - lastX);
        x_coeff[3] = (coords[4] - lastX) - dx32;  // A = P3 - P0 - 3 (P2 - P1) = (P3 - P0) + 3 (P1 - P2)
        x_coeff[2] = (dx32 - dx21);               // B = 3 (P2 - P1) - 3(P1 - P0) = 3 (P2 + P0) - 6 P1
        x_coeff[1] = dx21;                        // C = 3 (P1 - P0)
        x_coeff[0] = lastX;                       // D = P0

        x_deriv_coeff[0] = x_coeff[1];
        x_deriv_coeff[1] = 2.0 * x_coeff[2];
        x_deriv_coeff[2] = 3.0 * x_coeff[3];
         */
        final BigDecimal bx0 = new BigDecimal(x0);
        final BigDecimal bx1 = new BigDecimal(x1);
        final BigDecimal bx2 = new BigDecimal(x2);
        final BigDecimal bx3 = new BigDecimal(x3);

        final BigDecimal dx21 = THREE.multiply(bx2.subtract(bx1));
        final BigDecimal dx10 = THREE.multiply(bx1.subtract(bx0));
        coeff[3] = bx3.subtract(bx0).subtract(dx21);
        coeff[2] = dx21.subtract(dx10);
        coeff[1] = dx10;
        coeff[0] = bx0;

        deriv_coeff[0] = coeff[1];
        deriv_coeff[1] = coeff[2].multiply(TWO);
        deriv_coeff[2] = coeff[3].multiply(THREE);

        final int tExtremaCount = solveQuadratic(deriv_coeff, tExtrema);

        BigDecimal leftX = bx0;
        BigDecimal rightX = bx0;

        for (int i = 0; i < tExtremaCount; i++) {
            final BigDecimal t = tExtrema[i];

            if (t.compareTo(BigDecimal.ZERO) > 0 && t.compareTo(BigDecimal.ONE) < 0) {
                final BigDecimal x = coeff[0].add(t.multiply(coeff[1].add(t.multiply(coeff[2].add(t.multiply(coeff[3]))))));
                if (x.compareTo(leftX) < 0) {
                    leftX = x;
                }
                if (x.compareTo(rightX) > 0) {
                    rightX = x;
                }
            }
        }
        range[0] = leftX;
        range[1] = rightX;
    }

    private static int solveQuadratic(final BigDecimal[] eqn, final BigDecimal[] res) {
        BigDecimal a = eqn[2];
        BigDecimal b = eqn[1];
        BigDecimal c = eqn[0];
        int roots = 0;
        if (a.equals(BigDecimal.ZERO)) {
            // The quadratic parabola has degenerated to a line.
            if (b.equals(BigDecimal.ZERO)) {
                // The line has degenerated to a constant.
                return -1;
            }
            res[roots++] = c.negate().divide(b, RoundingMode.HALF_EVEN);
        } else {
            // From Numerical Recipes, 5.6, Quadratic and Cubic Equations
            BigDecimal d = b.multiply(b).add(new BigDecimal(-4.0).multiply(a).multiply(c));
            if (d.compareTo(BigDecimal.ZERO) < 0) {
                // If d < 0.0, then there are no roots
                return 0;
            }
            // UNLIMITED precision can not be achieved => use DECIMAL128
            /*
            java.lang.ArithmeticException: Computed square root not exact.
        	at java.base/java.math.BigDecimal.sqrt(BigDecimal.java:2264)
             */
            d = d.sqrt(MathContext.DECIMAL128);
            // For accuracy, calculate one root using:
            //     (-b +/- d) / 2a
            // and the other using:
            //     2c / (-b +/- d)
            // Choose the sign of the +/- so that b+d gets larger in magnitude
            if (b.compareTo(BigDecimal.ZERO) < 0) {
                d = d.negate();
            }
            final BigDecimal q = b.add(d).divide(TWO.negate(), RoundingMode.HALF_EVEN);
            // We already tested a for being 0 above
            res[roots++] = q.divide(a, RoundingMode.HALF_EVEN);

            if (!q.equals(BigDecimal.ZERO)) {
                res[roots++] = c.divide(q, RoundingMode.HALF_EVEN);
            }
        }
        return roots;
    }

    /**
     * Returns a high precision bounding box of the specified PathIterator.
     * <p>
     * This method provides a basic facility for implementors of the {@link Shape} interface to
     * implement support for the {@link Shape#getBounds2D()} method.
     * </p>
     * @param pi path iterator
     * @return an instance of {@code Rectangle2D} that is a high-precision bounding box of the
     *         {@code PathIterator}.
     * @see Shape#getBounds2D()
     */
    public static Rectangle2D getBounds2D(final PathIterator pi) {
        // define x and y parametric coefficients where:
        // x(t) = x_coeff[0] + x_coeff[1] * t + x_coeff[2] * t^2 + x_coeff[3] * t^3
        final double[] x_coeff = new double[4];
        final double[] y_coeff = new double[4];

        // define the derivative's coefficients
        final double[] x_deriv_coeff = new double[3];
        final double[] y_deriv_coeff = new double[3];

        final double[] coords = new double[6];
        final double[] tExtrema = new double[2];
        boolean isDefined = false;
        double leftX = 0.0;
        double rightX = 0.0;
        double topY = 0.0;
        double bottomY = 0.0;
        double lastX = 0.0;
        double lastY = 0.0;

        for (; !pi.isDone(); pi.next()) {
            int type = pi.currentSegment(coords);
            switch (type) {
                case PathIterator.SEG_MOVETO:
                    if (!isDefined) {
                        isDefined = true;
                        leftX = rightX = coords[0];
                        topY = bottomY = coords[1];
                    } else {
                        if (coords[0] < leftX) {
                            leftX = coords[0];
                        }
                        if (coords[0] > rightX) {
                            rightX = coords[0];
                        }
                        if (coords[1] < topY) {
                            topY = coords[1];
                        }
                        if (coords[1] > bottomY) {
                            bottomY = coords[1];
                        }
                    }
                    lastX = coords[0];
                    lastY = coords[1];
                    break;
                case PathIterator.SEG_LINETO:
                    if (coords[0] < leftX) {
                        leftX = coords[0];
                    }
                    if (coords[0] > rightX) {
                        rightX = coords[0];
                    }
                    if (coords[1] < topY) {
                        topY = coords[1];
                    }
                    if (coords[1] > bottomY) {
                        bottomY = coords[1];
                    }
                    lastX = coords[0];
                    lastY = coords[1];
                    break;
                case PathIterator.SEG_QUADTO:
                    if (coords[2] < leftX) {
                        leftX = coords[2];
                    }
                    if (coords[2] > rightX) {
                        rightX = coords[2];
                    }
                    if (coords[3] < topY) {
                        topY = coords[3];
                    }
                    if (coords[3] > bottomY) {
                        bottomY = coords[3];
                    }

                    if (coords[0] < leftX || coords[0] > rightX) {
                        final double dx21 = (coords[0] - lastX);
                        x_coeff[2] = (coords[2] - coords[0]) - dx21;  // A = P3 - P0 - 2 P2
                        x_coeff[1] = 2.0 * dx21;                      // B = 2 (P2 - P1)
                        x_coeff[0] = lastX;                           // C = P1

                        x_deriv_coeff[0] = x_coeff[1];
                        x_deriv_coeff[1] = 2.0 * x_coeff[2];

                        final double t = -x_deriv_coeff[0] / x_deriv_coeff[1];

                        if (t > 0.0 && t < 1.0) {
                            double x = x_coeff[0] + t * (x_coeff[1] + t * x_coeff[2]);
                            if (x < leftX) {
                                leftX = x;
                            }
                            if (x > rightX) {
                                rightX = x;
                            }
                        }
                    }
                    if (coords[1] < topY || coords[1] > bottomY) {
                        final double dy21 = (coords[1] - lastY);
                        y_coeff[2] = (coords[3] - coords[1]) - dy21;
                        y_coeff[1] = 2.0 * dy21;
                        y_coeff[0] = lastY;

                        y_deriv_coeff[0] = y_coeff[1];
                        y_deriv_coeff[1] = 2.0 * y_coeff[2];

                        final double t = -y_deriv_coeff[0] / y_deriv_coeff[1];

                        if (t > 0.0 && t < 1.0) {
                            double y = y_coeff[0] + t * (y_coeff[1] + t * y_coeff[2]);
                            if (y < topY) {
                                topY = y;
                            }
                            if (y > bottomY) {
                                bottomY = y;
                            }
                        }
                    }
                    lastX = coords[2];
                    lastY = coords[3];
                    break;
                case PathIterator.SEG_CUBICTO:
                    if (coords[4] < leftX) {
                        leftX = coords[4];
                    }
                    if (coords[4] > rightX) {
                        rightX = coords[4];
                    }
                    if (coords[5] < topY) {
                        topY = coords[5];
                    }
                    if (coords[5] > bottomY) {
                        bottomY = coords[5];
                    }

                    if (coords[0] < leftX || coords[0] > rightX || coords[2] < leftX || coords[2] > rightX) {
                        final double dx32 = 3.0 * (coords[2] - coords[0]);
                        final double dx21 = 3.0 * (coords[0] - lastX);
                        x_coeff[3] = (coords[4] - lastX) - dx32;  // A = P3 - P0 - 3 (P2 - P1) = (P3 - P0) + 3 (P1 - P2)
                        x_coeff[2] = (dx32 - dx21);               // B = 3 (P2 - P1) - 3(P1 - P0) = 3 (P2 + P0) - 6 P1
                        x_coeff[1] = dx21;                        // C = 3 (P1 - P0)
                        x_coeff[0] = lastX;                       // D = P0

                        x_deriv_coeff[0] = x_coeff[1];
                        x_deriv_coeff[1] = 2.0 * x_coeff[2];
                        x_deriv_coeff[2] = 3.0 * x_coeff[3];

                        final int tExtremaCount = QuadCurve2D.solveQuadratic(x_deriv_coeff, tExtrema);

                        for (int i = 0; i < tExtremaCount; i++) {
                            final double t = tExtrema[i];
                            if (t > 0.0 && t < 1.0) {
                                final double x = x_coeff[0] + t * (x_coeff[1] + t * (x_coeff[2] + t * x_coeff[3]));
                                if (x < leftX) {
                                    leftX = x;
                                }
                                if (x > rightX) {
                                    rightX = x;
                                }
                            }
                        }
                    }
                    if (coords[1] < topY || coords[1] > bottomY || coords[3] < topY || coords[3] > bottomY) {
                        final double dy32 = 3.0 * (coords[3] - coords[1]);
                        final double dy21 = 3.0 * (coords[1] - lastY);
                        y_coeff[3] = (coords[5] - lastY) - dy32;
                        y_coeff[2] = (dy32 - dy21);
                        y_coeff[1] = dy21;
                        y_coeff[0] = lastY;

                        y_deriv_coeff[0] = y_coeff[1];
                        y_deriv_coeff[1] = 2.0 * y_coeff[2];
                        y_deriv_coeff[2] = 3.0 * y_coeff[3];

                        final int tExtremaCount = QuadCurve2D.solveQuadratic(y_deriv_coeff, tExtrema);

                        for (int i = 0; i < tExtremaCount; i++) {
                            final double t = tExtrema[i];
                            if (t > 0.0 && t < 1.0) {
                                final double y = y_coeff[0] + t * (y_coeff[1] + t * (y_coeff[2] + t * y_coeff[3]));
                                if (y < topY) {
                                    topY = y;
                                }
                                if (y > bottomY) {
                                    bottomY = y;
                                }
                            }
                        }
                    }
                    lastX = coords[4];
                    lastY = coords[5];
                    break;
                case PathIterator.SEG_CLOSE:
                default:
            }
        }
        if (isDefined) {
            return new Rectangle2D.Double(leftX, topY, rightX - leftX, bottomY - topY);
        }

        // there's room to debate what should happen here, but historically we return a zeroed
        // out rectangle here. So for backwards compatibility let's keep doing that:
        return new Rectangle2D.Double();
    }
}