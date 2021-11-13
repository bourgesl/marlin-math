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
import static java.math.BigDecimal.ONE;
import static java.math.BigDecimal.ZERO;
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
public class FindExtremaAccuracyTest extends BaseTest {

    private final static CurveType TEST_TYPE = CurveType.CUBIC;
    private final static TestDiff TEST_MODE = TestDiff.EDGE_OR; // DIST

    private final static boolean DEBUG = false;
    private final static boolean TRACE = false;

    private final static boolean DEBUG_QUAD_SOLVER = false;

    private final static boolean USE_PATH_ITER = false;

    private final static int N = 1000000;
    private final static int M = N / 10;

    /** statistics on condition-number */
    private static final WelfordVariance condStats = new WelfordVariance();
    /** statistics on ratio(delta / ucond) */
    private static final WelfordVariance ratioCondStats = new WelfordVariance();
    /** statistics on n * ulp(x) */
    private static final WelfordVariance ulpXStats = new WelfordVariance();

    private static final EnumMap<Result, AtomicInteger> results = new EnumMap<Result, AtomicInteger>(Result.class);

    enum CurveType {
        QUAD, CUBIC;
    }

    enum TestDiff {
        DIST, EDGE_OR;
    }

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
        System.out.println("-------------------");
        System.out.println(FindExtremaAccuracyTest.class.getSimpleName() + ": start");
        System.out.println("USE_PATH_ITER: " + USE_PATH_ITER);
        System.out.println("N: " + N);

        CurveType ctype = TEST_TYPE;
        System.out.println("TYPE: " + ctype);
        System.out.println("TEST_MODE: " + TEST_MODE);

        if (DEBUG) {
            resetStats();

            // Bad cubic:
            // t is inacurrate != 0.161:
            /*
            tExtrema:   [0.7454713624448268522743715, 0.161]
            tExtrema_d: [0.7454713624448269, 0.16051708637880857]
            t: 0.7454713624448269 delta(t): -7.152438593426965E-17 nUlp: -0.6442343956829607
            t: 0.16051708637880857 delta(t): 4.829136211914289E-4 nUlp: 1.7398796835598848E13
             */
            test(CurveType.CUBIC, 0, 0.0, 0.0, 2.63209572625E8, 0.3, -4.700001119893978E8, 0.6, 0.0, 1.0);

            if (false) {
                test(ctype, 0, 1e9);
            }
            return;
        }

        // original:
        test(ctype, 0, 10);
        test(ctype, 0, 1e-6);
        test(ctype, 0, 1e6);
        test(ctype, 0, 1e9);

        test(ctype, 1e6, 10);
        test(ctype, 1e6, 1e-3);
        test(ctype, 1e6, 1e6);
        test(ctype, 1e6, 1e9);

        test(ctype, 1e9, 10);
        test(ctype, 1e9, 1e6);
        test(ctype, 1e9, 1e9);

        // test scaling:
        for (double scale = 1e6; scale <= 1e16; scale *= 10.0) {
            test(ctype, scale, scale);
        }
        System.out.println(FindExtremaAccuracyTest.class.getSimpleName() + ": done");
        System.out.println("-------------------");
    }

    private static void test(CurveType ctype, final double off, final double rng) {
        System.out.println("-------------------");
        System.out.println("offset value: " + off);
        System.out.println("random scale: " + rng);

        resetStats();

        // Initialize random once (seed fixed):
        final Random random = new Random(0L);

        for (int n = 0; n < N; n++) {
            test(ctype, n, random, off, rng);

            if (n % M == 0) {
                System.out.println("Iteration " + n + " ---");
                // dumpStats(false);
            }
        }
        System.out.println("Test done ---");
        dumpStats(true);
    }

    private static void dumpStats(boolean verbose) {
        System.out.println("maxCond:       " + condStats.max());
        final double ucondMax = Math.ulp(condStats.max());
        System.out.println("ulp(maxCond):  " + ucondMax);
        System.out.println("ratioCond:     " + ratioCondStats.max());

        if (verbose) {
            System.out.println("stats(cond) :  " + condStats);
            System.out.println("stats(ratio):  " + ratioCondStats);
            System.out.println("stats(ulp(x)): " + ulpXStats);

            System.out.println("result stats:");
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

    private static void test(CurveType ctype, int trial, Random random, final double off, final double rng) {
        final double half = rng / 2.0;
        double x0 = off;
        double cx1 = off + (random.nextDouble() * rng - half);
        double cx2 = off + (random.nextDouble() * rng - half);
        double x3 = off;

        double y0 = 0.0;
        double cy1 = 0.3;
        double cy2 = 0.6;
        double y3 = 1.0;

        test(ctype, trial, x0, y0, cx1, cy1, cx2, cy2, x3, y3);
    }

    private static void test(CurveType ctype, int trial,
                             double x0, double y0,
                             double cx1, double cy1,
                             double cx2, double cy2,
                             double x3, double y3) {

        // The incoming data from a PathIterator is always represented by doubles, so that needs
        // to be where we start. (That is: if there's machine error already baked into those
        // doubles, then that's not something we can control for or accommodate.)
        // ... but everything that follows can, technically be calculated in really high precision:
        final BigDecimal[] coeff = new BigDecimal[4];
        final BigDecimal[] deriv_coeff = new BigDecimal[3];
        final BigDecimal[] tExtrema = new BigDecimal[2];
        final BigDecimal[] xExtrema = new BigDecimal[2];

        // BigDecimal high-accuracy implementation (ref):
        switch (ctype) {
            case QUAD:
                findExtrema(x0, cx1, x3, coeff, deriv_coeff, tExtrema, xExtrema);
                break;
            default:
            case CUBIC:
                findExtrema(x0, cx1, cx2, x3, coeff, deriv_coeff, tExtrema, xExtrema);
                break;
        }

        final BigDecimal minX = xExtrema[0];
        final BigDecimal maxX = xExtrema[1];

        final double[] coeff_d = new double[4];
        final double[] deriv_coeff_d = new double[3];
        final double[] tExtrema_d = new double[2];
        final double[] xExtrema_d = new double[2];

        final double minX_d;
        final double maxX_d;

        if (!USE_PATH_ITER) {
            switch (ctype) {
                case QUAD:
                    findExtrema_d(x0, cx1, x3, coeff_d, deriv_coeff_d, tExtrema_d, xExtrema_d);
                    break;
                default:
                case CUBIC:
                    findExtrema_d(x0, cx1, cx2, x3, coeff_d, deriv_coeff_d, tExtrema_d, xExtrema_d);
                    break;
            }
            minX_d = xExtrema_d[0];
            maxX_d = xExtrema_d[1];
        } else {
            final Shape curve;
            switch (ctype) {
                case QUAD:
                    curve = new QuadCurve2D.Double(x0, y0, cx1, cy1, x3, y3);
                    break;
                default:
                case CUBIC:
                    curve = new CubicCurve2D.Double(x0, y0, cx1, cy1, cx2, cy2, x3, y3);
                    break;
            }
            final Rectangle2D r = getBounds2D(curve.getPathIterator(null));
            minX_d = r.getMinX();
            maxX_d = r.getMaxX();
        }

        final BigDecimal obsMinX = new BigDecimal(minX_d);
        final BigDecimal obsMaxX = new BigDecimal(maxX_d);

        final boolean badMin = obsMinX.compareTo(minX) > 0;
        final boolean badMax = obsMaxX.compareTo(maxX) < 0;

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
        final double cond;
        switch (ctype) {
            case QUAD:
                cond = coeff[2].abs().add(coeff[1].abs()).add(coeff[0].abs()).doubleValue();
                break;
            default:
            case CUBIC:
                cond = coeff[3].abs().add(coeff[2].abs()).add(coeff[1].abs()).add(coeff[0].abs()).doubleValue();
                break;
        }
        final double ucond = Math.ulp(cond);

        final double deltaMin = minX.subtract(obsMinX).abs().doubleValue();
        final double deltaMax = maxX.subtract(obsMaxX).abs().doubleValue();

        final double delta;
        switch (TEST_MODE) {
            case EDGE_OR:
                delta = Math.max((badMin ? deltaMin : 0.0), (badMax ? deltaMax : 0.0));
                break;
            default:
            case DIST:
                delta = Math.max(deltaMin, deltaMax);
        }

        // update stats:
        condStats.add(cond);
        if (delta > 0.0) {
            final double uratio = delta / ucond;

            ratioCondStats.add(uratio);

            final double nUlpMin = deltaMin / Math.ulp(minX_d);
            final double nUlpMax = deltaMax / Math.ulp(maxX_d);

            final double nUlp;
            switch (TEST_MODE) {
                case EDGE_OR:
                    nUlp = Math.max((badMin ? nUlpMin : 0.0), (badMax ? nUlpMax : 0.0));
                    break;
                default:
                case DIST:
                    nUlp = Math.max(nUlpMin, nUlpMax);
            }
            ulpXStats.add(nUlp);
        }

        if (!TRACE && (delta <= ucond)) {
            // test OK
            return;
        }

        System.out.println("Examining (trial #" + trial + "), " + result + ", "
                + (ctype == CurveType.QUAD ? toString(x0, y0, cx1, cy1, 0.0, 0.0, x3, y3)
                        : toString(x0, y0, cx1, cy1, cx2, cy2, x3, y3))
        );

        System.out.println("Cond number:\t" + cond + "\tulp: " + ucond);
        System.out.println("deltaMin:\t" + deltaMin);
        System.out.println("deltaMax:\t" + deltaMax);
        System.out.println("delta:\t" + delta);
        System.out.println("cond KO:\t" + (delta / ucond));

        System.out.println("coeff[3]: " + coeff[3]);
        System.out.println("coeff[2]: " + coeff[2]);
        System.out.println("coeff[1]: " + coeff[1]);
        System.out.println("coeff[0]: " + coeff[0]);

        System.out.println("dcoeff[2]: " + deriv_coeff[2]);
        System.out.println("dcoeff[1]: " + deriv_coeff[1]);
        System.out.println("dcoeff[0]: " + deriv_coeff[0]);

        System.out.println("tExtrema:   " + Arrays.toString(tExtrema));
        System.out.println("tExtrema_d: " + Arrays.toString(tExtrema_d));

        for (int i = 0; i < 2; i++) {
            final BigDecimal t = tExtrema[i];
            if (t != null) {
                final double deltaT = t.subtract(new BigDecimal(tExtrema_d[i])).doubleValue();
                final double nUlpT = deltaT / Math.ulp(tExtrema_d[i]);
                System.out.println("t: " + tExtrema_d[i] + " delta(t): " + deltaT + " nUlp: " + nUlpT);
            }
        }

        final String leftStr = toUniformString(minX);
        final String rightStr = toUniformString(maxX);
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
            case PASSING:
                System.out.println("Exp:\t" + leftStr + "\t" + rightStr);
                break;
            default:
                break;
        }

        final String leftStr2 = toComparisonString(obsMinX, leftStr);
        final String rightStr2 = toComparisonString(obsMaxX, rightStr);

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
            case PASSING:
                System.out.println("Orig:\t" + leftStr2 + "\t" + rightStr2);
                break;
            default:
                break;
        }
    }

    // Find Extrema implementations (High-Accuracy) based on BigDecimal
    private final static BigDecimal TWO = new BigDecimal(2.0);
    private final static BigDecimal THREE = new BigDecimal(3.0);
    private final static BigDecimal HALF_NEG = new BigDecimal(-0.5);

    private static void findExtrema(final double x1, final double x2, final double x3,
                                    final BigDecimal[] coeff, final BigDecimal[] deriv_coeff,
                                    final BigDecimal[] tExtrema, final BigDecimal[] range) {
        // Quad
        final BigDecimal bx1 = new BigDecimal(x1);
        final BigDecimal bx2 = new BigDecimal(x2);
        final BigDecimal bx3 = new BigDecimal(x3);

        final BigDecimal dx21 = bx2.subtract(bx1);
        coeff[3] = ZERO; // useless
        coeff[2] = bx3.subtract(bx2).subtract(dx21);  // A =   (P3 - P2) - (P2 - P1)
        coeff[1] = dx21.multiply(TWO);                // B = 2 (P2 - P1)
        coeff[0] = bx1;                               // C = P1

        deriv_coeff[2] = ZERO; // useless
        deriv_coeff[1] = coeff[2].multiply(TWO);
        deriv_coeff[0] = coeff[1];

//        final int tExtremaCount = QuadCurve2D.solveQuadratic(deriv_coeff, tExtrema);
        final BigDecimal t;

        // The quadratic parabola has degenerated to a line.
        if (deriv_coeff[1].equals(ZERO)) {
            // The line has degenerated to a constant.
            t = null;
        } else {
            t = deriv_coeff[0].negate().divide(deriv_coeff[1]);
        }
        tExtrema[0] = t;
        tExtrema[1] = null;

        BigDecimal minX = bx1;
        BigDecimal maxX = bx1;

        if ((t != null) && (t.compareTo(ZERO) > 0) && (t.compareTo(ONE) < 0)) {
            final BigDecimal x = coeff[0].add(t.multiply(coeff[1].add(t.multiply(coeff[2]))));
            if (x.compareTo(minX) < 0) {
                minX = x;
            }
            if (x.compareTo(maxX) > 0) {
                maxX = x;
            }
        }
        range[0] = minX;
        range[1] = maxX;
    }

    private static void findExtrema(final double x1, final double x2, final double x3, final double x4,
                                    final BigDecimal[] coeff, final BigDecimal[] deriv_coeff,
                                    final BigDecimal[] tExtrema, final BigDecimal[] range) {
        // Cubic
        final BigDecimal bx1 = new BigDecimal(x1);
        final BigDecimal bx2 = new BigDecimal(x2);
        final BigDecimal bx3 = new BigDecimal(x3);
        final BigDecimal bx4 = new BigDecimal(x4);

        final BigDecimal dx32 = THREE.multiply(bx3.subtract(bx2));
        final BigDecimal dx21 = THREE.multiply(bx2.subtract(bx1));
        coeff[3] = bx4.subtract(bx1).subtract(dx32);  // A =   (P4 - P1) - 3 (P3 - P2)
        coeff[2] = dx32.subtract(dx21);               // B = 3 (P3 - P2) - 3 (P2 - P1)
        coeff[1] = dx21;                              // C = 3 (P2 - P1)
        coeff[0] = bx1;                               // D = P1

        deriv_coeff[2] = coeff[3].multiply(THREE);
        deriv_coeff[1] = coeff[2].multiply(TWO);
        deriv_coeff[0] = coeff[1];

        tExtrema[0] = null;
        tExtrema[1] = null;
        final int tExtremaCount = solveQuadratic(deriv_coeff, tExtrema);

        BigDecimal minX = bx1;
        BigDecimal maxX = bx1;

        for (int i = 0; i < tExtremaCount; i++) {
            final BigDecimal t = tExtrema[i];

            if ((t.compareTo(ZERO) > 0) && (t.compareTo(ONE) < 0)) {
                final BigDecimal x = coeff[0].add(t.multiply(coeff[1].add(t.multiply(coeff[2].add(t.multiply(coeff[3]))))));
                if (x.compareTo(minX) < 0) {
                    minX = x;
                }
                if (x.compareTo(maxX) > 0) {
                    maxX = x;
                }
            }
        }
        range[0] = minX;
        range[1] = maxX;
    }

    // Derived from QuadCurve2D.solveQuadratic
    private static int solveQuadratic(final BigDecimal[] eqn, final BigDecimal[] res) {
        final BigDecimal a = eqn[2];
        final BigDecimal b = eqn[1];
        final BigDecimal c = eqn[0];

        int roots = 0;
        if (a.equals(ZERO)) {
            // The quadratic parabola has degenerated to a line.
            if (b.equals(ZERO)) {
                // The line has degenerated to a constant.
                return -1;
            }
            res[roots++] = c.negate().divide(b);
        } else {
            // From Numerical Recipes, 5.6, Quadratic and Cubic Equations
            BigDecimal d = b.multiply(b).add(new BigDecimal(-4.0).multiply(a).multiply(c));
            if (DEBUG_QUAD_SOLVER) {
                System.out.println("d: " + d);
            }
            if (d.compareTo(ZERO) < 0) {
                // If d < 0.0, then there are no roots
                return 0;
            }
            // sqrt: UNLIMITED precision can not be achieved => use DECIMAL128
            /*
            java.lang.ArithmeticException: Computed square root not exact.
        	at java.base/java.math.BigDecimal.sqrt(BigDecimal.java:2264)
             */
            d = d.sqrt(MathContext.DECIMAL128);
            if (DEBUG_QUAD_SOLVER) {
                System.out.println("sqrt(d): " + d);
                System.out.println("b:       " + b);
            }
            // For accuracy, calculate one root using:
            //     (-b +/- d) / 2a
            // and the other using:
            //     2c / (-b +/- d)
            // Choose the sign of the +/- so that b+d gets larger in magnitude
            if (b.compareTo(ZERO) < 0) {
                d = d.negate();
            }
            final BigDecimal q = HALF_NEG.multiply(b.add(d)); // no div
            // We already tested a for being 0 above
            res[roots++] = q.divide(a, RoundingMode.HALF_EVEN);

            if (!q.equals(ZERO)) {
                res[roots++] = c.divide(q, RoundingMode.HALF_EVEN);
            }
        }
        return roots;
    }

    // Find Extrema implementations (Low-Accuracy) based on double primitive type (64bits)
    private static void findExtrema_d(final double x1, final double x2, final double x3,
                                      final double[] coeff, final double[] deriv_coeff,
                                      final double[] tExtrema, final double[] range) {
        // Quad
        final double dx21 = (x2 - x1);
        coeff[3] = 0.0; // useless 
        coeff[2] = (x3 - x2) - dx21;  // A =   (P3 - P2) - (P2 - P1)
        coeff[1] = 2.0 * dx21;        // B = 2 (P2 - P1)
        coeff[0] = x1;                // C = P1

        deriv_coeff[2] = 0.0; // useless
        deriv_coeff[1] = 2.0 * coeff[2];
        deriv_coeff[0] = coeff[1];

//        final int tExtremaCount = QuadCurve2D.solveQuadratic(deriv_coeff, tExtrema);
        final double t;

        // The quadratic parabola has degenerated to a line.
        if (deriv_coeff[1] == 0.0) {
            // The line has degenerated to a constant.
            t = -1.0;
        } else {
            t = -deriv_coeff[0] / deriv_coeff[1];
        }
        tExtrema[0] = t;
        tExtrema[1] = Double.NaN;

        double minX = x1;
        double maxX = x1;

        if ((t > 0.0) && (t < 1.0)) {
            final double x = coeff[0] + t * (coeff[1] + t * (coeff[2]));
            if (x < minX) {
                minX = x;
            }
            if (x > maxX) {
                maxX = x;
            }
        }
        range[0] = minX;
        range[1] = maxX;
    }

    private static void findExtrema_d(final double x1, final double x2, final double x3, final double x4,
                                      final double[] coeff, final double[] deriv_coeff,
                                      final double[] tExtrema, final double[] range) {
        // Cubic
        final double dx32 = 3.0 * (x3 - x2);
        final double dx21 = 3.0 * (x2 - x1);
        coeff[3] = (x4 - x1) - dx32;  // A =   (P4 - P1) - 3 (P3 - P2)
        coeff[2] = (dx32 - dx21);     // B = 3 (P3 - P2) - 3 (P2 - P1)
        coeff[1] = dx21;              // C = 3 (P2 - P1)
        coeff[0] = x1;                // D = P1

        deriv_coeff[2] = 3.0 * coeff[3];
        deriv_coeff[1] = 2.0 * coeff[2];
        deriv_coeff[0] = coeff[1];

        tExtrema[0] = Double.NaN;
        tExtrema[1] = Double.NaN;
        final int tExtremaCount = solveQuadratic(deriv_coeff, tExtrema);

        double minX = x1;
        double maxX = x1;

        for (int i = 0; i < tExtremaCount; i++) {
            final double t = tExtrema[i];

            if ((t > 0.0) && (t < 1.0)) {
                final double x = coeff[0] + t * (coeff[1] + t * (coeff[2] + t * coeff[3]));
                if (x < minX) {
                    minX = x;
                }
                if (x > maxX) {
                    maxX = x;
                }
            }
        }
        range[0] = minX;
        range[1] = maxX;
    }

    // Copied from QuadCurve2D.solveQuadratic
    public static int solveQuadratic(double[] eqn, double[] res) {
        double a = eqn[2];
        double b = eqn[1];
        double c = eqn[0];

        int roots = 0;
        if (a == 0.0) {
            // The quadratic parabola has degenerated to a line.
            if (b == 0.0) {
                // The line has degenerated to a constant.
                return -1;
            }
            res[roots++] = -c / b;
        } else {
            // From Numerical Recipes, 5.6, Quadratic and Cubic Equations
            double d = b * b - 4.0 * a * c;
            if (DEBUG_QUAD_SOLVER) {
                System.out.println("d:      " + d);
                d += Math.ulp(d);
                System.out.println("d+1ulp: " + d);
            }
            /*
            double d_alt = diff_of_products(b, b, a * 4.0, c);
            double d_alt2 = diff_of_products(b, b, a, c * 4.0);

            System.out.println("d_alt : " + d_alt);
            System.out.println("d_alt2: " + d_alt2);
             */
            if (d < 0.0) {
                // If d < 0.0, then there are no roots
                return 0;
            }
            d = Math.sqrt(d);
            if (DEBUG_QUAD_SOLVER) {
                System.out.println("sqrt(d): " + d);
                System.out.println("b:       " + b);
            }
            // For accuracy, calculate one root using:
            //     (-b +/- d) / 2a
            // and the other using:
            //     2c / (-b +/- d)
            // Choose the sign of the +/- so that b+d gets larger in magnitude
            if (b < 0.0) {
                d = -d;
            }

            // double q = (b + d) / -2.0;
            double q = -0.5 * (b + d); // no div
            /*
            double q = (b + d) / -2.0;
            t: 0.7454713624448269 delta(t): -7.152438643426965E-17 nUlp: -0.6442344001865603
            t: 0.16051708637880857 delta(t): 4.829136211914289E-4 nUlp: 1.7398796835598848E

            double q = -0.5 * (b + d); // LBO
            t: 0.7454713624448269 delta(t): -7.152438643426965E-17 nUlp: -0.6442344001865603
            t: 0.16051708637880857 delta(t): 4.829136211914289E-4 nUlp: 1.7398796835598848
             */
            // We already tested a for being 0 above
            res[roots++] = q / a;
            if (q != 0.0) {
                res[roots++] = c / q;
            }
        }
        return roots;
    }

    /*
     * diff_of_products() computes a*b-c*d with a maximum error <= 1.5 ulp
     * 
     * Claude-Pierre Jeannerod, Nicolas Louvet, and Jean-Michel Muller, 
     * "Further Analysis of Kahan's Algorithm for the Accurate Computation 
     * of 2x2 Determinants". Mathematics of Computation, Vol. 82, No. 284, 
     * Oct. 2013, pp. 2245-2264
     */
    static double diff_of_products(double a, double b, double c, double d) {
        double w = d * c;
        double e = Math.fma(-d, c, w);
        double f = Math.fma(a, b, -w);
        return f + e;
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

    private FindExtremaAccuracyTest() {
        super();
    }
}
