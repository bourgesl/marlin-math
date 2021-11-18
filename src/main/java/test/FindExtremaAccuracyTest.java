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
    private final static TestDiff TEST_MODE = TestDiff.DIST; // EDGE_OR or DIST

    private final static boolean DEBUG = false;

    private final static boolean DEBUG_QUAD_SOLVER = false;

    private final static boolean QUIET_RUN = true & !DEBUG;

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
        System.out.println("N:         " + N);

        CurveType ctype = TEST_TYPE;
        System.out.println("TYPE:      " + ctype);
        System.out.println("TEST_MODE: " + TEST_MODE);

        System.out.println("QUIET_RUN: " + QUIET_RUN);
        System.out.println("DEBUG:     " + DEBUG);

        final FindExtremaContext ctx = new FindExtremaContext();

        if (DEBUG) {
            resetStats();

            // Bad cubic:
            // t is inacurrate != 0.161:
            test(CurveType.CUBIC, 0, ctx, 0.0, 2.63209572625E8, -4.700001119893978E8, 0.0);

            if (false) {
                test(ctype, ctx, 0, 1e9);
            }
            return;
        }

        // original:
        test(ctype, ctx, 0, 10);
        test(ctype, ctx, 0, 1e-6);
        test(ctype, ctx, 0, 1e6);
        test(ctype, ctx, 0, 1e9);

        test(ctype, ctx, 1e6, 10);
        test(ctype, ctx, 1e6, 1e-3);
        test(ctype, ctx, 1e6, 1e6);
        test(ctype, ctx, 1e6, 1e9);

        test(ctype, ctx, 1e9, 10);
        test(ctype, ctx, 1e9, 1e6);
        test(ctype, ctx, 1e9, 1e9);

        // test scaling:
        for (double scale = 1e6; scale <= 1e16; scale *= 10.0) {
            test(ctype, ctx, scale, scale);
        }
        System.out.println(FindExtremaAccuracyTest.class.getSimpleName() + ": done");
        System.out.println("-------------------");
    }

    private static void test(CurveType ctype, final FindExtremaContext ctx,
                             final double off, final double rng) {
        System.out.println("-------------------");
        System.out.println("offset value: " + off);
        System.out.println("random scale: " + rng);

        resetStats();

        // Initialize random once (seed fixed):
        final Random random = new Random(0L);

        for (int n = 0; n < N; n++) {
            test(ctype, ctx, n, random, off, rng);

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

    private static void test(CurveType ctype, final FindExtremaContext ctx, int trial,
                             Random random, final double off, final double rng) {
        final double half = rng / 2.0;
        double x0 = off;
        double cx1 = off + (random.nextDouble() * rng - half);
        double cx2 = off + (random.nextDouble() * rng - half);
        double x3 = off;

        test(ctype, trial, ctx, x0, cx1, cx2, x3);
    }

    private static void test(final CurveType ctype, final int trial,
                             final FindExtremaContext ctx,
                             final double x1, final double x2, final double x3, final double x4) {

        // The incoming data from a PathIterator is always represented by doubles, so that needs
        // to be where we start. (That is: if there's machine error already baked into those
        // doubles, then that's not something we can control for or accommodate.)
        // ... but everything that follows can, technically be calculated in really high precision:
        // BigDecimal high-accuracy implementation (ref):
        switch (ctype) {
            case QUAD:
                findExtrema_d(ctx, x1, x2, x4);

                findExtrema(ctx, x1, x2, x4);
                break;
            default:
            case CUBIC:
                findExtrema_d(ctx, x1, x2, x3, x4);

                if (DEBUG) {
                    findExtremaWithDblCoeffs(ctx, x1, x2, x3, x4);
                } else {
                    findExtrema(ctx, x1, x2, x3, x4);
                }
                break;
        }

        // Get results through the context:
        final BigDecimal[] coeff = ctx.coeff;
        final BigDecimal[] deriv_coeff = ctx.deriv_coeff;
        final BigDecimal[] tExtrema = ctx.tExtrema;
        final BigDecimal[] xExtrema = ctx.xExtrema;

        final double[] coeff_d = ctx.coeff_d;
        final double[] deriv_coeff_d = ctx.deriv_coeff_d;
        final double[] tExtrema_d = ctx.tExtrema_d;
        final double[] xExtrema_d = ctx.xExtrema_d;

        final BigDecimal minX = xExtrema[0];
        final BigDecimal maxX = xExtrema[1];

        final double minX_d = xExtrema_d[0];
        final double maxX_d = xExtrema_d[1];

        // Compare results:
        final BigDecimal obsMinX = new BigDecimal(minX_d);
        final BigDecimal obsMaxX = new BigDecimal(maxX_d);

        final boolean badMin = obsMinX.compareTo(minX) > 0;
        final boolean badMax = obsMaxX.compareTo(maxX) < 0;

        final Result result;
        if (DEBUG) {
            result = Result.FAIL_BOTH;
        } else {
            if (badMin && badMax) {
                result = Result.FAIL_BOTH;
            } else if (badMin) {
                result = Result.FAIL_LEFT;
            } else if (badMax) {
                result = Result.FAIL_RIGHT;
            } else {
                result = Result.PASSING;
            }
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

        if (QUIET_RUN || (!DEBUG && (delta <= ucond))) {
            // test OK
            return;
        }

        System.out.println("Examining (trial #" + trial + "), " + result + ", "
                + (ctype == CurveType.QUAD ? toString(x1, x2, 0.0, x4) : toString(x1, x2, x3, x4))
        );

        System.out.println("Cond number:\t" + cond + "\tulp: " + ucond);
        System.out.println("deltaMin:\t" + deltaMin);
        System.out.println("deltaMax:\t" + deltaMax);
        System.out.println("delta:\t" + delta);
        System.out.println("Cond[" + ((delta <= ucond) ? "OK" : " KO") + "]:\t" + (delta / ucond));

        // Check parameters:
        System.out.println("----");
        System.out.println("coeffs:   " + Arrays.toString(coeff));
        System.out.println("coeffs_d: " + Arrays.toString(coeff_d));

        for (int i = 0; i < 4; i++) {
            final BigDecimal bd = new BigDecimal(coeff_d[i]);
            final String bStr = toUniformString(coeff[i]);
            final String dStr = toComparisonString(bd, bStr);
            System.out.println("coeffs[" + i + "]: Exp: " + bStr);
            System.out.println("coeffs[" + i + "]: Dbl: " + dStr);

            final double diff = coeff[i].subtract(bd).doubleValue();
            final double nUlp = diff / Math.ulp(coeff_d[i]);
            System.out.println("delta: " + diff + " nUlp: " + nUlp);
        }

        System.out.println("dcoeff:   " + Arrays.toString(deriv_coeff));
        System.out.println("dcoeff_d: " + Arrays.toString(deriv_coeff_d));

        for (int i = 0; i < 3; i++) {
            final BigDecimal bd = new BigDecimal(deriv_coeff_d[i]);
            final String bStr = toUniformString(deriv_coeff[i]);
            final String dStr = toComparisonString(bd, bStr);
            System.out.println("dcoeff[" + i + "]: Exp: " + bStr);
            System.out.println("dcoeff[" + i + "]: Dbl: " + dStr);

            final double diff = deriv_coeff[i].subtract(bd).doubleValue();
            final double nUlp = diff / Math.ulp(deriv_coeff_d[i]);
            System.out.println("delta: " + diff + " nUlp: " + nUlp);
        }

        System.out.println("tExtrema:   " + Arrays.toString(tExtrema));
        System.out.println("tExtrema_d: " + Arrays.toString(tExtrema_d));

        for (int i = 0; i < 2; i++) {
            final BigDecimal t = tExtrema[i];
            if (t != null) {
                final BigDecimal bd = new BigDecimal(tExtrema_d[i]);
                final String bStr = toUniformString(t);
                final String dStr = toComparisonString(bd, bStr);
                System.out.println("tExtrema[" + i + "]: Exp: " + bStr);
                System.out.println("tExtrema[" + i + "]: Dbl: " + dStr);

                final double deltaT = t.subtract(bd).doubleValue();
                final double nUlpT = deltaT / Math.ulp(tExtrema_d[i]);
                System.out.println("t: " + tExtrema_d[i] + " delta(t): " + deltaT + " nUlp: " + nUlpT);
            }
        }

        System.out.println("xExtrema:   " + Arrays.toString(xExtrema));
        System.out.println("xExtrema_d: " + Arrays.toString(xExtrema_d));

        for (int i = 0; i < 2; i++) {
            final BigDecimal t = xExtrema[i];
            if (t != null) {
                final BigDecimal bd = new BigDecimal(xExtrema_d[i]);
                final String bStr = toUniformString(t);
                final String dStr = toComparisonString(bd, bStr);
                System.out.println("xExtrema[" + i + "]: Exp: " + bStr);
                System.out.println("xExtrema[" + i + "]: Dbl: " + dStr);

                final double deltaT = t.subtract(bd).doubleValue();
                final double nUlpT = deltaT / Math.ulp(xExtrema_d[i]);
                System.out.println("x: " + xExtrema_d[i] + " delta(t): " + deltaT + " nUlp: " + nUlpT);
            }
        }

        System.out.println("----");
        System.out.println("xExtrema bad:");

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

    // Find Extrema context
    private static final class FindExtremaContext {

        final double[] coeff_d = new double[4];
        final double[] deriv_coeff_d = new double[3];
        final double[] tExtrema_d = new double[2];
        final double[] xExtrema_d = new double[2];

        final BigDecimal[] coeff = new BigDecimal[4];
        final BigDecimal[] deriv_coeff = new BigDecimal[3];
        final BigDecimal[] tExtrema = new BigDecimal[2];
        final BigDecimal[] xExtrema = new BigDecimal[2];
    }

    // Find Extrema implementations (High-Accuracy) based on BigDecimal
    private final static BigDecimal TWO = new BigDecimal(2.0);
    private final static BigDecimal THREE = new BigDecimal(3.0);
    private final static BigDecimal HALF_NEG = new BigDecimal(-0.5);

    private static void findExtrema(final FindExtremaContext ctx,
                                    final double x1, final double x2, final double x3) {

        final BigDecimal[] coeff = ctx.coeff;
        final BigDecimal[] deriv_coeff = ctx.deriv_coeff;
        final BigDecimal[] tExtrema = ctx.tExtrema;
        final BigDecimal[] range = ctx.xExtrema;

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

    private static void findExtrema(final FindExtremaContext ctx,
                                    final double x1, final double x2, final double x3, final double x4) {

        final BigDecimal[] coeff = ctx.coeff;
        final BigDecimal[] deriv_coeff = ctx.deriv_coeff;
        final BigDecimal[] tExtrema = ctx.tExtrema;
        final BigDecimal[] range = ctx.xExtrema;

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

    private static void findExtremaWithDblCoeffs(final FindExtremaContext ctx,
                                                 final double x1, final double x2, final double x3, final double x4) {

        final BigDecimal[] coeff = ctx.coeff;
        final BigDecimal[] deriv_coeff = ctx.deriv_coeff;
        final BigDecimal[] tExtrema = ctx.tExtrema;
        final BigDecimal[] range = ctx.xExtrema;

        final double[] coeff_d = ctx.coeff_d;
        final double[] deriv_coeff_d = ctx.deriv_coeff_d;

        // Cubic
        final BigDecimal bx1 = new BigDecimal(x1);

        coeff[3] = new BigDecimal(coeff_d[3]);
        coeff[2] = new BigDecimal(coeff_d[2]);
        coeff[1] = new BigDecimal(coeff_d[1]);
        coeff[0] = new BigDecimal(coeff_d[0]);

        deriv_coeff[2] = new BigDecimal(deriv_coeff_d[2]);
        deriv_coeff[1] = new BigDecimal(deriv_coeff_d[1]);
        deriv_coeff[0] = new BigDecimal(deriv_coeff_d[0]);

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
        BigDecimal c = eqn[0];

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
            BigDecimal q = HALF_NEG.multiply(b.add(d)); // no div
            // ensure q has enough decimals for division:
            q = q.setScale(30, RoundingMode.HALF_EVEN);

            if (DEBUG_QUAD_SOLVER) {
                System.out.println("a (plain): " + a.toPlainString());
                System.out.println("c (plain): " + c.toPlainString());
                System.out.println("q (plain): " + q.toPlainString());
            }
            // We already tested a for being 0 above
            res[roots++] = q.divide(a, RoundingMode.HALF_EVEN);

            if (!q.equals(ZERO)) {
                // ensure c has enough decimals for division:
                c = c.setScale(30, RoundingMode.HALF_EVEN);
                res[roots++] = c.divide(q, RoundingMode.HALF_EVEN);
            }
        }
        return roots;
    }

    // Find Extrema implementations (Low-Accuracy) based on double primitive type (64bits)
    private static void findExtrema_d(final FindExtremaContext ctx,
                                      final double x1, final double x2, final double x3) {

        final double[] coeff = ctx.coeff_d;
        final double[] deriv_coeff = ctx.deriv_coeff_d;
        final double[] tExtrema = ctx.tExtrema_d;
        final double[] range = ctx.xExtrema_d;

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

    private static void findExtrema_d(final FindExtremaContext ctx,
                                      final double x1, final double x2, final double x3, final double x4) {

        final double[] coeff = ctx.coeff_d;
        final double[] deriv_coeff = ctx.deriv_coeff_d;
        final double[] tExtrema = ctx.tExtrema_d;
        final double[] range = ctx.xExtrema_d;

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

    private final static double[] disc_err = new double[2];

    // Copied from QuadCurve2D.solveQuadratic
    public static int solveQuadratic(final double[] eqn, final double[] res) {
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
            if (false && DEBUG) {
                /* LBO: Scaling equation to have (1).t^2 + (b/a).t + (c / a) = 0 */
                // bad idea => shift results by +1 ulp
                b /= a;
                c /= a;
                a = 1.0;
            }
            // From Numerical Recipes, 5.6, Quadratic and Cubic Equations
            double d = b * b - 4.0 * a * c;
            if (DEBUG_QUAD_SOLVER) {
                System.out.println("d:      " + d);
                d += Math.ulp(d);
                System.out.println("d+1ulp: " + d);

                diff_of_products(b, b, a * 4.0, c, disc_err);
                System.out.println("d_alt:     " + disc_err[0]);
                System.out.println("d_alt_err: " + disc_err[1]);
            }

            if (d < 0.0) {
                // If d < 0.0, then there are no roots
                return 0;
            }
            d = Math.sqrt(d);
            if (DEBUG_QUAD_SOLVER) {
                System.out.println("sqrt(d): " + d);
                if (false) {
                    final double fixed_d;
                    if (true) {
                        // SQRT(D + eps) = SQRT(D) + (eps / (2 SQRT(D) + 1) )
                        final double sqrt_eps = (0.5 * disc_err[1] / (1.0 + d));
                        System.out.println("sqrt_eps: " + sqrt_eps);
                        fixed_d = d + sqrt_eps;
                    } else {
                        // no precision
                        final double eps = (disc_err[1] / disc_err[0]);
                        System.out.println("eps: " + eps);

                        final double epsp1 = 1.0 + eps;
                        System.out.println("epsp1: " + epsp1);

                        final double err_sqrt = Math.sqrt(epsp1);
                        System.out.println("err_sqrt: " + err_sqrt);
                        fixed_d = d * err_sqrt;
                    }
                    System.out.println("fixed_d: " + fixed_d);
                    System.out.println("delta_d: " + (fixed_d - d));
                }
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
            if (DEBUG_QUAD_SOLVER) {
                System.out.println("q: " + q);

                System.out.println("a (dbl): " + a);
                System.out.println("c (dbl): " + c);
                System.out.println("q (dbl): " + q);
            }

            // We already tested a for being 0 above
            final double r = q / a; // first root (bigger in magnitude)
            res[roots++] = r;
            if (q != 0.0) {
                final double r2 = c / q;
                res[roots++] = r2;

                if (DEBUG_QUAD_SOLVER) {
                    System.out.println("r2: " + r2);

                    System.out.println("b2: " + (b * b));
                    System.out.println("b2-4ac: " + (b * b - 4.0 * a * c));
                    System.out.println("Same magnitude ?");

                    // use vieta formula: x1*x2 = c / a
                    final double r2_alt = (c / a) / r;
                    System.out.println("r2_alt: " + r2_alt);
                    System.out.println("delta_r: " + (r2_alt - r2));

                    final double r2_alt2 = c / (a * r);
                    System.out.println("r2_alt2: " + r2_alt2);
                    System.out.println("delta_r2: " + (r2_alt2 - r2));
                }
            }

            if (DEBUG_QUAD_SOLVER) {
                checkRoots(eqn, res, roots);
            }
        }
        return roots;
    }

    static void checkRoots(final double[] eqn, final double[] res, final int roots) {
        // Check precision:
        for (int i = 0; i < roots; i++) {
            final double t = res[i];
            final double xt = quadraticPolynom(eqn, t);
            System.out.println("t[" + i + "]: " + t + " x(t) = " + xt);

            final double eps = xt - 0.0;
            System.out.println("eps(x): " + eps);

            if (true || (eps != +0.0) && (eps != -0.0)) {
                final double dxt = quadraticDerivedPolynom(eqn, t);
                System.out.println("t[" + i + "]: " + t + " x(t) = " + xt + " dx(t) = " + dxt);

                final double u = Math.ulp(t);

                for (int j = -5; j <= 5; j++) {
                    final double te = t + u * j;
                    final double xte = quadraticPolynom(eqn, te);
                    System.out.println("t[" + j + " ulp]: " + te + " x(t) = " + xte + " (xt > 0.0): " + (xte > 0.0));
                }
            }
        }
    }

    static double quadraticPolynom(final double[] eqn, final double t) {
        return eqn[0] + t * (eqn[1] + t * (eqn[2]));
    }

    static double quadraticDerivedPolynom(final double[] eqn, final double t) {
        return eqn[1] + t * (2.0 * eqn[2]);
    }

    /*
     * diff_of_products() computes a*b-c*d with a maximum error <= 1.5 ulp
     * 
     * Claude-Pierre Jeannerod, Nicolas Louvet, and Jean-Michel Muller, 
     * "Further Analysis of Kahan's Algorithm for the Accurate Computation 
     * of 2x2 Determinants". Mathematics of Computation, Vol. 82, No. 284, 
     * Oct. 2013, pp. 2245-2264
     */
    static void diff_of_products(final double a, final double b,
                                 final double c, final double d,
                                 final double[] val_err) {
        double w = d * c;
        double e = Math.fma(-d, c, w);
        double f = Math.fma(a, b, -w);
        val_err[0] = f + e;
        // kahan sum:
        // final double t = f + e;
        // e = (t - f) - e;
        // val_err[1] = (val_err[0] - f) - e;
        val_err[1] = e; // + or - e ?

        if (DEBUG_QUAD_SOLVER) {
            System.out.println("diff_of_products: " + Arrays.toString(val_err));
        }
    }

    private FindExtremaAccuracyTest() {
        super();
    }
}
