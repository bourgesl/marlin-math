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
package org.marlin.math.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.Arrays;

/**
 * Simple class to provide utility statistical functions
 */
@SuppressWarnings("UseOfSystemOutOrSystemErr")
public final class StatUtils {

    private StatUtils() {
        // forbidden constructor
    }

    public static double min(final double[] array) {
        double min = Double.POSITIVE_INFINITY;

        for (int n = 0; n < array.length; n++) {
            if (array[n] < min) {
                min = array[n];
            }
        }
        return min;
    }

    public static double max(final double[] array) {
        double max = Double.NEGATIVE_INFINITY;

        for (int n = 0; n < array.length; n++) {
            if (array[n] > max) {
                max = array[n];
            }
        }
        return max;
    }

    public static double mean(final double[] array) {
        double sample, sum = 0.0;
        int ns = 0;
        for (int n = 0; n < array.length; n++) {
            sample = array[n];
            // No Compensated-summation (double):
            if (!Double.isNaN(sample)) {
                sum += sample;
                ns++;
            }
        }
        return (ns != 0) ? (sum / ns) : 0.0;
    }

    public static double[] moments(final double[] array) {
        final double[] moments = new double[4];
        moments(array, moments);
        return moments;
    }

    public static void moments(final double[] array, final double[] moments) {
        final double mean = mean(array);

        double sample, diff;
        double sum_diff = 0.0;
        double sum_diff2 = 0.0;

        for (int n = 0; n < array.length; n++) {
            sample = array[n];
            // Compensated-summation variant for better numeric precision:
            diff = sample - mean;
            sum_diff += diff;
            sum_diff2 += diff * diff;
        }

        // variance(norm):
        // note: this algorithm ensures correctness (stable) even if the mean used in diff is wrong !
        final double variance = (sum_diff2 - (sum_diff * sum_diff) / array.length) / (array.length - 1);

        final double stddev = Math.sqrt(variance);

        // Moments ordre 3 et 4:
        double sum_diff3 = 0.0;
        double sum_diff4 = 0.0;

        for (int n = 0; n < array.length; n++) {
            sample = array[n];
            // Compensated-summation variant for better numeric precision:
            diff = (sample - mean) / stddev;
            sum_diff3 += diff * diff * diff;
            sum_diff4 += diff * diff * diff * diff;
        }

        final double asymetry = sum_diff3 / array.length;
        final double kurtosis = (sum_diff4 / array.length) - 3.0; // normalised

        // output:
        moments[0] = mean;
        moments[1] = variance;
        moments[2] = asymetry;
        moments[3] = kurtosis;
    }

    public static double naiveSum(double[] values) {
        final double[] state = new double[1]; // sum
        state[0] = 0.0;
        for (int i = 0; i < values.length; i++) {
            state[0] += values[i];
        }
        return state[0];
    }

    public static double kahanSum(double[] values) {
        final double[] state = new double[2]; // sum | error
        state[0] = 0.0;
        state[1] = 0.0;
        for (int i = 0; i < values.length; i++) {
            final double y = values[i] - state[1];
            final double t = state[0] + y;
            state[1] = (t - state[0]) - y;
            state[0] = t;
        }
        return state[0];
    }

    // --- TEST ---
    public static void main(String[] args) throws IOException {
        final boolean TEST_SUM = false;
        final boolean TEST_DIST = true;
        final boolean DO_DUMP = false;

        // Test kahan sum:
        if (TEST_SUM) {
            final double[] values = new double[10 * 1024 * 1024];
            testSum(values, 1.0e-8);
            testSum(values, 1.0);
            testSum(values, 1.0e8);
        }
    }

    // sum tests
    private static void testSum(final double[] values, final double val) {
        Arrays.fill(values, val);
        values[0] = 1.0;

        final double naive = naiveSum(values);
        System.out.println("naiveSum[1 + " + val + " x " + values.length + "]: " + naive);
        final double kahan = kahanSum(values);
        System.out.println("kahanSum[1 + " + val + " x " + values.length + "]: " + kahan);
        System.out.println("delta: " + (naive - kahan));
    }

    // --- utility functions ---
    public static void saveArray(final File file, final double[] array) throws IOException {
        if (file != null) {
            final int len = array.length;
            final int capacity = (len + 1) * 40;

            final StringBuilder sb = new StringBuilder(capacity);
            sb.append("# X\n");

            for (int i = 0; i < len; i++) {
                sb.append(array[i]).append('\n');
            }

            System.out.println("Writing file: " + file.getAbsolutePath());
            writeFile(file, sb.toString());
        }
    }

    private static void writeFile(final File file, final String content) throws IOException {
        final Writer w = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file), "UTF-8"));
        try {
            w.write(content);
        } finally {
            try {
                w.close();
            } catch (IOException ioe) {
                System.err.println("Failed writing file: " + file.getAbsolutePath());
                ioe.printStackTrace(System.err);
            }
        }
    }
}
