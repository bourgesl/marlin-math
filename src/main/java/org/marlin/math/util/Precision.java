package org.marlin.math.util;

/**
 * Basic ulp scaling test
 */
public class Precision {

    /* epsilon value to roughly evaluate Math.ulp((double) x) ~ 1e-15 */
    private final static double EPS = 2.5e-16d;

    public static void main(String[] unused) {
        final double dbl26bits = Math.pow(2.0, 26);
        System.out.println("dbl26bits: " + dbl26bits);
        
        final double dbl53bits = Math.pow(2.0, 53);
        System.out.println("dbl53bits: " + dbl53bits);

        for (double d = 1e-6; d < 1e20; d = 10.0 * d) {
            System.out.println("---");
            System.out.println("ulp(" + d + "d) = " + Math.ulp(d));
            System.out.println("rough ulp(" + (d) + ") = " + Math.abs(EPS * d));
            System.out.println("ratio rough = " + (Math.abs(EPS * d) / Math.ulp(d)));
        }
    }

    private Precision() {
        super();
    }
}
