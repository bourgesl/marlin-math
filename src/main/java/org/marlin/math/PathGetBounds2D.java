/*******************************************************************************
 * JMMC project ( http://www.jmmc.fr ) - Copyright (C) CNRS.
 ******************************************************************************/
package org.marlin.math;

import java.awt.Shape;
import java.awt.geom.PathIterator;
import java.awt.geom.QuadCurve2D;
import java.awt.geom.Rectangle2D;

/**
 *
 * @author @mickleness
 */
public class PathGetBounds2D {

    private PathGetBounds2D() {
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
        final double[] coeff = new double[4];

        // define the derivative's coefficients
        final double[] deriv_coeff = new double[3];

        final double[] coords = new double[6];
        final double[] tExtrema = new double[2];
        boolean started = false;
        double leftX = 0.0;
        double rightX = 0.0;
        double topY = 0.0;
        double bottomY = 0.0;
        double lastX = 0.0;
        double lastY = 0.0;

        // Whenever we have to examine cubic or quadratic extrema that change
        // our bounding box: we run the risk of machine error that may produce
        // a box that is slightly too small. But the contract of this method
        // says we should err on the side of being too large.
        // So to address this: we take using the upper limit of numerical error
        // caused by the polynomial evaluation (horner scheme).

        for (; !pi.isDone(); pi.next()) {
            final int type = pi.currentSegment(coords);
            switch (type) {
                case PathIterator.SEG_MOVETO:
                    if (!started) {
                        started = true;
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
                        coeff[2] = (coords[2] - coords[0]) - dx21;  // A = P3 - P0 - 2 P2
                        coeff[1] = 2.0 * dx21;                      // B = 2 (P2 - P1)
                        coeff[0] = lastX;                           // C = P1

                        deriv_coeff[0] = coeff[1];
                        deriv_coeff[1] = 2.0 * coeff[2];

                        double t = -deriv_coeff[0] / deriv_coeff[1];
                        if (t > 0.0 && t < 1.0) {
                            double x = coeff[0] + t * (coeff[1] + t * coeff[2]);

                            // error condition = sum ( abs (coeff) ):
                            final double margin = Math.ulp( Math.abs(coeff[0])
                                    + Math.abs(coeff[1]) + Math.abs(coeff[2]));

                            if (x - margin < leftX) {
                                leftX = x - margin;
                            }
                            if (x + margin > rightX) {
                                rightX = x + margin;
                            }
                        }
                    }
                    if (coords[1] < topY || coords[1] > bottomY) {
                        final double dy21 = (coords[1] - lastY);
                        coeff[2] = (coords[3] - coords[1]) - dy21;
                        coeff[1] = 2.0 * dy21;
                        coeff[0] = lastY;

                        deriv_coeff[0] = coeff[1];
                        deriv_coeff[1] = 2.0 * coeff[2];

                        double t = -deriv_coeff[0] / deriv_coeff[1];
                        if (t > 0.0 && t < 1.0) {
                            double y = coeff[0] + t * (coeff[1] + t * coeff[2]);
                            
                            // error condition = sum ( abs (coeff) ):
                            final double margin = Math.ulp( Math.abs(coeff[0])
                                    + Math.abs(coeff[1]) + Math.abs(coeff[2]));
                            
                            if (y - margin < topY) {
                                topY = y - margin;
                            }
                            if (y + margin > bottomY) {
                                bottomY = y + margin;
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
                        coeff[3] = (coords[4] - lastX) - dx32;  // A = P3 - P0 - 3 (P2 - P1) = (P3 - P0) + 3 (P1 - P2)
                        coeff[2] = (dx32 - dx21);               // B = 3 (P2 - P1) - 3(P1 - P0) = 3 (P2 + P0) - 6 P1
                        coeff[1] = dx21;                        // C = 3 (P1 - P0)
                        coeff[0] = lastX;                       // D = P0

                        deriv_coeff[0] = coeff[1];
                        deriv_coeff[1] = 2.0 * coeff[2];
                        deriv_coeff[2] = 3.0 * coeff[3];

                        // solveQuadratic should be improved to get correct t extrema (1 ulp):
                        final int tExtremaCount = QuadCurve2D.solveQuadratic(deriv_coeff, tExtrema);
                        if (tExtremaCount > 0) {
                            // error condition = sum ( abs (coeff) ):
                            final double margin = Math.ulp(Math.abs(coeff[0])
                                    + Math.abs(coeff[1]) + Math.abs(coeff[2])
                                    + Math.abs(coeff[3]));

                            for (int i = 0; i < tExtremaCount; i++) {
                                final double t = tExtrema[i];
                                if (t > 0.0 && t < 1.0) {
                                    double x = coeff[0] + t * (coeff[1] + t * (coeff[2] + t * coeff[3]));
                                    if (x - margin < leftX) {
                                        leftX = x - margin;
                                    }
                                    if (x + margin > rightX) {
                                        rightX = x + margin;
                                    }
                                }
                            }
                        }
                    }
                    if (coords[1] < topY || coords[1] > bottomY || coords[3] < topY || coords[3] > bottomY) {
                        final double dy32 = 3.0 * (coords[3] - coords[1]);
                        final double dy21 = 3.0 * (coords[1] - lastY);
                        coeff[3] = (coords[5] - lastY) - dy32;
                        coeff[2] = (dy32 - dy21);
                        coeff[1] = dy21;
                        coeff[0] = lastY;

                        deriv_coeff[0] = coeff[1];
                        deriv_coeff[1] = 2.0 * coeff[2];
                        deriv_coeff[2] = 3.0 * coeff[3];

                        int tExtremaCount = QuadCurve2D.solveQuadratic(deriv_coeff, tExtrema);
                        if (tExtremaCount > 0) {
                            // error condition = sum ( abs (coeff) ):
                            final double margin = Math.ulp(Math.abs(coeff[0])
                                    + Math.abs(coeff[1]) + Math.abs(coeff[2])
                                    + Math.abs(coeff[3]));

                            for (int i = 0; i < tExtremaCount; i++) {
                                double t = tExtrema[i];
                                if (t > 0.0 && t < 1.0) {
                                    double y = coeff[0] + t * (coeff[1] + t * (coeff[2] + t * coeff[3]));
                                    if (y - margin < topY) {
                                        topY = y - margin;
                                    }
                                    if (y + margin > bottomY) {
                                        bottomY = y + margin;
                                    }
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
        if (started) {
            return new Rectangle2D.Double(leftX, topY, rightX - leftX, bottomY - topY);
        }

        // there's room to debate what should happen here, but historically we return a zeroed
        // out rectangle here. So for backwards compatibility let's keep doing that:
        return new Rectangle2D.Double();
    }

}
