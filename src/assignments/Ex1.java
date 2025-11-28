package assignments;

import java.util.Arrays;
import java.util.Spliterator;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
        double ans = 0;
		if (poly == null || poly.length == 0) {return ans;}
        double exponent = 1;
        for(int i = 0; i < poly.length; i++) {
            ans += poly[i] * exponent;
            exponent *= x;
        }
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);
        double f2 = f(p, x2);
        double x12 = (x1 + x2) / 2.0;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {return x12;}
        if (Math.abs(f1) < eps) {return x1;}
        if (Math.abs(f2) < eps) {return x2;}
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        }
        else {
            return root_rec(p, x12, x2, eps);
        }
    }
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
		if (lx == 2) {
            if (xx[0] == xx[1]) {
            return null;
            }
            double m = (yy[0] - yy[1]) /  (xx[0] - xx[1]);
            double b = yy[0] - m * (xx[0]);
            ans = new double[] {m, b};
        }
        else {
            double denom = (xx[0] - xx[1]) * (xx[0] - xx[2]) * (xx[1] - xx[2]);
            if (denom == 0.0) {return null;}
            double a = (xx[2] * (yy[1] - yy[0]) + xx[1] * (yy[0] - yy[2]) + xx[0] * (yy[2] - yy[1])) / denom;
            double b = (xx[2]*xx[2] * (yy[0] - yy[1]) + xx[1]*xx[1] * (yy[2] - yy[0]) + xx[0]*xx[0] * (yy[1] - yy[2])) / denom;
            double c = (xx[1] * xx[2] * (xx[1] - xx[2]) * yy[0] + xx[2] * xx[0] * (xx[2] - xx[0]) * yy[1] + xx[0] * xx[1] * (xx[0] - xx[1]) * yy[2]) / denom;
            if (a == 0) {
                ans = new double[] {c, b};
            }
            else {
                ans = new double[]{c, b, a};
            }
        }
		}
		return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        int maxLength = Math.max(p1.length,p2.length);
        for (int i = 0; i < maxLength; i++) {
        if (Math.abs(f(p1, i) - f(p2, i)) >= EPS) {
            return false;
        }
        }
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
    public static String poly(double[] poly) {
        String ans = "";
        if (poly.length == 0) return "0";
        for (int i = poly.length -1; i >= 0; i--) {
            if (poly[i] == 0) continue;
            if (i == 0) ans += poly[i];
            else if (i == 1) ans += poly[i] + "x";
            else ans += poly[i] + "x^" + i;
            int nextIndex = -1;
            for (int j = i -1; j >= 0; j--) {
                if (poly[j] != 0) {
                    nextIndex = j;
                    break;
                }
            }
            if (nextIndex != -1) {
                if (poly[nextIndex] > 0) ans += " + ";
                else ans += " - ";
            }
        }
        return ans;
    }
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
		double ans = x1;
        double[] negative;
        double[] positive;
        if (p1.length >= p2.length) {
            negative = Arrays.copyOf(p2, p2.length);
            positive = Arrays.copyOf(p1, p1.length);
        }
        else {
            negative = Arrays.copyOf(p1, p1.length);
            positive = Arrays.copyOf(p2, p2.length);
        }
        for (int i = 0; i < negative.length; i++) {
            negative[i] *= (-1);
        }
        double[] result = add(positive, negative);
        if (Math.abs(x2 - x1) < eps) {
            return (x1 + x2) / 2;
        }
        ans = root_rec(result,x1,x2,eps);
		return ans;
	}
	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        if (numberOfSegments == 0) {return ans;}
        ans = 0;
        double step = Math.abs(x2 - x1)/numberOfSegments;
        for (int i = 0; i < numberOfSegments; i++) {
            double xi = x1 + i * step;
            double xiNext = x1 + (i + 1) * step;
            double yi = f(p, xi);
            double yiNext = f(p, xiNext);
            double segmentLength = Math.sqrt(((yiNext - yi) * (yiNext - yi)) + ((xiNext - xi) * (xiNext - xi)));
            ans += segmentLength;
        }
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0;
        double h = (x2 - x1) / numberOfTrapezoid;
        for (int i = 0; i < numberOfTrapezoid; i++) {
            double xi = x1 + i * h;
            double xiNext = xi + h;
            double d1 = f(p1, xi) - f(p2, xi);
            double d2 = f(p1, xiNext) - f(p2, xiNext);
            if (d1 * d2 < 0) {
                int maxLength = Math.max(p1.length, p2.length);
                double[] negative = new double[maxLength];
                for (int j = 0; j < p2.length; j++) {
                    negative[j] = -p2[j];
                }
                double[] diff = add(p1, negative);
                double root = root_rec(diff, xi, xiNext, EPS);
                double dRoot = f(p1, root) - f(p2, root);
                ans += (Math.abs(d1) + Math.abs(dRoot)) * (root - xi) / 2.0;
                ans += (Math.abs(dRoot) + Math.abs(d2)) * (xiNext - root) / 2.0;
            } else {
                ans += (Math.abs(d1) + Math.abs(d2)) * h / 2.0;
            }
        }
        return roundTo4Decimal(ans);
    }

    public static double roundTo4Decimal(double num) {
        double round = ((int)(num * 10000)) / 10000.0;
        int decimal3 = (int)(round * 1000) % 10;
        if (decimal3 <= 2) {
            round = ((int)(round * 100)) / 100.0;
        }
        return round;
    }

    /**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
    public static double[] getPolynomFromString(String p) {
        double[] ans = ZERO;
        if (p == null) return ans;
        String[] parts = p.trim().split("\\s+");
        int maxPower = 0;
        for (String part : parts) {
            if (part.contains("x^")) {
                int power = Integer.parseInt(part.split("x\\^")[1]);
                if (power > maxPower) maxPower = power;
            }
            else if (part.contains("x")) {
                if (1 > maxPower) maxPower = 1;
            }
        }
        ans = new double[maxPower + 1];
        for (String part : parts) {
            double coeff = 0.0;
            int power = 0;
            if (part.contains("x^")) {
                String[] split = part.split("x\\^");
                String coeffStr = split[0];
                if (coeffStr.isEmpty() || coeffStr.equals("+")) coeff = 1.0;
                else if (coeffStr.equals("-")) coeff = -1.0;
                else coeff = Double.parseDouble(coeffStr);
                power = Integer.parseInt(split[1]);
            }
            else if (part.contains("x")) {
                String coeffStr = part.replace("x", "");
                if (coeffStr.isEmpty() || coeffStr.equals("+")) coeff = 1.0;
                else if (coeffStr.equals("-")) coeff = -1.0;
                else coeff = Double.parseDouble(coeffStr);
                power = 1;
            }
            else {
                if (part.isEmpty() || part.equals("+")) {
                    coeff = 1.0;
                    power = 0;
                }
                else if (part.equals("-")) {
                    coeff = -1.0;
                    power = 0;
                }
                else {
                    try {
                        coeff = Double.parseDouble(part);
                        power = 0;
                    }
                    catch (NumberFormatException e) {return ZERO;}
                }
            }
            ans[power] += coeff;
        }
        return ans;
    }
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
    public static double[] add(double[] p1, double[] p2) {
        int maxLength = Math.max(p1.length,p2.length);
        double[] result = new double[maxLength];
        for(int i=0;i < maxLength;i++) {
            double a = 0, b = 0;
            if(i< p1.length) a=p1[i];
            if(i< p2.length) b=p2[i];
            result[i] = a + b;
        }
        return result;
    }
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
    public static double[] mul(double[] p1, double[] p2) {
        double[] result = new double[p1.length + p2.length -1];
        for (int i = 0; i < p1.length; i++) {
            for (int j = 0; j < p2.length; j++) {
                result[i + j] += p1[i] * p2[j];
            }
        }
        return result;
    }

	/**
     * This function calculates the derivative of a polynomial.
     * The polynomial is stored in an array where po[i] is the
     * coefficient of x^i.
     * @param po the array of polynomial coefficients.
	 * @return ans, a new array with the coefficients of the derivative.
	 */
	public static double[] derivative (double[] po) {
        double [] ans = ZERO; // Start with a default result.
        int length = po.length; // Get the number of coefficients in the polynomial.
        // Checks if the exponent's coefficient is one or less and then returns the default value.
        if (length <= 1) return ans;
        // Otherwise the answer will be equal along the array of coefficients of the power minus one.
        else ans = new double[length -1];
        for(int i = 0; i < length -1; i++) { // A loop that runs over the size of the new array.
            ans[i]=po[i + 1] * (i + 1); // Calculating the array coefficients.
        }
        return ans; // Returning the derivative of a polynomial.
	}
}
