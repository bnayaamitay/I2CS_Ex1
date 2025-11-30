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
		if (poly == null || poly.length == 0) {return ans;} // If array is null or empty return 0.
        double exponent = 1; // Begin with x^0.
        for(int i = 0; i < poly.length; i++) { // Loop through coefficients.
            ans += poly[i] * exponent; // Add term to result.
            exponent *= x; // Update power of x.
        }
		return ans; // Return final value.
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
        double f1 = f(p, x1); // Calculate function value at x1.
        double f2 = f(p, x2); // Calculate function value at x2.
        double x12 = (x1 + x2) / 2.0; // Find midpoint between x1 and x2.
        double f12 = f(p, x12); // Calculate function value at midpoint.
        if (Math.abs(f12) < eps) {return x12;} // If midpoint value is close enough to 0 return midpoint.
        if (Math.abs(f1) < eps) {return x1;} // If value at x1 is close enough to 0 return x.
        if (Math.abs(f2) < eps) {return x2;} // If value at x2 is close enough to 0 return x2.
        // If root lies between x1 and midpoint search left side.
        if (f12 * f1 <= 0) {return root_rec(p, x1, x12, eps);}
        // Otherwise search right side.
        else {return root_rec(p, x12, x2, eps);}
    }

	/**
     * Computes a polynomial from up to 3 given points.
     * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
     * Works only for 2 or 3 points, otherwise returns null.
     *
     * @param xx x-coordinates
     * @param yy y-coordinates
     * @return polynomial coefficients
     */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null; // Default result.
        int lx = xx.length; // Length of x array.
		int ly = yy.length; // Length of y array
        // Check arrays are not null and same length and size is 2 or 3.
        if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
            if (lx == 2) {
            if (xx[0] == xx[1]) { // If x values are equal return null.
                return null;
            }
                double m = (yy[0] - yy[1]) /  (xx[0] - xx[1]); // Calculate slope.
                double b = yy[0] - m * (xx[0]); // Calculate intercept.
                ans = new double[] {m, b}; // Store line coefficients.
            }
        else {
                // Calculate denominator for quadratic.
                double denom = (xx[0] - xx[1]) * (xx[0] - xx[2]) * (xx[1] - xx[2]);
            if (denom == 0.0) {return null;}
                // Calculate quadratic coefficients.
                double a = (xx[2] * (yy[1] - yy[0]) + xx[1] * (yy[0] - yy[2]) + xx[0] * (yy[2] - yy[1])) / denom;
                double b = (xx[2]*xx[2] * (yy[0] - yy[1]) + xx[1]*xx[1] * (yy[2] - yy[0]) + xx[0]*xx[0] * (yy[1] - yy[2])) / denom;
                double c = (xx[1] * xx[2] * (xx[1] - xx[2]) * yy[0] + xx[2] * xx[0] * (xx[2] - xx[0]) * yy[1] + xx[0] * xx[1] * (xx[0] - xx[1]) * yy[2]) / denom;
                if (a == 0) { // If a is zero return linear coefficients.
                    ans = new double[] {c, b};
            }
            else {
                ans = new double[]{c, b, a};
            }
        }
		}
		return ans; // Return result.
    }

    /**
     * Checks if two polynomials are equal.
     * They are equal if f(x) values match for n+1 points,
     * where n is the maximum degree, within tolerance EPS.
     *
     * @param p1 first polynomial
     * @param p2 second polynomial
     * @return true if p1 and p2 represent the same polynomial
     */
    public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true; // Default answer is true.
        int maxLength = Math.max(p1.length,p2.length); // Use the longer polynomial length.
        for (int i = 0; i < maxLength; i++) { // Check values at each point from 0 to maxLength-1.
            // If difference is greater than EPS, return false.
            if (Math.abs(f(p1, i) - f(p2, i)) >= EPS) {return false;}
        }
		return ans; // Return true if all checks passed.
    }

	/**
     * Builds a String representation of a polynomial.
     * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
    public static String poly(double[] poly) {
        String ans = ""; // Start with empty string.
        if (poly.length == 0) return "0"; // If no coefficients, return 0.
        for (int i = poly.length -1; i >= 0; i--) { // Loop from highest power down to 0.
            if (poly[i] == 0) continue; // Skip 0 coefficients.
            if (i == 0) ans += poly[i]; // Constant term.
            else if (i == 1) ans += poly[i] + "x"; // Linear term.
            else ans += poly[i] + "x^" + i; // Higher powers.
            // Look ahead for next non-zero coefficient.
            int nextIndex = -1;
            for (int j = i -1; j >= 0; j--) {
                if (poly[j] != 0) {
                    nextIndex = j;
                    break;
                }
            }
            // Add + or - between terms.
            if (nextIndex != -1) {
                if (poly[nextIndex] > 0) ans += " + ";
                else ans += " - ";
            }
        }
        return ans; // Return final polynomial string.
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
		double ans = x1; // Default answer.
        // Arrays for positive and negative polynomial.
        double[] negative;
        double[] positive;
        if (p1.length >= p2.length) { // Choose longer polynomial as positive, shorter as negative
            negative = Arrays.copyOf(p2, p2.length);
            positive = Arrays.copyOf(p1, p1.length);
        }
        else {
            negative = Arrays.copyOf(p1, p1.length);
            positive = Arrays.copyOf(p2, p2.length);
        }
        // Negate coefficients of the "negative" polynomial
        for (int i = 0; i < negative.length; i++) {
            negative[i] *= (-1);
        }
        double[] result = add(positive, negative); // Compute difference polynomial.
        if (Math.abs(x2 - x1) < eps) { // If range is smaller than tolerance, return midpoint.
            return (x1 + x2) / 2;
        }
        ans = root_rec(result,x1,x2,eps); // Otherwise, find root of difference polynomial.
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
		double ans = x1; // Default.
        if (numberOfSegments == 0) {return ans;} // If no segments, return x1.
        ans = 0; // Reset answer to 0.
        double step = Math.abs(x2 - x1)/numberOfSegments; // Step size along x.
        for (int i = 0; i < numberOfSegments; i++) { // Loop through each segment.
            // Current x and next x.
            double xi = x1 + i * step;
            double xiNext = x1 + (i + 1) * step;
            // Polynomial values at xi and xiNext.
            double yi = f(p, xi);
            double yiNext = f(p, xiNext);
            // Segment length using distance formula.
            double segmentLength = Math.sqrt(((yiNext - yi) * (yiNext - yi)) + ((xiNext - xi) * (xiNext - xi)));
            ans += segmentLength; // Add segment length to total.
        }
		return ans; // Return total arc length.
    }

    /**
     * Approximates the area between two polynomials (p1, p2)
     * over the range [x1, x2] using the trapezoid rule.
     *
     * @param p1 first polynomial
     * @param p2 second polynomial
     * @param x1 start of range
     * @param x2 end of range
     * @param numberOfTrapezoid number of trapezoids (samples)
     * @return approximated area between p1 and p2
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0; // Initialize result.
        double h = (x2 - x1) / numberOfTrapezoid; // Step size.
        for (int i = 0; i < numberOfTrapezoid; i++) { // Loop over trapezoids.
            // Current x and next x
            double xi = x1 + i * h;
            double xiNext = xi + h;
            // Difference between polynomials at xi and xiNext
            double d1 = f(p1, xi) - f(p2, xi);
            double d2 = f(p1, xiNext) - f(p2, xiNext);
            if (d1 * d2 < 0) { // If the functions cross between xi and xiNext.
                int maxLength = Math.max(p1.length, p2.length); // Build array for -p2.
                double[] negative = new double[maxLength];
                for (int j = 0; j < p2.length; j++) {
                    negative[j] = -p2[j];
                }
                double[] diff = add(p1, negative); // Compute difference polynomial p1 - p2.
                double root = root_rec(diff, xi, xiNext, EPS); // Find root (intersection point) between xi and xiNext.
                double dRoot = f(p1, root) - f(p2, root); // Difference at root.
                ans += (Math.abs(d1) + Math.abs(dRoot)) * (root - xi) / 2.0;  // Add area from xi to root.
                ans += (Math.abs(dRoot) + Math.abs(d2)) * (xiNext - root) / 2.0; // Add area from root to xiNext.
            } else { // Normal trapezoid area if no crossing.
                ans += (Math.abs(d1) + Math.abs(d2)) * h / 2.0;
            }
        }
        return roundTo4Decimal(ans); // Round result to 4 decimals.
    }

    /**
     * Rounds a number to 4 decimal places.
     * If the 3rd decimal digit is 4 or less, it further rounds to 2 decimal places.
     *
     * @param num the input number
     * @return the rounded number
     */
    public static double roundTo4Decimal(double num) {
        double round = ((int)(num * 10000)) / 10000.0; // Truncate to 4 decimals.
        int decimal3 = (int)(round * 1000) % 10; // Get 3rd decimal digit.
        if (decimal3 <= 4) { // If digit is less than or equal to 4, round to 2 decimals.
                round = ((int)(round * 100)) / 100.0;
        }
        return round; // Return result
    }

    /**
     * Converts a polynomial string into an array of coefficients.
     * Each array index represents the power of x, and the value is its coefficient.
     *
     * @param p a String representing a polynomial function
     * @return a double array of coefficients for the polynomial
     */
    public static double[] getPolynomFromString(String p) {
        double[] ans = ZERO;  // Start with a default result.
        if (p == null) return ans; // If the input string is null, return ZERO.
        String[] parts = p.trim().split("\\s+"); // Split the string into parts by whitespace.
        // Track the highest power of x found
        int maxPower = 0;
        for (String part : parts) { // Loop1- find the maximum power in the polynomial.
            if (part.contains("x^")) {  // Extract the power after "x^".
                int power = Integer.parseInt(part.split("x\\^")[1]);
                if (power > maxPower) maxPower = power;
            }
            else if (part.contains("x")) {
                if (1 > maxPower) maxPower = 1; // If it only contains "x", the power is 1.
            }
        }
        ans = new double[maxPower + 1]; // Create the result array with size = maxPower + 1.
        for (String part : parts) { // Loop2- parse each part and fill coefficients.
            double coeff = 0.0; // coefficient value.
            int power = 0; // power of x.
            if (part.contains("x^")) { // If the term like "3x^2" or "-x^4".
                String[] split = part.split("x\\^");
                String coeffStr = split[0];
                // Handle missing or sign-only coefficients.
                if (coeffStr.isEmpty() || coeffStr.equals("+")) coeff = 1.0;
                else if (coeffStr.equals("-")) coeff = -1.0;
                else coeff = Double.parseDouble(coeffStr);
                power = Integer.parseInt(split[1]); // Parse the power after "^".
            }
            else if (part.contains("x")) { // If the term is "x" or "-x".
                String coeffStr = part.replace("x", "");
                if (coeffStr.isEmpty() || coeffStr.equals("+")) coeff = 1.0;
                else if (coeffStr.equals("-")) coeff = -1.0;
                else coeff = Double.parseDouble(coeffStr);
                power = 1;
            }
            else { // If the term is a constant (no x).
                if (part.isEmpty() || part.equals("+")) {
                    coeff = 1.0;
                    power = 0;
                }
                else if (part.equals("-")) {
                    coeff = -1.0;
                    power = 0;
                }
                else {
                    try { // Try parsing a number, otherwise return ZERO.
                        coeff = Double.parseDouble(part);
                        power = 0;
                    }
                    catch (NumberFormatException e) { return ZERO; }
                }
            }
            ans[power] += coeff; // Add the coefficient to the correct power index.
        }
        return ans; // Return the final array of coefficients.
    }

	/**
     * Adds two polynomials p1 and p2.
     * Each polynomial is stored in an array where p[i] is the coefficient of x^i.
     * The function returns a new array with the coefficients of the sum polynomial.
     *
     * @param p1 coefficients of the first polynomial
     * @param p2 coefficients of the second polynomial
     * @return coefficients of the resulting polynomial after addition
     */
    public static double[] add(double[] p1, double[] p2) {
        int maxLength = Math.max(p1.length,p2.length); // Find the maximum length.
        double[] result = new double[maxLength]; // Create a result array with that maximum length.
        for(int i=0;i < maxLength;i++) { // Loop through all positions up to the longest polynomial.
            double a = 0, b = 0; // Default coefficients are 0 if index is out of bounds.
            if(i< p1.length) a=p1[i]; // Take coefficient from p1 if it exists.
            if(i< p2.length) b=p2[i]; // Take coefficient from p2 if it exists.
            result[i] = a + b; // Add the two coefficients together.
        }
        return cleanFromZeros(result); // Return the coefficients of the sum polynomial.
    }

	/**
	 * Multiplies two polynomials (p1 and p2).
     * Each polynomial is stored in an array where p[i] is the coefficient of x^i.
     * The function returns a new array with the coefficients of the product polynomial.
     *
     * @param p1 coefficients of the first polynomial
     * @param p2 coefficients of the second polynomial
     * @return coefficients of the resulting polynomial after multiplication
	 */
    public static double[] mul(double[] p1, double[] p2) {
        double[] result = new double[p1.length + p2.length - 1]; // The array size is p1.length + p2.length - 1.
        for (int i = 0; i < p1.length; i++) { // Loop through all coefficients of the first polynomial.
            for (int j = 0; j < p2.length; j++) { // Loop through all coefficients of the second polynomial.
                result[i + j] += p1[i] * p2[j]; // Multiply coefficients and add to the correct position (i+j)
            }
        }
        return cleanFromZeros(result); // Return the coefficients of the product polynomial.
    }

    /**
     * Removes unnecessary zeros from the end of a polynomial array.
     * Returns the shortened array, or {0} if all elements are zero.
     */
    public static double[] cleanFromZeros(double[] p) {
        int i = p.length - 1;
        while (i >= 0 && p[i] == 0) { // A loop that runs backwards until a non-zero term is found.
            i--;
        }
        if (i < 0) {return ZERO;} // If all terms are zero
        return Arrays.copyOf(p, i + 1); // Otherwise, return the array up to the last non-zero element.
    }

	/**
     * This function calculates the derivative of a polynomial.
     * The polynomial is stored in an array where po[i] is the
     * coefficient of x^i.
     *
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
