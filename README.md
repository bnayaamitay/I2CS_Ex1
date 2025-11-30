# Ex1
This is my first project I have undertaken as part of my academic studies at the university.
The project focuses on building and testing basic functions in Java, using arrays to represent polynomial coefficients.

## Features

- `root_rec(double[] p, double x1, double x2, double eps)`  
  Recursively finds a root of the polynomial within the range `[x1, x2]`.  
  Returns an `x` value such that `|p(x)| < eps`, assuming `p(x1) * p(x2) <= 0`.

- `f(double[] poly, double x)`  
  Computes the value of the polynomial function at a given `x`.  
  Returns `f(x)`.

- `add(double[] p1, double[] p2)`  
  Adds two polynomials.  
  Each polynomial is stored in an array where `p[i]` is the coefficient of `x^i`.  
  Returns a new array with the coefficients of the sum polynomial.

- `area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid)`  
  Approximates the area between two polynomials over the range `[x1, x2]` using the trapezoid rule.  
  Returns the estimated area based on the given number of trapezoids (samples).

- `derivative(double[] po)`  
  Calculates the derivative of a polynomial.  
  The polynomial is stored in an array where `po[i]` is the coefficient of `x^i`.  
  Returns a new array with the coefficients of the derivative polynomial.

- `roundTo4Decimal(double num)`  
  Rounds a number to 4 decimal places.  
  If the 3rd decimal digit is 4 or less, it further rounds to 2 decimal places.

- `poly(double[] p)`  
  Converts a polynomial array into a human-readable string.  
  Each array index represents the power of `x`, and the value is its coefficient.  
  For example, `{2, 0, 3}` becomes `"3.0x^2 + 2.0"`.

- `getPolynomFromString(String p)`  
  Converts a polynomial string into an array of coefficients.  
  Each array index represents the power of `x`, and the value is its coefficient.  
  For example, `"3.0x^2 + 2.0"` becomes `{2, 0, 3}`.

- `mul(double[] p1, double[] p2)`  
  Multiplies two polynomials.  
  Each polynomial is stored in an array where `p[i]` is the coefficient of `x^i`.  
  Returns a new array with the coefficients of the product polynomial.

- `cleanFromZeros(double[] p)`  
  Removes unnecessary zeros from the end of a polynomial array.  
  Returns the shortened array, or `{0}` if all elements are zero.

- `PolynomFromPoints(double[] xx, double[] yy)`  
  Computes a polynomial from 2 or 3 given points.  
  Returns the polynomial coefficients, or `null` if the number of points is not valid.

- `equals(double[] p1, double[] p2)`  
  Checks if two polynomials are equal.  
  They are considered equal if their values match for `n+1` test points (where `n` is the maximum degree), within a tolerance `EPS`.

- `sameValue(double[] p1, double[] p2, double x1, double x2, double eps)`
  Finds an `x` value in the range `[x1, x2]` where the two polynomials have approximately the same value.  
  Returns an `x` such that `|p

## Image result
<img width="1002" height="1010" alt="image" src="https://github.com/user-attachments/assets/37b020b4-75d2-4cd3-80ce-4b116a08469e" />

