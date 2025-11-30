# Ex1
This is my first project I have undertaken as part of my academic studies at the university.
The project focuses on building and testing basic functions in Java, using arrays to represent polynomial coefficients.

## Features

- `add(double[] p1, double[] p2)`  
  Adds two polynomials.  
  Each polynomial is stored in an array where `p[i]` is the coefficient of `x^i`.  
  Returns a new array with the coefficients of the sum polynomial.

- `sub(double[] p1, double[] p2)`  
  Subtracts one polynomial from another.  
  Returns a new array with the coefficients of the difference polynomial.

- `mul(double[] p1, double[] p2)`  
  Multiplies two polynomials.  
  Each polynomial is stored in an array where `p[i]` is the coefficient of `x^i`.  
  Returns a new array with the coefficients of the product polynomial.

- `cleanFromZeros(double[] p)`  
  Removes unnecessary zeros from the end of a polynomial array.  
  Returns the shortened array, or `{0}` if all elements are zero.

- `f(double[] poly, double x)`  
  Computes the value of the polynomial function at a given `x`.  
  Returns `f(x)`.

- `root_rec(double[] p, double x1, double x2, double eps)`  
  Recursively finds a root of the polynomial within the range `[x1, x2]`.  
  Returns an `x` value such that `|p(x)| < eps`, assuming `p(x1) * p(x2) <= 0`.

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

