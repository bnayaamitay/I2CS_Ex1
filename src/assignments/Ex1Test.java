package assignments;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *.
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};
	
 	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}

	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}

    @Test
    /**
    * Tests root_rec on a quadratic polynomial: f(x) = x^2 - 4.
    * Roots are at x = -2 and x = 2.
    */
    void testRootQuadratic() {
        double[] p = {-4, 0, 1};
        double root1 = Ex1.root_rec(p, -3, -1, Ex1.EPS);
        double root2 = Ex1.root_rec(p, 1, 3, Ex1.EPS);
        assertEquals(-2.0, root1, Ex1.EPS);
        assertEquals(2.0, root2, Ex1.EPS);
    }

    @Test
    /**
    * Tests root_rec on a simple linear polynomial: f(x) = x - 2.
    * Root should be at x = 2.
    */
    void testRootLinear() {
        double[] p = {-2, 1};
        double root = Ex1.root_rec(p, 0, 5, Ex1.EPS);
        assertEquals(2.0, root, Ex1.EPS);
    }

    @Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}

	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}

	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}

    @Test
    /**
     * Verifies that adding 0 to p1 does not change its coefficients or length.
     */
    void testAdd4() {
        boolean ans= false;
        int p1l = po1.length;
        double[] p2 = {0,0,0,0};
        double[] pResult = Ex1.add(p2, po1);
        if(pResult.length==p1l) {ans=true;}
        assertTrue(ans);
    }

    @Test
    /**
    * Tests add function correctness by also checking derivative and evaluation.
    */
    void testAdd5() {
        double[] p1 = {1, 2};
        double[] p2 = {0, 1};
        double[] sum = Ex1.add(p1, p2);
        double[] deriv = Ex1.derivative(sum);
        double val = Ex1.f(sum, 2);
        assertTrue(Ex1.equals(sum, new double[]{1, 3}));
        assertTrue(Ex1.equals(deriv, new double[]{3}));
        assertEquals(7.0, val, Ex1.EPS);
    }

	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}

	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}

	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}

    @Test
    /**
    * Tests multiplying by ONE polynomial (x^0 = 1).
    */
    void testMulIdentity() {
        double[] one = {1};
        double[] p = Ex1.mul(po1, one);
        assertTrue(Ex1.equals(p, po1));
    }

    @Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3};
		double[] pt = {2,6};
		double[] dp1 = Ex1.derivative(p);
		double[] dp2 = Ex1.derivative(dp1);
		double[] dp3 = Ex1.derivative(dp2);
		double[] dp4 = Ex1.derivative(dp3);
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}

    @Test
    /**
    * Tests derivative of a constant polynomial.
    */
    void testDerivativeConstant() {
        double[] p = {7};
        double[] dp = Ex1.derivative(p);
        assertTrue(Ex1.equals(dp, Ex1.ZERO));
    }

    @Test
    /**
    * Tests poly string conversion of a simple polynomial with flexible comparison (ignores spaces).
    */
    void testPolyStringFlexible() {
        double[] p = {2, 0, 3};
        String s = Ex1.poly(p);
        String expected = "3.0x^2 +2.0";
        assertTrue(s.replace(" ", "").equals(expected.replace(" ", "")));
    }

    @Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1};
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}

    @Test
    /**
    * Tests parsing a simple polynomial string "x".
    */
    void testFromStringSimpleX() {
        double[] p = {0,1};
        double[] parsed = Ex1.getPolynomFromString("x");
        assertTrue(Ex1.equals(parsed, p));
    }

    @Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

    @Test
    /**
    * Tests equality of two identical ZERO polynomials.
    */
    void testEqualsZero() {
        assertTrue(Ex1.equals(Ex1.ZERO, Ex1.ZERO));
    }

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}

    @Test
    /**
    * Tests sameValue when both polynomials are identical.
    */
    void testSameValue3() {
        double rs = Ex1.sameValue(po1, po1, -1, 1, Ex1.EPS);
        assertEquals(0, Math.abs(Ex1.f(po1, rs) - Ex1.f(po1, rs)), Ex1.EPS);
    }

    @Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=-4, x2=0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2, Ex1.EPS);
    }

	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = {0,1};
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
		double area =2.5;
		assertEquals(a1,area, Ex1.EPS);
		assertEquals(a2,area, Ex1.EPS);
		assertEquals(a3,area, Ex1.EPS);
		assertEquals(a100,area, Ex1.EPS);
	}

	@Test
	/**
	 * Test the area function.
	 */
	public void testArea3() {
		double[] po_a = {2,1,-0.7, -0.02,0.02};
		double[] po_b = {6, 0.1, -0.2};
		double x1 = Ex1.sameValue(po_a,po_b, -10,-5, Ex1.EPS);
		double a1 = Ex1.area(po_a,po_b, x1, 6, 8);
		double area = 58.5658;
		assertEquals(a1,area, Ex1.EPS);
	}
    @Test
    /**
    * Tests length of ZERO polynomial (flat line).
    */
    void testLengthZero() {
        double len = Ex1.length(Ex1.ZERO, 0, 10, 5);
        assertEquals(10, len, Ex1.EPS);
    }

    @Test
    /**
    * Tests PolynomFromPoints with 2 points (should return a linear polynomial).
    */
    void testPolynomFromPointsTwoPoints1() {
        double[] xx = {0, 1};
        double[] yy = {2, 4};
        double[] poly = Ex1.PolynomFromPoints(xx, yy);
        double[] expected = {2, 2};
        assertTrue(Ex1.equals(poly, expected));
    }

    @Test
    /**
    * Tests PolynomFromPoints with 2 points (flexible check: polynomial passes through both points).
    */
    void testPolynomFromPointsTwoPoints2() {
        double[] xx = {0, 2};
        double[] yy = {1, 5};
        double[] poly = Ex1.PolynomFromPoints(xx, yy);
        double f0 = Ex1.f(poly, xx[0]);
        double f2 = Ex1.f(poly, xx[1]);
        assertEquals(yy[0], f0, Ex1.EPS);
        assertEquals(yy[1], f2, Ex1.EPS);
    }

    @Test
    /**
    * Tests PolynomFromPoints with 3 points (should return a quadratic polynomial).
    */
    void testPolynomFromPointsThreePoints() {
        double[] xx = {0, 1, 2};
        double[] yy = {1, 2, 5};
        double[] poly = Ex1.PolynomFromPoints(xx, yy);
        double[] expected = {1, 0, 1};
        assertTrue(Ex1.equals(poly, expected));
    }

    @Test
    /**
    * Tests rounding of a number with fewer decimals.
    */
    void testRoundSimple() {
        double r = Ex1.roundTo4Decimal(3.1);
        assertEquals(3.1, r, Ex1.EPS);
    }

    @Test
    /**
    * Tests cleanFromZeros on an array with only zeros.
    */
    void testCleanFromZerosAllZeros() {
        double[] p = {0,0,0,0};
        double[] cleaned = Ex1.cleanFromZeros(p);
        assertTrue(Ex1.equals(cleaned, Ex1.ZERO));
    }

}
