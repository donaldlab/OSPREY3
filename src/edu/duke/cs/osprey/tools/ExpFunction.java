/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

///////////////////////////////////////////////////////////////////////////////////////////////
//ExpFunction.java
//
//Version:           2.1 beta
//
//
//  authors:
//	  initials    name                 organization                email
// ---------   -----------------    ------------------------    ----------------------------
//  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////
/**
 * Written by Ivelin Georgiev (2004-2009)
 * 
 */
package edu.duke.cs.osprey.tools;


import java.io.Serializable;
import java.math.*;

/**
 * Manages the computation of exp(x) for large values of x, using BigDecimal;
 * 		For large values of x, the standard Math.exp(x) cannot be used, since it only has double precision;
 * Implements pow() for integer powers of a BigDecimal number and an approximation to the natural logarithm of a BigDecimal number
 * 
 */
@SuppressWarnings("serial")
public class ExpFunction implements Serializable {

	public static final BigDecimal exp = new BigDecimal("2.71828182845904523536"); //Euler's number to 20 decimal digits
	public static MathContext mc = new MathContext(100, RoundingMode.HALF_EVEN);

	public final int maxPrecision; //the number of decimal digits to which the BigDecimal numbers must be accurate

	public ExpFunction() {
		this.mathContext = null;
		this.maxPrecision = 8;
	}

	private final MathContext mathContext;

	public ExpFunction(MathContext mathContext) {
		this.mathContext = mathContext;
		this.maxPrecision = mathContext.getPrecision();
	}

	//Computes exp(x) using BigDecimal arithmetic for large x or the standard Math.exp() function for small x;
	//		If x is large, it is divided into its integer and fractional part; BigDecimal is used to compute
	//			exp to the power of the integer-part of x; the standard Math.exp() function is used to compute
	//			exp(fractional power), since the fractional part < 1; the two results are then multiplied to
	//			obtain the BigDecimal exp(x)
	public BigDecimal exp(double x){

		final double minX = 0.0; //only use BigDecimal arithmetic for (x>=minX)

		BigDecimal expX = null;

		if (x<minX){ //x is small enough, so use the standard Math.exp() function
			expX = new BigDecimal(Math.exp(x));
		}
		else { //x is large, so use the BigDecimal arithmetic approach
			int integerPart = (int)Math.floor(x);
			double fractionPart = x - integerPart;

			BigDecimal intPartExp = pow(exp,integerPart);
			BigDecimal fractPartExp = new BigDecimal(Math.exp(fractionPart));

			expX = intPartExp.multiply(fractPartExp);
		}

		expX = expX.setScale(maxPrecision,4); //rounding is ROUND_HALF_UP (standard rounding: up for next digit >=5, down otherwise)

		return expX;
	}

	// this is probably enough precision for us, but there's more if we need it:
	// http://www-history.mcs.st-and.ac.uk/HistTopics/e_10000.html
	private BigDecimal Eprecise = new BigDecimal(
		"2.718281828459045235360287471352662497757247093699959574966967627724076630353"
		+ "547594571382178525166427427466391932003059921817413596629043572900334295260"
		+ "595630738132328627943490763233829880753195251019011573834187930702154089149"
		+ "934884167509244761460668082264800168477411853742345442437107539077744992069"
		+ "551702761838606261331384583000752044933826560297606737113200709328709127443"
		+ "747047230696977209310141692836819025515108657463772111252389784425056953696"
		+ "770785449969967946864454905987931636889230098793127736178215424999229576351"
		+ "482208269895193668033182528869398496465105820939239829488793320362509443117"
	);

	/** no really, implement exp(x) without loss of precision */
	public BigDecimal expPrecise(double x) {

		if (x == 0.0) {
			// exp(0) = 1, easy peasy
			return BigDecimal.ONE;
		}

		int intPart;
		if (x > 0) {
			intPart = (int)Math.floor(x);
		} else {
			intPart = (int)Math.ceil(x);
		}
		double fracPart = x - intPart;

		BigDecimal intPartExp = Eprecise.pow(intPart, mathContext);
		BigDecimal fractPartExp = new BigDecimal(Math.exp(fracPart));

		return intPartExp.multiply(fractPartExp, mathContext);
	}

	//Returns the BigDecimal number num to the power a	
	BigDecimal pow(BigDecimal num, int a) {
		return num.pow(a);
	}

	//Returns an approximation to the natural logarithm of the BigDecimal number num
	public BigDecimal log(BigDecimal num){

		if (num.compareTo(new BigDecimal("0.0"))<0){ //num is negative
			throw new IllegalArgumentException("log of a negative number: " + num);
		}

		BigDecimal sum = new BigDecimal("0.0");
		BigDecimal x = num;

		if (num.compareTo(new BigDecimal(Math.pow(10, 38)))<0){ //num is small, so use the standard Math.log() function
			if (num.compareTo(new BigDecimal("0.00001"))<0)
				sum = new BigDecimal("0.0");
			else
				sum = new BigDecimal(Math.log(num.doubleValue()));
		}
		else { //num is large, so compute an approximation to the natural logarithm

			double t = 0.0;

			boolean done = false;
			while (!done){
				if (x.compareTo(exp)>0){
					t += 1.0;
				}
				else {
					sum = sum.add(new BigDecimal(t+Math.log(x.doubleValue())));
					done = true;
				}
				x = x.divide(exp,4);
			}
		}

		return sum;
	}
	
	// Returns the natural logarithm of a big decimal
	public double logToDouble(BigDecimal num){
		double eDoub = 2.71828182845904523536;
		double log10ofE = Math.log10(eDoub);
		if (num.compareTo(new BigDecimal("0.0"))<0){ //num is negative
			throw new IllegalArgumentException("log of a negative number: " + num);
		}

		int powerOfTen = num.round(mc).scale() * -1;
		double fract = num.movePointLeft(powerOfTen).doubleValue();
		double returnVal = (double) (powerOfTen + Math.log10(fract)) / log10ofE;

		//double realVal = Math.log(num.doubleValue());
		//System.out.println(returnVal-realVal);
		return returnVal;

	}
	
	// Returns the logarithm in base 10 of a BigDecimal.
	public double log10(BigDecimal num){
		if (num.compareTo(new BigDecimal("0.0"))<0){ //num is negative
			throw new IllegalArgumentException("log of a negative number: " + num);
		}

		int powerOfTen = num.round(mc).scale() * -1;
		double fract = num.movePointLeft(powerOfTen).doubleValue();
		double returnVal = (double) (powerOfTen + Math.log10(fract));

		return returnVal;

	}
}
