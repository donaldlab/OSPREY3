/*
 * TODO:
 * 1) implement some simulated annealing / monte carlo steps to find better minima
 * 2) maybe use levenberg-marquardt to update the jacobian matrix? pg 14 of
 * Improvements to the Levenberg-Marquardt algorithm for nonlinear least-squares minimization
 */

package edu.duke.cs.osprey.minimization;

import java.util.ArrayList;
import java.util.Arrays;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.math.Functions;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class BFGSMinimizer extends CCDMinimizer implements Minimizer {

	DoubleMatrix1D g = null;//grad objFcn at x
	DoubleMatrix2D hessian = null;
	ArrayList<Integer> hessIdx = null; //indices for hessian rows

	DoubleMatrix1D xnew = null;
	DoubleMatrix1D xi = null; // new direction
	DoubleMatrix1D dg = null;
	DoubleMatrix1D hdg = null;
	DoubleMatrix1D zero = null;

	final double EPS = Math.pow(2, -7);
	final double TOLX = 100 * EPS;
	final double STPMX = 100.0;

	int n = 0; // degrees of freedom
	double sum, stpmax, den, fac, fae, sumdg, sumxi, hij;

	final double gtol = Math.pow(2, -7);

	double fp = Double.POSITIVE_INFINITY;
	double fret = Double.POSITIVE_INFINITY;
	
	int maxIter = 30;


	public BFGSMinimizer(ObjectiveFunction ofn, boolean useCorners) {
		super(ofn, useCorners);
		setMaxIter(numIter);
	}

	
	public void setMaxIter( int val ) {
		maxIter = val;
	}
	
	
	public int getMaxIter() {
		return maxIter;
	}
	

	public DoubleMatrix1D minimize() {

		long minStartTime = System.currentTimeMillis();

		// update constraints at each iteration?????
		DoubleMatrix1D constr[] = objFcn.getConstraints();
		DOFmin = constr[0];
		DOFmax = constr[1];

		if(!compInitVals())//Initialize x
			return null;//No initial values found (this is currently only for IBEX)

		firstStep = new double[numDOFs];
		lastStep = new double[numDOFs];
		Arrays.fill(firstStep, 1);
		Arrays.fill(lastStep, 1);

		objFcn.setDOFs(x);

		initBFGS();
		for( int it = 0; it < getMaxIter(); it++ ) {
			gradLnSrch(it);
			rippleJumper(xnew);
			fret = func(xnew);
			boolean converged = BFGSIter(it);
			if(converged) break;
		}

		minTime = System.currentTimeMillis() - minStartTime;
		//System.out.println(minTime);

		//System.out.println("Iter: " + its + " E: " + objFcn.getValue(x));

		return x;
	}


	protected boolean BFGSIter(int it) {

		//System.out.println("Iter: " + it + " Eold: " + fp + " Enew: " + fret);

		// update the line direction
		xi = subtract(xnew, x);

		// update the current point
		x.assign(xnew);
		objFcn.setDOFs(x);

		// test for convergence on delta x
		double test = 0.0;
		for( int i = 0; i < n; ++i ) {
			double temp = Math.abs(xi.getQuick(i)) / Math.max(Math.abs(x.getQuick(i)), 1.0);
			if( temp > test ) test = temp;
		}
		if( Math.abs(fp - fret) < EConvTol && test < TOLX )
			return true;

		// save old gradient
		dg = g;

		// get new gradient
		g = funcDF();

		// test for convergence on zero gradient
		test = 0.0;
		den = Math.max(fret, 1.0);
		for( int i = 0; i < n; ++i ) {
			double temp = Math.abs(g.getQuick(i)) * Math.max(Math.abs(x.getQuick(i)), 1.0) / den;
			if( temp > test ) test = temp;
		}
		if( Math.abs(fp - fret) < EConvTol && test < gtol )
			return true;

		fp = fret;

		// compute difference of gradients
		dg = subtract(g, dg);

		// compute gradient difference * hessian
		hdg = hessian.zMult(dg, hdg);

		// calculate dot products for denominators
		fac = dg.zDotProduct(xi);
		fae = dg.zDotProduct(hdg);
		sumdg = dg.zDotProduct(dg);
		sumxi = xi.zDotProduct(xi);

		// skip update if fac not sufficiently positive
		if( fac > Math.sqrt(EPS*sumdg*sumxi) ) {
			fac = 1.0 / fac;
			final double fad = 1.0 / fae;

			for( int i = 0; i < n; ++i ) dg.setQuick( i, fac*xi.getQuick(i)-fad*hdg.getQuick(i) );

			// bfgs hessian updating formula
			/*
			for( int i = 0; i < n; ++i ) for ( int j = i; j < n; ++j ) {
				hij = hessian.getQuick(i,j);

				hessian.setQuick( i, j, hij + fac*xi.getQuick(i)*xi.getQuick(j)
						-fad*hdg.getQuick(i)*hdg.getQuick(j)
						+fae*dg.getQuick(i)*dg.getQuick(j) );

				hessian.setQuick( j, i, hessian.getQuick(i,j) );
			}
			*/
			
			hessIdx.parallelStream().forEach((i) -> {
				for ( int j = i; j < n; ++j ) {
					hij = hessian.getQuick(i,j);

					hessian.setQuick( i, j, hij + fac*xi.getQuick(i)*xi.getQuick(j)
							-fad*hdg.getQuick(i)*hdg.getQuick(j)
							+fae*dg.getQuick(i)*dg.getQuick(j) );

					hessian.setQuick( j, i, hessian.getQuick(i,j) );
				}
			});
			
		}

		// calculate the next direction to follow
		for( int i = 0; i < n; ++i ) {
			xi.setQuick(i, 0.0);

			for( int j = 0; j < n; ++j ) {
				xi.setQuick(i, xi.getQuick(i) - hessian.getQuick(i, j)*g.getQuick(j) );
			}
		}

		if(rescalingGC != null)
			rescaleValues();

		return false;
	}


	protected double func(DoubleMatrix1D input) {
		return objFcn.getValue(input);
	}


	protected void initBFGS() {
		n = numDOFs;

		// get function value at x    
		fp = func(x);

		// get gradient
		g = funcDF();

		dg = x.like();
		hdg = x.like();
		xi = x.like();
		xnew = x.like();
		zero = x.like();

		hessian = DoubleFactory2D.dense.identity(n);

		// set indices for hessian row indexes
		hessIdx = new ArrayList<>();
		for(int i = 0; i < n; ++i) hessIdx.add(hessIdx.size());
		hessIdx.trimToSize();

		xi = subtract(zero, g);
		sum = x.zDotProduct(x);

		stpmax = STPMX * Math.max(Math.sqrt(sum), n);
	}


	protected void gradLnSrch(int it) {

		double ALF = 1e-4, a, alam, alam2=0.0, alamin, b, disc, f2=0.0;
		double rhs1, rhs2, slope=0.0,sum=0.0,temp, test, tmplam;
		int i;
		@SuppressWarnings("unused")
		boolean check = false;

		sum = Math.sqrt( xi.zDotProduct(xi) );
		if( sum > stpmax ) {
			for( i = 0; i < n; ++i ) xi.setQuick(i, xi.getQuick(i)*stpmax/sum);
		}

		slope = g.zDotProduct(xi);
		if( Double.isInfinite(sum) || Double.isNaN(sum) || Double.isNaN(slope) || slope >= 0.0 ) {

			// we have reached limits of numerical precision, 
			// so bfgs is not useful in this iteration
			xnew.assign(x);
			fret = fp;
			return;

			// throw new RuntimeException("ERROR: roundoff problem in gradlnsrch");
		}

		// compute lambda_min
		test = 0.0;
		for(i = 0; i < n; ++i) {
			temp = Math.abs(xi.getQuick(i)) / Math.max(Math.abs(x.getQuick(i)), 1.0);
			if(temp > test) test = temp;
		}

		alamin = TOLX/test;
		alam = 1.0;

		// always try the full newton step first
		// try iteration of loop
		while( true ) {

			xnew = addProd(x, xi, alam);

			// make sure xnew is within bounds
			for( int dof = 0; dof < n; ++dof ) {
				double dof_base = xnew.get(dof);
				if( isOutOfRange(dof_base,dof) )
					xnew.setQuick(dof, getEdgeDOFVal(dof_base,dof));
			}

			fret = func(xnew);

			if( alam < alamin ) {

				// convergence on delta x. for zero finding, the calling function
				// should verify the convergence
				if(fret > fp) {
					fret = fp;
					xnew.assign(x);
				}
				check = true;
				return;
			} else if( fret <= fp+ALF*alam*slope ) {
				//sufficient function decrease
				return;
			} else {
				// backtrack
				if( alam == 1.0 ) {
					// first time
					tmplam = -slope/(2.0*(fret-fp-slope));
				} else {
					// subsequent backtracks
					rhs1 = fret - fp - alam * slope;
					rhs2 = f2 - fp - alam2 * slope;
					a = (rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
					b = (-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
					if(a == 0.0) tmplam = -slope/(2.0*b);
					else {
						disc = b*b-3.0*a*slope;
						if(disc < 0.0) tmplam = 0.5*alam;
						else if(b <= 0.0) tmplam = (-b + Math.sqrt(disc))/(3.0*a);
						else tmplam = -slope/(b + Math.sqrt(disc));
					}

					if(tmplam > 0.5 * alam)
						tmplam = 0.5 * alam;

				}
			}

			alam2 = alam;
			f2 = fret;
			alam = Math.max(tmplam, 0.1*alam);

		}
	}


	protected DoubleMatrix1D funcDF() {

		DoubleMatrix1D grad = x.like();

		for(int dof=0; dof<numDOFs; dof++)
			grad.set(dof, partDerivDx(dof));

		return grad;
	}


	protected double partDeriv2Dx(int dof) {

		double x_dof = x.get(dof);

		double temp = x_dof;
		double dx = EPS*Math.abs(temp); 
		if(dx == 0.0) dx = EPS;
		// trick to reduce finite precision error
		x_dof = temp+dx;
		dx=x_dof-temp;
		x_dof = temp;

		//Values up ( f(x+dx) ) or down ( f(x-dx) ) a step for this DOF: initialized to inf (indicating infeasibility)
		//in case the DOF value is out of range
		double f_xpdx = Double.POSITIVE_INFINITY;
		double f_xmdx = Double.POSITIVE_INFINITY;

		if( ! isOutOfRange(x_dof+dx,dof) ) {
			f_xpdx = objFcn.getValForDOF(dof,x_dof+dx);
			// restore dof value, since getvalfordof is not read only
			objFcn.setDOF(dof, x_dof);
		}

		if( ! isOutOfRange(x_dof-dx,dof) ) {
			f_xmdx = objFcn.getValForDOF(dof,x_dof-dx);
			// restore dof value, since getvalfordof is not read only
			objFcn.setDOF(dof, x_dof);
		}

		return (f_xpdx-f_xmdx)/(2.0*dx);
	}


	protected double partDerivDx(int dof) {
		double x_dof = x.get(dof);
		double fx_dof = objFcn.getValForDOF(dof,x_dof);

		double dx = EPS*Math.abs(x_dof); 
		if(dx == 0.0) dx = EPS;

		double xpdx = x_dof + dx;

		//Value up ( f(x+dx) ) a step for this DOF: initialized to inf (indicating infeasibility)
		double f_xpdx = Double.POSITIVE_INFINITY;

		if( ! isOutOfRange(xpdx,dof) ) {
			f_xpdx = objFcn.getValForDOF(dof,xpdx);
			// restore dof value, since getvalfordof is not read only
			// objFcn.setDOF(dof, x_dof);
		}

		return (f_xpdx-fx_dof)/dx;
	}


	//Convenience functions
	static DoubleMatrix1D add(DoubleMatrix1D a, DoubleMatrix1D b){
		return a.copy().assign(b, Functions.plus);
	}

	static DoubleMatrix1D addProd(DoubleMatrix1D a, DoubleMatrix1D b, double coeff){
		// Return a + coeff*b
		return b.copy().assign( Functions.mult(coeff) ).assign( a, Functions.plus );
	}

	static DoubleMatrix1D subtract(DoubleMatrix1D a, DoubleMatrix1D b){
		return a.copy().assign(b, Functions.minus);
	}


	protected void rippleJumper(DoubleMatrix1D vector) {
		for(int dof = 0; dof < n; ++dof) {
			rippleJumper(vector, dof);
		}
	}

	protected void rippleJumper(DoubleMatrix1D vector, int dof) {
		//RIPPLE JUMPER!!
		//search up and down some, with the goal of being robust for rugged surfaces
		//optimized for dihedrals...
		double downRipple = vector.get(dof)-1;
		double upRipple = vector.get(dof)+1;
		double curVal = objFcn.getValForDOF(dof, vector.get(dof));
		if( ! isOutOfRange( downRipple, dof) ){
			if( objFcn.getValForDOF(dof,downRipple) < curVal )
				vector.set(dof, downRipple);
		}
		if( ! isOutOfRange( upRipple, dof ) ){
			if( objFcn.getValForDOF(dof,upRipple) < curVal )
				vector.set(dof, upRipple);
		}
		//END RIPPLE JUMPER

		// objFcn.setDOF(dof, vector.get(dof));
	}
}
