package edu.duke.cs.osprey.minimization;

import java.util.Arrays;

import cern.colt.matrix.DoubleMatrix1D;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
public class HBFGSCCDMinimizer extends BFGSMinimizer implements Minimizer {

	public HBFGSCCDMinimizer(ObjectiveFunction ofn, boolean useCorners) {
		super(ofn, useCorners);
		// TODO Auto-generated constructor stub
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

		fp = Double.POSITIVE_INFINITY;
		fret = Double.POSITIVE_INFINITY;
		objFcn.setDOFs(x);

		/*
		// ccd loop
		for( int it = 0; it < numIter; it++ ) {
			double fret = objFcn.getValue(x);
			if(Double.isInfinite(fret) && Double.isInfinite(fp)) {
				if(useRandMinCheck) fret = doRandMinCheck(fret);
				if(Double.isInfinite(fret)) break;
			}
			//Considered to be converged if energy improvement is small enough
			if( fp - fret < EConvTol ){//or if numerical error is causing the energy to rise
				if(useRandMinCheck) fret = doRandMinCheck(fret);
				if( fp - fret < EConvTol ) break;
			}
			fp = fret;
			dihedralsLnSrch(it, x);
			if(rescalingGC != null) rescaleValues();
		}
		 */

		//bfgs loop
		initBFGS();
		for( int it = 0; it < getMaxIter(); it++ ) {
			dihedralsLnSrch(it, x);
			gradLnSrch(it);
			boolean converged = BFGSIter(it);
			if(converged) break;
		}

		minTime = System.currentTimeMillis() - minStartTime;
		//System.out.println(minTime);

		//System.out.println("Iter: " + its + " E: " + objFcn.getValue(x));

		return x;
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


	protected void dihedralsLnSrch(int it, DoubleMatrix1D vector) {
		for( int dof = 0; dof < numDOFs; ++dof ) 
			DOFLnSrch(it, dof, vector);
	}


	protected void DOFLnSrch(int it, int dof, DoubleMatrix1D vector) {
		// this code is heavily borrowed from mark's ccd minimizer

		//Don't minimize for this DOF if it is constrained to a single value
		if( DOFmax.get(dof) == DOFmin.get(dof) )
			return;

		//Do line search for this dof

		double dof_base = vector.get(dof);
		double curVal = objFcn.getValForDOF(dof,dof_base);//Get value for DOF, shifted up or down as needed
		double step = getStepSize(dof,it);//Step size for the given degree of freedom, adjusted for iteration


		//Values up or down a step for this DOF: initialized to inf (indicating infeasibility)
		//in case the DOF value is out of range
		double upVal = Double.POSITIVE_INFINITY;
		double downVal = Double.POSITIVE_INFINITY;

		if( ! isOutOfRange(dof_base+step,dof) )
			upVal = objFcn.getValForDOF(dof,dof_base+step);
		if( ! isOutOfRange(dof_base-step,dof) )
			downVal = objFcn.getValForDOF(dof,dof_base-step);


		double a = (upVal+downVal-2*curVal)/(2*step*step);
		double estmin=0;


		//TRYING DIFFERENT VERSIONS
		if( ( a <= 0 ) || Double.isNaN(a) || Double.isInfinite(a) ){
			//negative a is not a good sign...use the lesser of upVal, downVal
			//infinite or nan a means we're hitting a constraint or impossible conformation
			if( upVal < downVal )
				estmin = dof_base+step;
			else
				estmin = dof_base-step;
		}

		else{
			double b = (upVal-curVal)/step-a*step;
			estmin = dof_base-b/(2*a);//2*a*(dof_val-dof_base)+b=0 here, with positive a
		}


		if( isOutOfRange(estmin,dof) )
			estmin = getEdgeDOFVal(estmin,dof);


		//if( Math.abs(estmin-dof_base) < numTol )//This can happen if both upVal and downVal are infinite (perhaps due to a loop closure failure)
		//    break;

		double estminVal = objFcn.getValForDOF(dof,estmin);
		double estminValOld = curVal;

		if(estminVal < curVal){
			while(true) {//Can break either on hitting a constraint or on estminVal starting to increase
				estmin = dof_base + 2*(estmin-dof_base);

				if( isOutOfRange(estmin,dof) ){
					double edge = getEdgeDOFVal(estmin,dof);

					double edgeVal = objFcn.getValForDOF(dof,edge);

					//137DEBUG!!
					GVCountEdge++;

					if(edgeVal<estminVal)
						x.set(dof, edge);
					else
						x.set(dof, dof_base+0.5*(estmin-dof_base) );

					break;
				}

				estminValOld = estminVal;
				estminVal = objFcn.getValForDOF(dof,estmin);

				//137DEBUG!!
				GVCountBigger++;

				double tol = numTol * Math.max(1,Math.abs(estminVal));

				if( !(estminVal < estminValOld - tol) ){//No noticeable improvement in the last step
					x.set(dof,dof_base+0.5*(estmin-dof_base));
					break;
				}
			}
		}
		else if(estminVal>curVal+numTol) {//need to backtrack
			//won't hit a constraint 
			while(true) {//Can break on estminVal starting to increase, or decreasing negligibly
				estmin = dof_base + 0.5*(estmin-dof_base);

				estminValOld = estminVal;
				estminVal = objFcn.getValForDOF(dof,estmin);


				//137DEBUG!!
				GVCountSmaller++;

				double tol = numTol * Math.max(1,Math.abs(estminVal));
				//need to avoid getting stuck at the same estmin, with numerically trivial improvements...

				if( estminValOld < estminVal + tol ){//No significant improvement in the last step
					if(estminValOld<curVal)//have improvement over curVal at least
						x.set(dof,dof_base+2*(estmin-dof_base));
					break;
				}
			}
		}


		lastStep[dof] = vector.get(dof) - dof_base;
		if(it==0)
			firstStep[dof] = lastStep[dof];


		//RIPPLE JUMPER!!
		//search up and down some, with the goal of being robust for rugged surfaces
		//optimized for dihedrals...
		double downRipple = vector.get(dof)-1;
		double upRipple = vector.get(dof)+1;
		curVal = objFcn.getValForDOF(dof, vector.get(dof));
		if( ! isOutOfRange( downRipple, dof) ){
			if( objFcn.getValForDOF(dof,downRipple) < curVal )
				vector.set(dof, downRipple);
		}
		if( ! isOutOfRange( upRipple, dof ) ){
			if( objFcn.getValForDOF(dof,upRipple) < curVal )
				vector.set(dof, upRipple);
		}
		//END RIPPLE JUMPER

		objFcn.setDOF(dof, vector.get(dof));
	}

}
