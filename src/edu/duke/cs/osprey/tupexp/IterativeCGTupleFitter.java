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

package edu.duke.cs.osprey.tupexp;

import java.util.ArrayList;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.ConjugateGradient;
import org.apache.commons.math3.linear.RealLinearOperator;
import org.apache.commons.math3.linear.RealVector;

/**
 *
 * Modified least-squares fitter for LUTE coefficients
 * Uses iterations of conjugate gradient
 * 
 * @author mhall44
 */
public class IterativeCGTupleFitter extends CGTupleFitter {
    
    /*RealLinearOperator AtA;
    RealVector Atb;
    int numSamp, numTup;
    TupleIndexMatrix tupIndMat;
    ArrayList<int[]> samples;*/
    ArrayList<double[]> goodRegionBounds;//bounds on the "good" region of fit values
    //for each sample
    
    double[] curFitVals = null;//value of fit at each sample before current iteration
    RealVector curCoeffs = null;
    double curResid = Double.POSITIVE_INFINITY;
    
    //Initial fit: A_samp,tup = 1_tupinsamp * 1_sampoverbounds
    //Hence AtA = sum_samp_in_bounds 1_tup_in_samp 1_tup2_in_samp
    //and Atb = sum_samp (tup_in_samp * samp_constr_broken)
    //indeed any fit is that way  hahahaha
    
    double damperLambda = 1e-4;//Some tuples will have no equality, only inequality
    //restraints, and thus may revert to 0 when moved to the 
    //inactive restraint set, preventing convergence.  A slight penalty on coeff changes
    //will prevent this
    
    
    
    public IterativeCGTupleFitter(TupleIndexMatrix tim, ArrayList<int[]> samp, int numTuples, ArrayList<double[]> goodRegionBounds){
        //We'll fit the specified (sample,trueVal) pairs to an expansion in the tuples in tim
        
        samples = samp;
        numSamp = samples.size();
        numTup = numTuples;
        tupIndMat = tim;
        
        this.goodRegionBounds = goodRegionBounds;        
        
        
        AtA = new RealLinearOperator(){

            @Override
            public int getRowDimension() {
                return numTup;
            }

            @Override
            public int getColumnDimension() {
                return numTup;
            }

            @Override
            public RealVector operate(RealVector rv) throws DimensionMismatchException {
                return applyAtA(rv);
            }
        };
    }
    
    
    double[] calcFitVals(RealVector rv){
        double[] fitVals = new double[numSamp];
        for(int s=0; s<numSamp; s++){
            //if(isSampleRestrained(s)){//we need these for validation, if not for AtA aplication
                ArrayList<Integer> sampTup = tupIndMat.calcSampleTuples(samples.get(s));
                for(int t : sampTup)
                    fitVals[s] += rv.getEntry(t);
            //}
        }
        
        return fitVals;
    }
    
    
    RealVector applyAtA(RealVector rv){
        //first apply A
        double Arv[] = calcFitVals(rv);

        //then apply A^T to Arv
        double ans[] = new double[numTup];
        for(int s=0; s<numSamp; s++){
            if(isSampleRestrained(s)){
                ArrayList<Integer> sampTup = tupIndMat.calcSampleTuples(samples.get(s));
                for(int t : sampTup)
                    ans[t] += Arv[s];
            }
        }
        
        //damping
        if(curCoeffs!=null){//not first iteration
            for(int t=0; t<numTup; t++)
                ans[t] += damperLambda * rv.getEntry(t);
        }

        return new ArrayRealVector(ans,false);//make RealVector without copying ans
    }
    
    
    
    
    
    double getCurTarget(int s){//get current "target" value (for iteration of modified lsq)
        //for the given sample.  NaN if not applicable (restraints inactive for sample)
        double bounds[] = goodRegionBounds.get(s);
        if(curFitVals==null){//First iteration.  Assume
            //samples with any leeway don't have a target
            if(bounds[0]==bounds[1])
                return bounds[0];
            else
                return Double.NaN;
        }
        else {
            double curFitVal = curFitVals[s];
            if(curFitVal<bounds[0])//lower restraint active
                return bounds[0];
            else if(curFitVal>bounds[1])//upper
                return bounds[1];
            else//neither
                return Double.NaN;
        }
    }
    
    
    boolean isSampleRestrained(int s){
        double bounds[] = goodRegionBounds.get(s);
        if(curFitVals==null){//First iteration.  Assume
            //samples with any leeway don't have a target
            return (bounds[0]==bounds[1]);
        }
        else {
            double curFitVal = curFitVals[s];
            return (curFitVal<bounds[0] || curFitVal>bounds[1]);
        }
    }
    
    
    RealVector calcRHS(){//Calculate right-hand side vector of normal equations
        double atb[] = new double[numTup];
        //apply A^T to true vals
        for(int s=0; s<numSamp; s++){
            double curTarget = getCurTarget(s);
            if(!Double.isNaN(curTarget)){//restraint active for sample
                ArrayList<Integer> sampTup = tupIndMat.calcSampleTuples(samples.get(s));
                for(int t : sampTup)
                    atb[t] += curTarget;
            }
        }
        
        //damping.  Slightly penalizes changes from curCoeffs
        if(curCoeffs!=null){
            for(int t=0; t<numTup; t++)
                atb[t] += damperLambda * curCoeffs.getEntry(t);
        }
        
        Atb = new ArrayRealVector(atb);
        return Atb;
    }
    
    
    double calcResidual(double[] fitVals){
        double resid=0;
        
        for(int s=0; s<numSamp; s++){
            double dev = 0;
            double bounds[] = goodRegionBounds.get(s);
            
            if(fitVals[s]<bounds[0])
                dev = fitVals[s] - bounds[0];
            else if(fitVals[s]>bounds[1])
                dev = fitVals[s] - bounds[1];
            
            resid += dev*dev;
        }
        
        return resid/numSamp; 
    }
    
    
    boolean checkDone(double[] oldFitVals, double[] newFitVals){
        //see if we crossed any restraint boundaries between two consecutive iterations
        //(by a numerically significant amount)
                        
        
        if(oldFitVals==null)//first iteration, can't be done (no comparison to do)
            return false;
        
        double tol = 1e-6;
        
        for(int s=0; s<numSamp; s++){
            double bounds[] = goodRegionBounds.get(s);
            if(bounds[0]<bounds[1]){//there are boundaries for this sample
                //(i.e. not pure quadratic objective function term)
                if(oldFitVals[s]<bounds[0]){//check if lower restraint disappeared
                    if(newFitVals[s]>bounds[0]+tol)
                        return false;
                }
                else if(oldFitVals[s]>bounds[1]){//upper
                    if(newFitVals[s]<bounds[1]-tol)
                        return false;
                }
                else {//no restraint now...see if one appeared
                    if(newFitVals[s]<bounds[0]-tol || newFitVals[s]>bounds[1]+tol)
                        return false;
                }
            }
        }
        
        return true;
    }
    
    
    @Override
    double[] doFit(){
        //return fit tuple coefficients
        
        ConjugateGradient cg = new ConjugateGradient(100000,1e-6,false);//max_iter; delta; whether to check pos def
        //delta is target ratio of residual norm to true vals norm
        
        long startTime = System.currentTimeMillis();
        
        
        while(true){
            double iterStartTime = System.currentTimeMillis();
            
            Atb = calcRHS();
            RealVector ans = cg.solve(AtA, Atb);
            double[] newFitVals = calcFitVals(ans);
            
            System.out.println( "Conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-iterStartTime) );
            
            //boolean done = checkDone(curFitVals, newFitVals);
            double resid = calcResidual(newFitVals);
            System.out.println("Step residual: "+resid);
            
            if(resid>curResid){//gotten worse...use previous vals
                System.out.println( "Iterative conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );
                return curCoeffs.toArray();
            }
            else if(resid>curResid-1e-4){//basically converged
                System.out.println( "Iterative conjugate gradient fitting time (ms): " + (System.currentTimeMillis()-startTime) );
                return ans.toArray();
            }
            else{//keep going
                curCoeffs = ans;
                curFitVals = newFitVals;
                curResid = resid;
            }
        }
    }
}
