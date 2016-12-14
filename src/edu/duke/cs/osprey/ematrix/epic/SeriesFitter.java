package edu.duke.cs.osprey.ematrix.epic;

/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
 */

///////////////////////////////////////////////////////////////////////////////////////////////
//	SeriesFitter.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  MAH           Mark A. Hallen	  Duke University               mah43@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import cern.jet.math.Functions;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.TreeSet;
import org.apache.commons.math3.linear.ConjugateGradient;


//This class contains methods to fit a series to data
//Used to find low-order Taylor series approximations for the energy in a voxel
public class SeriesFitter {
    
    
    static boolean useLineSearch = false;//generally doesn't provide significant benefit,
    //and numerically questionable.  


    public static double[] fitSeries(DoubleMatrix1D[] samp, double trueVals[], double weights[],
            double lambda, boolean includeConst, int order){

        return fitSeries(samp,trueVals,weights,lambda,includeConst,order,order,null,false,null,null);
    }


    static double[] fitSeries(DoubleMatrix1D[] samp, double trueVals[], double weights[],
            double lambda, boolean includeConst, int order, int PCOrder, boolean isPC[],
            boolean update, DoubleMatrix1D c, DoubleMatrix2D M){
            //Default: quadratic with 0 constant term
        //includeConst adds the constant term and quartic adds 3rd-4th terms
        
        //If(update), then we are updating a previous fit with the same samp, trueVals
        //but different weights (and weights contains the difference in weights in this case)
        //c and M should be null unless we're updating or making a fit that we will update later

        long startTime = System.currentTimeMillis();

        int nd = samp[0].size();//number of degrees of freedom
        int numParams = getNumParams(nd,includeConst,order);

        if(PCOrder > order){//add in parameters for PC orders
            int numPCs = countTrue(isPC);

            for(int n=order+1; n<=PCOrder; n++)
                numParams += getNumParamsForOrder(numPCs,n);
        }
        

        int numSamples = samp.length;

        System.out.println("Fit has "+numSamples+" samples for "+numParams+" parameters");


        if(c==null)
            c = DoubleFactory1D.dense.make(numParams);
        if(M==null)
            M = DoubleFactory2D.dense.make(numParams,numParams);


        //scratch for terms of c, M
        DoubleMatrix1D cScratch = DoubleFactory1D.dense.make(numParams);
        DoubleMatrix2D MScratch = DoubleFactory2D.dense.make(numParams,numParams);


        if(!update){
            for(int p=0; p<numParams; p++)//deal with regularization
                M.set(p, p, lambda);
        }


        //long setupDoneTime = System.currentTimeMillis();
        //System.out.println("fitSeries setup time (ms): "+(setupDoneTime-startTime));

        for(int s=0; s<numSamples; s++){//summing each of these terms over the samples

            if( (!update) || (weights[s]!=0) ) {
                //if updating, can ignore unchanged weights (weights[s]==0)

                double weight = 1;//weight for samples
                if(weights!=null)
                    weight = weights[s];


                //we want cScratch^T params = trueVals in the best-fit sense
                //Least-squares equations are then (sum_samples weight * MScratch) * params = sum_samples weight * c_scratch *trueVals
                //where MScratch is cScratch*cScratch^T
                calcSampParamCoeffs(cScratch,samp[s],nd,includeConst,order,PCOrder,isPC);

                Algebra.DEFAULT.multOuter(cScratch, cScratch, MScratch);
                MScratch.assign(Functions.mult(weight));
                M.assign(MScratch, Functions.plus);

                cScratch.assign(Functions.mult(trueVals[s]*weight));
                c.assign(cScratch, Functions.plus);
            }
        }


        DoubleMatrix2D C = DoubleFactory2D.dense.make(c.toArray(),numParams);//c as a column vector


        //long MCTime = System.currentTimeMillis();
        //System.out.println("fitSeries M, C calc time (ms): "+(MCTime-setupDoneTime));

        double[] v = null;//the best-fit parameters, as a column vector
        try {
            v = Algebra.DEFAULT.solve(M,C).viewColumn(0).toArray();
        }
        catch(IllegalArgumentException e){//indices singular M
            //solve in a way robust to singularities
            //basically pseudoinverse(M)*C
            SingularValueDecomposition svd = new SingularValueDecomposition(M);
            //M = U*S*V', so M^-1 = V*invS*U'
            DoubleMatrix2D invS = svd.getS().copy();
            
            //this tolerance is from SingularValueDecomposition.java (used there to compute rank)
            //here we use the special case that M is square
            double eps = Math.pow(2.0,-52.0);
            double tol = invS.rows()*invS.get(0,0)*eps;
            
            for(int i=0; i<invS.rows(); i++){
                double singVal = invS.get(i,i);
                if(singVal>tol)
                    invS.set(i, i, 1./singVal);
                else
                    invS.set(i, i, 0);
            }
            
            DoubleMatrix2D ansCol = Algebra.DEFAULT.mult( Algebra.DEFAULT.mult(
                    Algebra.DEFAULT.mult( svd.getV(), invS ), Algebra.DEFAULT.transpose(svd.getU())), C );
            
            
            //DEBUG!!!!
            //DoubleMatrix2D checkC = Algebra.DEFAULT.mult(M, ansCol);
            
            v = ansCol.viewColumn(0).toArray();
        }

        //long vTime = System.currentTimeMillis();
        //System.out.println("fitSeries matrix solution time (ms): "+(vTime-MCTime));



        //Checking fit accuracy (done in fitSeriesIterative if updating)
        if(!update){//Calculating the mean residual like this doesn't make sense for an update
            //since we don't have the actual weights
            double meanResidual = 0;
            double weightSum = 0;

            for(int s=0; s<numSamples; s++){
                double bv = evalSeries(v,samp[s],nd,includeConst,order,PCOrder,isPC);


                double weight = 1;//weight for samples
                if(weights!=null)
                    weight = weights[s];

                double residTerm = (trueVals[s]-bv)*(trueVals[s]-bv);

                if(Double.isInfinite(residTerm) || Double.isNaN(residTerm)){
                    System.err.println("ERROR: SeriesFitter.fitSeries gives infinite residual term: "+residTerm);

                    System.out.print("Sample: ");
                    for(int dof=0; dof<nd; dof++)
                        System.err.print(samp[s].get(dof)+" ");
                    System.out.println();
                    
                    System.err.println(" TRUEVAL="+trueVals[s]+" BV="+bv);
                    System.err.println("params:");
                    for(double param : v)
                        System.out.println(param);

                    throw new RuntimeException("Infinite or nan residual");
                }

                //if want sample-by-sample output
                //System.out.println("TRAININGSET TRUE:"+trueVals[s]+" G2:"+bv);
                
                meanResidual += weight*residTerm;
                weightSum += weight;
            }

            meanResidual /= weightSum;
            System.out.println("TRAINING SET MEAN RESIDUAL:"+meanResidual);
        }


        long doneTime = System.currentTimeMillis();
        //System.out.println("fitSeries checking time (ms): "+(doneTime-vTime));
        System.out.println("fitSeries time (ms): "+(doneTime-startTime));

        return v;
    }


    //Fit series with one-sided constraints
    //The idea here is that if we have a set of points P such that the fit
    //and true value are >= bCutoff for points in P and the fit is an optimal least-squares
    //fit with P excluded, then the fit is locally optimal with one-sided constraints
    //because the gradient of the objective function with respect to the coefficients
    //is the same in the two-sided fit with P excluded and the full one-sided fit
    //we also support bCutoff2, to have lower-bounding but not upper-bounding restraints
    //This is intended for EPIC
    static double[] fitSeriesIterative(DoubleMatrix1D[] samp, double trueVals[], double weights[],
            double lambda, boolean includeConst, int order, double bCutoffs[], double bCutoffs2[],
            int PCOrder, boolean isPC[]){
        
        long startTime = System.currentTimeMillis();
        System.out.println("Starting fitSeriesIterative...");


        int numSamples = samp.length;
        int nd = samp[0].size();
        
        if(bCutoffs.length==1){//single bCutoff to be used
            double bCutoff = bCutoffs[0];
            bCutoffs = new double[numSamples];
            Arrays.fill(bCutoffs,bCutoff);
        }
        //now bCutoffs has a cutoff for each sample
        
        int numParams = getNumParams(nd,includeConst,order);

        if(PCOrder > order){//add in parameters for PC orders
            int numPCs = countTrue(isPC);

            for(int n=order+1; n<=PCOrder; n++)
                numParams += getNumParamsForOrder(numPCs,n);
        }


        //now set up the data for the iterative fits
        //samples can be turned on and off using fitWeights
        //we will need to create two entries for samples with trueVals between bCutoff and bCutoff2
        //since they may be turned on either to penalize deviation from bCutoff or from trueVal
        
        //first entry for each entry will penalize deviation from bCutoff if trueVal>=bCutoff
        //secondEntries are for trueVals between bCutoff and bCutoff2
        ArrayList<Integer> secondEntries = new ArrayList<>();//list of samples needing second entry
        HashMap<Integer,Integer> revSecondEntries = new HashMap<>();//reverse lookup
        for(int s=0; s<numSamples; s++){
            if( (trueVals[s]>=bCutoffs[s]) && (trueVals[s]<bCutoffs2[s]) ){//as in isRestraintActive
                revSecondEntries.put(s,secondEntries.size());
                secondEntries.add(s);
            }
        }
        
        int numRestraints = numSamples+secondEntries.size();
        
        //data for basic least-squares fits
        DoubleMatrix1D[] fitSamp = new DoubleMatrix1D[numRestraints];
        double fitTrueVals[] = new double[numRestraints];
        double fitWeights[] = new double[numRestraints];
        
        
        for(int s=0; s<numSamples; s++){//"normal" entries
            fitSamp[s] = samp[s];
            
            if(trueVals[s]>=bCutoffs[s]){
                fitWeights[s] = 0;
                fitTrueVals[s] = bCutoffs[s];
            }
            else {
                fitWeights[s] = weights[s];
                fitTrueVals[s] = trueVals[s];
            }
        }
        for(int s2=0; s2<secondEntries.size(); s2++){
            fitSamp[numSamples+s2] = samp[secondEntries.get(s2)];
            fitWeights[numSamples+s2] = 0;
            fitTrueVals[numSamples+s2] = trueVals[secondEntries.get(s2)];
        }
        
        
        //Initial guess of set P is all points with trueVals[s] >= bCutoff
        //that is, all points that have possible series values that make the restraint inactive

        boolean done = false;
        double coeffs[] = null;
        double meanResidual=0, weightSum=0;
        double prevResid = Double.POSITIVE_INFINITY;
        double oldCoeffs[] = null;

        double oldSerVals[] = new double[numSamples];//values of series at each sample, for previous iteration
        //preallocating to all infinity because all trueVals[s]>=bCutoff points start outside P
        Arrays.fill(oldSerVals,Double.POSITIVE_INFINITY);

        //for updating
        boolean firstFit = true;//first fit is not an update
        DoubleMatrix1D c = DoubleFactory1D.dense.make(numParams);//matrices we update (used in fit)
        DoubleMatrix2D M = DoubleFactory2D.dense.make(numParams,numParams);
        double oldFitWeights[] = null;
        
        
        //double fitWeightsCheck[] = fitWeights.clone();//DEBUG!!!

        while(!done){

            if(firstFit){
                coeffs = fitSeries(fitSamp, fitTrueVals, fitWeights, lambda,
                    includeConst, order, PCOrder, isPC, false, c, M);

                firstFit = false;
            }
            else{
                double weightDiffs[] = fitWeights.clone();
                for(int s=0; s<numRestraints; s++)
                    weightDiffs[s] -= oldFitWeights[s];

                coeffs = fitSeries(fitSamp, fitTrueVals, weightDiffs, lambda,
                    includeConst, order, PCOrder, isPC, true, c, M);
                
                
                
                //DEBUG!!!
                /*
                for(int s=0; s<numRestraints; s++){
                    fitWeightsCheck[s] += weightDiffs[s];
                    if(fitWeightsCheck[s] != fitWeights[s]){
                        int cefAO = 1111;
                    }
                }
                
                double checkCoeffs[] = fitSeries(fitSamp, fitTrueVals, fitWeights, lambda,
                    includeConst, order, PCOrder, isPC, false, null, null);
                
                for(int a=0; a<coeffs.length; a++){
                    if(Math.abs(checkCoeffs[a]-coeffs[a])>1e-10){
                        int abc=123;
                    }
                }
                
                DoubleMatrix1D c2 = DoubleFactory1D.dense.make(numParams);
                DoubleMatrix2D M2 = DoubleFactory2D.dense.make(numParams,numParams);
                
                double checkCoeffs2[] = fitSeries(fitSamp, fitTrueVals, fitWeights, lambda,
                    includeConst, order, PCOrder, isPC, false, c2, M2);
                
                for(int a=0; a<coeffs.length; a++){
                    if(Math.abs(checkCoeffs2[a]-coeffs[a])>1e-10){
                        int abc=123;
                    }
                }*/
                //DEBUG!!!
                
            }
            
            oldFitWeights = fitWeights.clone();


            done = true;
            ArrayList<SampleCutoffCrossing> scc = new ArrayList<SampleCutoffCrossing>();
            
            meanResidual = 0;
            weightSum = 0;


            
            //boolean doneNoTol = true, done2 = true;//DEBUG!!!
            
            //values of series at each sample, based on coeffs
            double serVals[] = new double[numSamples];
            for(int s=0; s<numSamples; s++){
                serVals[s] = evalSeries(coeffs,samp[s],nd,includeConst,order,PCOrder,isPC);

                //If each series value is below or above bcutoff according to whether the
                //coeffs were generated by fitting with or without that value included
                //(respectively), then we have found a local and thus the global minimum
                //in this quadratic piece of the objective function, so we're done


                //check for doneness first, using tolerance
                //i.e. are coeffs (derived from fitWeights) consistent with fitWeights,
                //within numerical error?  If so we have a global minimum

                
                if( trueVals[s]>=bCutoffs[s] ){
                    if(fitWeights[s]>0){
                        if(!isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],false,-1e-6))
                            done = false;//fitWeights penalizing deviation from bCutoff, and this isn't right at coeffs
                    }
                    else {
                        boolean secondRestraintOn = revSecondEntries.containsKey(s);
                        if(secondRestraintOn)
                            secondRestraintOn = (fitWeights[numSamples+revSecondEntries.get(s)]>0);
                        
                        if(secondRestraintOn){
                            if(!isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],true,-1e-6))
                                done = false;//fitWeights penalizing deviation from trueVal, and this isn't right at coeffs
                        }
                        else {//restraints currently off
                            if(isRestraintActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],1e-6))
                                done = false;//a restraint should be on at coeffs
                        }
                    }
                }
                //for trueVals below bCutoff, restraints don't turn on and off
                
                
                
                //DEBUG!!!!!
                //trying to calculate done w/o tol
                /*
                if( trueVals[s]>=bCutoffs[s] ){
                    if(fitWeights[s]>0){
                        if(!isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],false,0))
                            doneNoTol = false;//fitWeights penalizing deviation from bCutoff, and this isn't right at coeffs
                    }
                    else {
                        boolean secondRestraintOn = revSecondEntries.containsKey(s);
                        if(secondRestraintOn)
                            secondRestraintOn = (fitWeights[numSamples+revSecondEntries.get(s)]>0);
                        
                        if(secondRestraintOn){
                            if(!isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],true,0))
                                doneNoTol = false;//fitWeights penalizing deviation from trueVal, and this isn't right at coeffs
                        }
                        else {//restraints currently off
                            if(isRestraintActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],0))
                                doneNoTol = false;//a restraint should be on at coeffs
                        }
                    }
                }
                
                //OK now done2 will be calculated and should be the same but will be calculated like 
                //the weight changes below
                if( trueVals[s]>=bCutoffs[s] ){
                    if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],false)){
                        //activate penalty for deviating from bCutoff
                        if(fitWeights[s]!=weights[s])
                            done2 = false;
                        if(revSecondEntries.containsKey(s)){
                            if(fitWeights[numSamples+revSecondEntries.get(s)]!=0)
                                done2 = false;
                        }
                    }
                    else if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],true)){
                        //activate penalty for deviating from trueVal
                        if(fitWeights[s] != 0)
                            done2 = false;
                        if(revSecondEntries.containsKey(s)){
                            if(fitWeights[numSamples+revSecondEntries.get(s)] != weights[s])
                                done2 = false;
                        }
                        else
                            throw new RuntimeException("ERROR: should have second entry for restraint but don't!!");
                    }
                    else {
                        //deactivate all penalties
                        if(fitWeights[s] != 0)
                            done2 = false;
                        if(revSecondEntries.containsKey(s)){
                            if(fitWeights[numSamples+revSecondEntries.get(s)] != 0)
                                done2 = false;
                        }
                        //no contribution to residual
                    }
                }*/
                
                //DEBUG!!!
                
                
                
                
                
                //Now calculate mean residual and crossing points, and update fitWeights
                double residTerm = 0;
                
                                
                if( trueVals[s]>=bCutoffs[s] ){
                    if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],false)){
                        //activate penalty for deviating from bCutoff
                        fitWeights[s] = weights[s];
                        if(revSecondEntries.containsKey(s))
                            fitWeights[numSamples+revSecondEntries.get(s)] = 0;
                        residTerm = (serVals[s]-bCutoffs[s])*(serVals[s]-bCutoffs[s]);
                    }
                    else if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],true)){
                        //activate penalty for deviating from trueVal
                        fitWeights[s] = 0;
                        if(revSecondEntries.containsKey(s))
                            fitWeights[numSamples+revSecondEntries.get(s)] = weights[s];
                        else
                            throw new RuntimeException("ERROR: should have second entry for restraint but don't!!");
                        residTerm = (serVals[s]-trueVals[s])*(serVals[s]-trueVals[s]);
                    }
                    else {
                        //deactivate all penalties
                        fitWeights[s] = 0;
                        if(revSecondEntries.containsKey(s))
                            fitWeights[numSamples+revSecondEntries.get(s)] = 0;
                        //no contribution to residual
                    }
                }
                else //normal least-squares penalty.  fitWeights[s] will stay at weights[s] 
                    residTerm = (serVals[s]-trueVals[s])*(serVals[s]-trueVals[s]);
               
                
                meanResidual += weights[s]*residTerm;
                

                //If want sample-by-sample output...
                //System.out.println("TRAININGSET TRUE: "+trueVals[s]+" SER: "+serVals[s]);

                weightSum += weights[s];
            }

            meanResidual /= weightSum;

            if(meanResidual==prevResid)
                System.out.println();



            
            
            
            
            
            
            //DEBUG!!!
            /*
            if(done!=doneNoTol || done2!=done){
                //Let's see what happens if we remove the tolerance...
                done = doneNoTol;
            }
            
            
            
            if(done){
                double checkCoeffs[] = fitSeries(fitSamp, fitTrueVals, fitWeights, lambda,
                    includeConst, order, PCOrder, isPC, false, null, null);
                
                for(int a=0; a<coeffs.length; a++){
                    if(Math.abs(checkCoeffs[a]-coeffs[a])>1e-10){
                        int abc=123;
                    }
                }
                
            }*/
            //DEBUG!!!
            
            
            


            if( (!done) && (meanResidual>=prevResid) ) {
                //Did not obtain a decrease using the Newton step
                //Let's do an exact line search to rectify the situation

                if(!useLineSearch){
                    System.out.println("Skipping line search, returning with residual "+prevResid);
                    return oldCoeffs;
                }
                
                
                System.out.println("LINE SEARCH");
                
                
                for(int s=0; s<numSamples; s++){
                    
                    //If we go in or out of either type of restraint between serVals and oldSerVals, we create
                    //a SampleCutoffCrossing of the appropriate type (upper or lower (ordinary) restraint)
                    if( (isRestraintTypeActive(trueVals[s],oldSerVals[s],bCutoffs[s],bCutoffs2[s],false)) !=
                            (isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],false)) ){
                        //If the restraint disappears at one end we know trueVal>=bCutoff here
                        //create lower restraint SampleCutoffCrossing
                        double crossingPoint = (bCutoffs[s]-oldSerVals[s]) / (serVals[s]-oldSerVals[s]);
                        scc.add( new SampleCutoffCrossing(s,crossingPoint,false) );
                    }
                    
                    if( (isRestraintTypeActive(trueVals[s],oldSerVals[s],bCutoffs[s],bCutoffs2[s],true)) !=
                            (isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],true)) ){
                        
                        //upper
                        double crossingPoint2 = (trueVals[s]-oldSerVals[s]) / (serVals[s]-oldSerVals[s]);
                        scc.add( new SampleCutoffCrossing(s,crossingPoint2,true) );
                    }
                    
                }
                

                int changeCount = scc.size();
                Collections.sort(scc);

                TreeSet<Integer> crossingIndices = new TreeSet<Integer>();
                for(SampleCutoffCrossing cr : scc){
                    int s = cr.sampleIndex;
                    crossingIndices.add(s);
                    if(cr.upperRestraint){//penalizing difference from trueVal
                        cr.quadTerm = weights[s]*(serVals[s]-oldSerVals[s])*(serVals[s]-oldSerVals[s]);
                        cr.linTerm = 2*weights[s]*(serVals[s]-oldSerVals[s])*(oldSerVals[s]-trueVals[s]);
                        cr.constTerm =  weights[s]*(oldSerVals[s]-trueVals[s])*(oldSerVals[s]-trueVals[s]);
                    }
                    else{//penalizing difference from lesser of bCutoff, trueVal
                        cr.quadTerm = weights[s]*(serVals[s]-oldSerVals[s])*(serVals[s]-oldSerVals[s]);
                        double baseVal = Math.min(trueVals[s], bCutoffs[s]);
                        cr.linTerm = 2*weights[s]*(serVals[s]-oldSerVals[s])*(oldSerVals[s]-baseVal);
                        cr.constTerm =  weights[s]*(oldSerVals[s]-baseVal)*(oldSerVals[s]-baseVal);
                        //cr.linTerm = 2*weights[s]*(serVals[s]-oldSerVals[s])*(oldSerVals[s]-trueVals[s]);
                        //cr.constTerm =  weights[s]*(oldSerVals[s]-trueVals[s])*(oldSerVals[s]-trueVals[s]);
                    }
                }

                //Set up quadratic function
                double quadTerm = 0;
                double linTerm = 0;
                double constTerm = 0;
                //Add in contributions from all non-cutoff-crossing points
                for(int s=0; s<numSamples; s++){
                    if(!crossingIndices.contains(s)){
                        if(isRestraintTypeActive(trueVals[s],oldSerVals[s],bCutoffs[s],bCutoffs2[s],false)){
                        //if(trueVals[s]<bCutoffs[s]||serVals[s]<bCutoffs[s]){//penalizing difference from lesser of bCutoff, trueVal
                            quadTerm += weights[s]*(serVals[s]-oldSerVals[s])*(serVals[s]-oldSerVals[s]);
                            double baseVal = Math.min(trueVals[s], bCutoffs[s]);
                            linTerm += 2*weights[s]*(serVals[s]-oldSerVals[s])*(oldSerVals[s]-baseVal);
                            constTerm +=  weights[s]*(oldSerVals[s]-baseVal)*(oldSerVals[s]-baseVal);
                            //linTerm += 2*weights[s]*(serVals[s]-oldSerVals[s])*(oldSerVals[s]-trueVals[s]);
                            //constTerm +=  weights[s]*(oldSerVals[s]-trueVals[s])*(oldSerVals[s]-trueVals[s]);
                        }
                        else if(isRestraintTypeActive(trueVals[s],oldSerVals[s],bCutoffs[s],bCutoffs2[s],true)){
                        //else if(isRestraintActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],0)){//penalizing difference from trueVal
                            quadTerm += weights[s]*(serVals[s]-oldSerVals[s])*(serVals[s]-oldSerVals[s]);
                            linTerm += 2*weights[s]*(serVals[s]-oldSerVals[s])*(oldSerVals[s]-trueVals[s]);
                            constTerm +=  weights[s]*(oldSerVals[s]-trueVals[s])*(oldSerVals[s]-trueVals[s]);                            
                        }
                    }
                }
                
                //contributions from cutoff-crossing points at the beginning of the interval
                //(i.e. at coeffs)
                for(SampleCutoffCrossing cr : scc){
                    int s = cr.sampleIndex;
                    if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],cr.upperRestraint)){
                        //check if this particular restraint (rather than either restraint for this s) is active
                    //if(serVals[s]<bCutoffs[s]){
                        quadTerm += cr.quadTerm;
                        linTerm += cr.linTerm;
                        constTerm += cr.constTerm;
                    }
                }


                //double checkMeanResid = (quadTerm+linTerm+constTerm)/weightSum;//evaluate objective function at a=1
                //should match meanResidual!


                double prevNodeResid = Double.POSITIVE_INFINITY;
                //The first increase we may consider is from node 0 to 1

                //Now walk back until we get an increase
                //then the minimum will be in one of the last two quadratic pieces
                int lowestNodeIndex = 0;
                for(int curChangeIndex=changeCount-1; curChangeIndex>=0; curChangeIndex--){
                    
                    SampleCutoffCrossing cr = scc.get(curChangeIndex);

                    double a = cr.crossingPoint;
                    double curNodeResid = (a*a*quadTerm + a*linTerm + constTerm)/weightSum;
                    

                    int s = cr.sampleIndex;


                    if(curNodeResid>prevNodeResid){
                        lowestNodeIndex = curChangeIndex + 1;
                        break;
                    }

                    //if the restraint is being removed going from serVals to oldSerVals, remove its coefficients from the quadratic function
                    //if it's being added, add them
                    if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],cr.upperRestraint)){
                        quadTerm -= cr.quadTerm;
                        linTerm -= cr.linTerm;
                        constTerm -= cr.constTerm;
                    }
                    else  {
                        quadTerm += cr.quadTerm;
                        linTerm += cr.linTerm;
                        constTerm += cr.constTerm;
                    }

                    prevNodeResid = curNodeResid;
                }


                //At this point, we know our minimum is in either the piece with the
                //current quad, lin, constTerms or the one we looked at right before
                //(where lowestNodeIndex is the node separating them)
                double a_min = -linTerm/(2*quadTerm);
                SampleCutoffCrossing cr = scc.get(lowestNodeIndex);
                if(a_min>cr.crossingPoint){//minimum must be in previous piece
                    //revert quad and linTerms
                    int s = cr.sampleIndex;
                    //change back quadratic-function coefficients
                    if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],cr.upperRestraint)){
                        quadTerm += cr.quadTerm;
                        linTerm += cr.linTerm;
                        constTerm += cr.constTerm;
                    }
                    else  {
                        quadTerm -= cr.quadTerm;
                        linTerm -= cr.linTerm;
                        constTerm -= cr.constTerm;
                    }
                    a_min = -linTerm/(2*quadTerm);
                }

                //double minResid = (a_min*a_min*quadTerm + a_min*linTerm + constTerm)/weightSum;


                for(int p=0; p<numParams; p++)
                    coeffs[p] = coeffs[p]*a_min + (1-a_min)*oldCoeffs[p];


                double minResid = 0;
                for(int s=0; s<numSamples; s++){
                    serVals[s] = evalSeries(coeffs,samp[s],nd,includeConst,order,PCOrder,isPC);

                    if( trueVals[s]>=bCutoffs[s] ){
                        if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],false)){
                            //activate penalty for deviating from bCutoff
                            fitWeights[s] = weights[s];
                            if(revSecondEntries.containsKey(s))
                                fitWeights[numSamples+revSecondEntries.get(s)] = 0;
                            minResid += weights[s]*(serVals[s]-bCutoffs[s])*(serVals[s]-bCutoffs[s]);
                        }
                        else if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],true)){
                            minResid += weights[s]*(serVals[s]-trueVals[s])*(serVals[s]-trueVals[s]);
                            //activate penalty for deviating from trueVal
                            fitWeights[s] = 0;
                            fitWeights[numSamples+revSecondEntries.get(s)] = weights[s];
                        }
                        else {
                            //deactivate all penalties
                            fitWeights[s] = 0;
                            if(revSecondEntries.containsKey(s))
                                fitWeights[numSamples+revSecondEntries.get(s)] = 0;
                        }
                    }
                    else
                        minResid += weights[s]*(serVals[s]-trueVals[s])*(serVals[s]-trueVals[s]);
                }

                minResid /= weightSum;

                if(minResid>=prevResid){
                    //consider to have converged
                    //NOTE THIS CAN HAPPEN IF THE QUADRATIC APPROXIMATION AT OLDCOEFFS HAS SOLUTION
                    //FAR FROM THE EXACT VALUE (ASSUMING EXACT FITSERIES) OF 1
                    //THIS CAN HAPPEN IF WE'RE GETTING BELOW THE NUMERICAL PRECISION OF FITSERIES
                    System.out.println("TRAINING SET MEAN RESIDUAL:"+prevResid);
                    System.out.println("CONVERGED IN LINE SEARCH, line search min: "+minResid);
                    System.out.println("fitSeriesIterative time (ms): "+(System.currentTimeMillis()-startTime));

                    return oldCoeffs;
                }

                meanResidual = minResid;
            }
                
                
/*                
                //trying backtracking line search, hoping it'll be more numerically stable
                double backtrackRate = 0.5;
                double a = 1;//how far we want to be on the line from oldCoeffs (which should allow some descent
                //along the line) to coeffs (which overshoots)
                double minResid;
                double aCoeffs[] = new double[coeffs.length];//coefficients at our point on the line
                
                do {
                    
                    a *= backtrackRate;
                    
                    for(int p=0; p<numParams; p++)
                        aCoeffs[p] = coeffs[p]*a + (1-a)*oldCoeffs[p];


                    minResid = 0;
                    for(int s=0; s<numSamples; s++){
                        serVals[s] = evalSeries(aCoeffs,samp[s],nd,includeConst,order,PCOrder,isPC);

                        if( trueVals[s]>=bCutoffs[s] ){
                            if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],false)){
                                minResid += weights[s]*(serVals[s]-bCutoffs[s])*(serVals[s]-bCutoffs[s]);
                                weights2[s] = weights[s];
                            }
                            else if(isRestraintTypeActive(trueVals[s],serVals[s],bCutoffs[s],bCutoffs2[s],true)){
                                minResid += weights[s]*(serVals[s]-trueVals[s])*(serVals[s]-trueVals[s]);
                                weights2[s] = weights[s];
                            }
                            else
                                weights2[s] = 0;
                        }
                        else
                            minResid += weights[s]*(serVals[s]-trueVals[s])*(serVals[s]-trueVals[s]);

                    }

                    minResid /= weightSum;
                    
                } while(minResid>prevResid);
                
                if(a<1e-4)
                    System.out.println("Warning: line search a got down to "+a);

                
                meanResidual = minResid;
            }
*/

            oldCoeffs = coeffs;
            System.out.println("STEP RESIDUAL: "+meanResidual);
            prevResid = meanResidual;
            oldSerVals = serVals;
        }

        System.out.println("TRAINING SET MEAN RESIDUAL:"+meanResidual);
        System.out.println("fitSeriesIterative time (ms): "+(System.currentTimeMillis()-startTime));

        
        
        
        
        
        //DEBUG!!!
        //Note this assumes weights are all 1!
        /*
        int numCoeffs = coeffs.length;
        DoubleMatrix1D residGrad = DoubleFactory1D.dense.make(numCoeffs);
        double checkResid = 0;
        for(int s=0; s<samp.length; s++){
            double diff = 0;
            double fitVal = evalSeries(coeffs,samp[s],nd,includeConst,order,PCOrder,isPC);
            if(trueVals[s]>=bCutoffs[s] 
                    && fitVal<bCutoffs[s]){
                diff = fitVal-bCutoffs[s];
                
                if(fitWeights[s]==0){
                    int allosaurus = 15;
                }
                if(revSecondEntries.containsKey(s)){
                    if(fitWeights[numSamples+revSecondEntries.get(s)]>0){
                        int allosaurus=15;
                    }
                }
                    
            }
            else if( trueVals[s]<bCutoffs[s] 
                    || (trueVals[s]<bCutoffs2[s] && fitVal>trueVals[s]) ){
             
                diff = fitVal-trueVals[s];
                
                if(trueVals[s]<bCutoffs[s]){
                    if(fitWeights[s]==0){
                        int allosaurus = 15;
                    }
                }
                else if(fitWeights[numSamples+revSecondEntries.get(s)]==0 || fitWeights[s]>0){
                    int allosaurus=15;
                }
            }
            else if(fitWeights[s]>0){
                int allosaurus=15;
            }
            else if(revSecondEntries.containsKey(s)){
                if(fitWeights[numSamples+revSecondEntries.get(s)]>0){
                    int allosaurus=15;
                }
            }
            
            checkResid += diff*diff;
            
            DoubleMatrix1D termMonomials = DoubleFactory1D.dense.make(numCoeffs);
            SeriesFitter.calcSampParamCoeffs( termMonomials, samp[s],
                    nd, includeConst, order, PCOrder, isPC );
            residGrad.assign( termMonomials, Functions.plusMult(2*diff) );
        }
        
        checkResid /= samp.length;
        residGrad.assign( Functions.mult(1./samp.length) );
        
        
        
        double checkCoeffs[] = fitSeries(fitSamp, fitTrueVals, fitWeights, lambda,
                    includeConst, order, PCOrder, isPC, false, null, null);
        
        double checkResid2 = 0;
        DoubleMatrix1D residGrad2 = DoubleFactory1D.dense.make(numCoeffs);

        for(int s=0; s<samp.length; s++){
            
            double diff = 0;
            double fitVal = evalSeries(checkCoeffs,samp[s],nd,includeConst,order,PCOrder,isPC);
            if(trueVals[s]>=bCutoffs[s] 
                    && fitVal<bCutoffs[s]){
                diff = fitVal-bCutoffs[s];
                
                if(fitWeights[s]==0){
                    int allosaurus = 15;
                }
                if(revSecondEntries.containsKey(s)){
                    if(fitWeights[numSamples+revSecondEntries.get(s)]>0){
                        int allosaurus=15;
                    }
                }
                    
            }
            else if( trueVals[s]<bCutoffs[s] 
                    || (trueVals[s]<bCutoffs2[s] && fitVal>trueVals[s]) ){
             
                diff = fitVal-trueVals[s];
                
                if(trueVals[s]<bCutoffs[s]){
                    if(fitWeights[s]==0){
                        int allosaurus = 15;
                    }
                }
                else if(fitWeights[numSamples+revSecondEntries.get(s)]==0 || fitWeights[s]>0){
                    int allosaurus=15;
                }
            }
            else if(fitWeights[s]>0){
                int allosaurus=15;
            }
            else if(revSecondEntries.containsKey(s)){
                if(fitWeights[numSamples+revSecondEntries.get(s)]>0){
                    int allosaurus=15;
                }
            }
            
            checkResid2 += diff*diff;
            
            DoubleMatrix1D termMonomials = DoubleFactory1D.dense.make(numCoeffs);
            SeriesFitter.calcSampParamCoeffs( termMonomials, samp[s],
                    nd, includeConst, order, PCOrder, isPC );
            residGrad2.assign( termMonomials, Functions.plusMult(2*diff) );
        }
        
        checkResid2 /= samp.length;
        residGrad2.assign( Functions.mult(1./samp.length) );
        
        
        
        
        if(residGrad.zDotProduct(residGrad)>1e-10){
            int struthiomimus = 23;
        }
        int cubist = -1;
        */
        //DEBUG!!
        
        
        
        
        
        
        
        
        return coeffs;
    }
    
    
    static boolean isRestraintActive(double trueVal, double serVal, double bCutoff, double bCutoff2, double tol){
        //Figure out if the restraint for the given sample is active at serVal
        //If tol>0, then err on the side of being inactive; if <0, active; if ==0, neutral
        if(trueVal<bCutoff || serVal<bCutoff-tol)//if lower or 2-sided constraint is active
            return true;
        if(trueVal>=bCutoff2 || serVal<trueVal+tol)//if upper constraint is inactive
            return false;
        return true;
    }
    
    
    static boolean isRestraintTypeActive(double trueVal, double serVal, double bCutoff, double bCutoff2, boolean type, double tol){
        //If type is false, check if the ordinary restraint (2-sided if below trueVal<bCutoff; lower restraint centered at bCutoff
        //otherwise) is active
        //If type is true, check if the upper restraint (penalizing serVal>trueVal for bCutoff<=trueVal<bCutoff2) is active
        if(type)
            return (trueVal>=bCutoff) && (trueVal<bCutoff2) && (serVal>trueVal+tol);
        else
            return (trueVal<bCutoff || serVal<bCutoff-tol);
    }
    
    static boolean isRestraintTypeActive(double trueVal, double serVal, double bCutoff, double bCutoff2, boolean type){
        //default tolerance is 0
        return isRestraintTypeActive(trueVal,serVal,bCutoff,bCutoff2,type,0);
    }
    
    

    private static class SampleCutoffCrossing implements Comparable {
        //These are samples whose trueVal is above bCutoff
        //and whose series approximation crosses bCutoff
        //in the line segment between the previous and current Newton steps'
        //coefficients
        
        int sampleIndex;
        double crossingPoint;//where on the line segment the crossing occurs
        //(0=previous, 1=current coefficients, linearly interpolated)

        double quadTerm, linTerm, constTerm;//contributions of sample
        //to terms of quadratic specifying objective function (in pieces)
        //along line search line--not in constructor, we'll fill them in if needed for
        //line search
        
        boolean upperRestraint = false;//upperRestraint indicates
        //that we're restraining the series to stay below a true value between the 2 bCutoffs

        public SampleCutoffCrossing(int sampleIndex, double crossingPoint, boolean upperRestraint) {
            this.sampleIndex = sampleIndex;
            this.crossingPoint = crossingPoint;
            this.upperRestraint = upperRestraint;
        }

        
        public int compareTo(Object o) {//We'll want to sort by crossingPoint
            return Double.valueOf(crossingPoint).compareTo( ((SampleCutoffCrossing)o).crossingPoint );
        }
    }




    //default is w/o PCs
    public static double evalSeries(double[] coeffs, DoubleMatrix1D x, int nd, boolean includeConst,
            int order){

        return evalSeries(coeffs,x,nd,includeConst,order,order,null);
    }


    
    static double evalSeries(double[] coeffs, DoubleMatrix1D x, int nd, boolean includeConst, 
            int order, int PCOrder, boolean isPC[]){

        
        //we support orders 1-6
        //and PCOrder may range up to 6 (but has no effect if <=order)
        if(order<1||order>6||PCOrder>6){
            throw new RuntimeException("ERROR: SeriesFitter.evalSeries does not support order "+order+" and/or PCOrder "+PCOrder);
        }
        
        
        if(order==1 && PCOrder==2)
            throw new RuntimeException("ERROR: Order 1 and PCOrder 2 not supported");
        
            //evaluate series at point x from flat-array coefficients
            //assuming up to quadratic is included with all DOFs
        int count = 0;
        double ans = 0;

        if(includeConst){
            ans += coeffs[0];
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            ans += coeffs[count]*x.get(dof);
            count++;
        }

        
        if(order>=2){
            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<dof; dof2++){
                    ans += coeffs[count]*x.get(dof)*x.get(dof2);
                    count++;
                }
                ans += coeffs[count]*x.get(dof)*x.get(dof);
                count++;
            }
        }


        if(order>=3){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){

                        ans += coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3);
                        count++;
                    }
                }
            }
        }
        else if(PCOrder>=3){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    ans += coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3);
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=4){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            ans += coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4);
                            count++;
                        }
                    }
                }
            }
        }
        else if(PCOrder>=4){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            ans += coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4);
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        if(order>=5){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    ans += coeffs[count]*x.get(dof)*x.get(dof2)*
                                                            x.get(dof3)*x.get(dof4)*x.get(dof5);
                                                    count++;
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=5){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    ans += coeffs[count]*x.get(dof)*x.get(dof2)*
                                                            x.get(dof3)*x.get(dof4)*x.get(dof5);
                                                    count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=6){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                            ans += coeffs[count]*x.get(dof)*x.get(dof2)*
                                                                    x.get(dof3)*x.get(dof4)*x.get(dof5)*
                                                                    x.get(dof6);
                                                            count++;
                                                    }
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=6){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        if(isPC[dof6]){

                                                            ans += coeffs[count]*x.get(dof)*x.get(dof2)*
                                                                    x.get(dof3)*x.get(dof4)*x.get(dof5)*
                                                                    x.get(dof6);
                                                            count++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return ans;
    }




    static double[] evalSeriesByDegree(double[] coeffs, DoubleMatrix1D x, int nd, boolean includeConst, int order) {
        //like evalSeries, but return terms of each degree in the polynomial separately
        //(so ans[i] is the i-th degree monomials)
        //this only goes up to 4

        if(order>4){
            throw new Error("SeriesFitter.evalSeriesByDegree doesn't support order "+order);
        }
        
        int count = 0;//for counting through coeffs

        int topDegree = order;

        double ans[] = new double[topDegree+1];

        if(includeConst){
            ans[0] = coeffs[0];
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            ans[1] += coeffs[count]*x.get(dof);
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            for(int dof2=0; dof2<dof; dof2++){
                ans[2] += coeffs[count]*x.get(dof)*x.get(dof2);
                count++;
            }
            ans[2] += coeffs[count]*x.get(dof)*x.get(dof);///2;
            count++;
        }


        if(order>=3){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){

                        ans[3] += coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3);
                                //*distingPerm(dof,dof2,dof3)/6;
                        count++;
                    }
                }
            }
        }

        if(order>=4){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            ans[4] += coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4);
                                    //*distingPerm(dof,dof2,dof3,dof4)/24;
                            count++;
                        }
                    }
                }
            }
        }

        return ans;
    }






    /*private static int distingPerm(int... a){
        //given a set of integers in nondescending order
        //count how many distinguishable permutations of them there are
        int n=a.length;
        int ans = factorial(n);
        int curNum = a[0], count=1;
        for(int b=1; b<n; b++){
            if(a[b]==curNum)
                count++;
            else{
                ans /= factorial(count);
                count = 1;
                curNum = a[b];
            }
        }
        ans /= factorial(count);
        return ans;
    }
*/

    static int factorial(int a){
        if(a<=1)
            return 1;
        return a*factorial(a-1);
    }

    
    //A polynomial fit is a linear least-squares fit with respect to the parameters
    //This method calculates the contribution of a given sample to the coefficients in this linear fit
    //for each parameter
    //nd = number of dimensions
    static void calcSampParamCoeffs(DoubleMatrix1D v, DoubleMatrix1D sample, int nd, 
            boolean includeConst, int order, int PCOrder, boolean isPC[]){

        if(order<1||order>6||PCOrder>6){
            throw new Error("SeriesFitter.calcSampParamCoeffs does not support order "
                    +order+" and/or PCOrder "+PCOrder);
        }
        
        if(order==1 && PCOrder==2)
            throw new RuntimeException("ERROR: Order 1 and PCOrder 2 not supported");
        
        int count=0;

        if(includeConst){
            v.set(count, 1);
            count++;
        }

        for(int i=0; i<nd; i++){
            v.set(count, sample.get(i));
            count++;
        }


        if(order>=2){
            for(int a=0; a<nd; a++){
                for(int b=0; b<=a; b++){
                    double term = sample.get(a)*sample.get(b);
                    v.set(count, term);
                    count++;
                }
            }
        }

        
        if(order>=3){
            
            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){

                        v.set( count, sample.get(dof)*sample.get(dof2)*sample.get(dof3) );
                        count++;
                    }
                }
            }
        }
        else if(PCOrder>=3){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    v.set( count, sample.get(dof)*sample.get(dof2)*sample.get(dof3) );
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        }

        
        if(order>=4){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            
                            v.set( count, sample.get(dof)*sample.get(dof2)*sample.get(dof3)
                                    *sample.get(dof4) );
                            count++;
                        }
                    }
                }
            }
        }
        else if(PCOrder>=4){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            v.set( count, sample.get(dof)*sample.get(dof2)*sample.get(dof3)
                                                *sample.get(dof4) );
                                            
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        
        if(order>=5){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    v.set( count, sample.get(dof)*sample.get(dof2)*
                                                            sample.get(dof3)*sample.get(dof4)*sample.get(dof5) );
                                                    count++;
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=5){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    v.set( count, sample.get(dof)*sample.get(dof2)*
                                                            sample.get(dof3)*sample.get(dof4)*sample.get(dof5) );
                                                    count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        
        
        if(order>=6){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        
                                                            v.set( count, sample.get(dof)*sample.get(dof2)*
                                                                sample.get(dof3)*sample.get(dof4)*
                                                                sample.get(dof5)*sample.get(dof6) );
                                                            
                                                            count++;
                                                    }
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=6){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        if(isPC[dof6]){

                                                            v.set( count, sample.get(dof)*sample.get(dof2)*
                                                                sample.get(dof3)*sample.get(dof4)*
                                                                sample.get(dof5)*sample.get(dof6) );
                                                            
                                                            count++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


    }








    public static int getNumParamsForOrder(int nd, int order){
        //Get the number of coefficients for terms of the specified order
        //in a series expansion with nd variables
        if(order==0)
            return 1;
        else if(order==1)
            return nd;
        else if(order==2)
            return nd*(nd+1)/2;
        else if(order==3)
            return nd + nd*(nd-1) + nd*(nd-1)*(nd-2)/6;
        else if(order==4)
            return nd + 3*nd*(nd-1)/2 + nd*(nd-1)*(nd-2)/2 + nd*(nd-1)*(nd-2)*(nd-3)/24;
        else if(order==5 || order==6){

            //WOULD BE NICE TO DO THE COMBINATORICS SOMETIME
            int count = 0;
            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            for(int dof5=0; dof5<=dof4; dof5++){

                                if(order==6){
                                    for(int dof6=0; dof6<=dof5; dof6++)
                                            count++;
                                }
                                else
                                    count++;
                            }
                        }
                    }
                }
            }

            return count;
        }

        throw new Error("ORDER NOT SUPPORTED IN SeriesFitter.getNumParamsForOrder: "+order);
    }
    
    
    
    public static int getNumParams(int nd, boolean includeConst, int order){
        
        int ans = 0;
        if(includeConst)
            ans++;
        
        for(int ord=1; ord<=order; ord++)
            ans += getNumParamsForOrder(nd,ord);
        
        return ans;
    }

    
    static int countTrue(boolean[] a){
        int count = 0;

        for(boolean b : a)
            if(b)
                count++;

        return count;
    }
    
    
    
    static DoubleMatrix2D getHessian(double coeffs[], int numDOFs, boolean includeConst){
        //Extract the Hessian from the given series coeffs, 
        //which are expected to have order>=2
        
        int count = numDOFs;//start at first index in coeffs for a second-order coefficient
        if(includeConst)
            count++;
        
        DoubleMatrix2D hess = DoubleFactory2D.dense.make(numDOFs,numDOFs);

        for(int a=0; a<numDOFs; a++){
            for(int b=0; b<=a; b++){
                double hessVal = coeffs[count];
                if(b==a)
                    hessVal *= 2;
                hess.set(a, b, hessVal);
                hess.set(b, a, hessVal);
                count++;
            }
        }
        
        return hess;
    }
    
    
    
    static DoubleMatrix1D evalSeriesGradient(double[] coeffs, DoubleMatrix1D x, int nd, boolean includeConst, 
            int order, int PCOrder, boolean isPC[]){

        //like evalSeries but evaluating the gradient
        //this function is structured like evalSeries
        //counting each term's contribution to the gradient
        
        
        
        int count = 0;
        if(includeConst)
            count++;
        
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(nd);

        for(int dof=0; dof<nd; dof++){
            ans.set( dof, coeffs[count] );
            count++;
        }

        for(int dof=0; dof<nd; dof++){
            for(int dof2=0; dof2<dof; dof2++){
                ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2) );
                ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof) );
                count++;
            }
            ans.set(dof, ans.get(dof)+2*coeffs[count]*x.get(dof) );
            count++;
        }


        if(order>=3){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){

                        ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3) );
                        ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3) );
                        ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2) );

                        count++;
                    }
                }
            }
        }
        else if(PCOrder>=3){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3) );
                                    ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3) );
                                    ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2) );

                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=4){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            
                            ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3)*x.get(dof4) );
                            ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3)*x.get(dof4) );
                            ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof4) );
                            ans.set(dof4, ans.get(dof4)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3) );

                            count++;
                        }
                    }
                }
            }
        }
        else if(PCOrder>=4){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){
                                            
                                            ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3)*x.get(dof4) );
                                            ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3)*x.get(dof4) );
                                            ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof4) );
                                            ans.set(dof4, ans.get(dof4)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3) );

                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        //ONLY SUPPORTING QUINTIC, HEXIC IN PC FORM
        if(order>=5){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                                                                        
                                                    ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof5) );
                                                    ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3)*x.get(dof4)*x.get(dof5) );
                                                    ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof4)*x.get(dof5) );
                                                    ans.set(dof4, ans.get(dof4)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof5) );
                                                    ans.set(dof5, ans.get(dof5)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4) );

                                                    count++;
                                                }
                                        }
                                }
                        }
                }
        }
        else if(PCOrder>=5){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){
                                                                                                        
                                                    ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof5) );
                                                    ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3)*x.get(dof4)*x.get(dof5) );
                                                    ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof4)*x.get(dof5) );
                                                    ans.set(dof4, ans.get(dof4)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof5) );
                                                    ans.set(dof5, ans.get(dof5)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4) );

                                                    count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        if(order>=6){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                            
                                                            ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3)*x.get(dof4)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof4)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof4, ans.get(dof4)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof5, ans.get(dof5)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof6) );
                                                            ans.set(dof6, ans.get(dof6)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof5) );

                                                            count++;
                                                    }
                                            }
                                    }
                            }
                    }
            }
            
        }
        else if(PCOrder>=6){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        if(isPC[dof6]){
                                                            
                                                            
                                                            ans.set(dof, ans.get(dof)+coeffs[count]*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof2, ans.get(dof2)+coeffs[count]*x.get(dof)*x.get(dof3)*x.get(dof4)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof3, ans.get(dof3)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof4)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof4, ans.get(dof4)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof5)*x.get(dof6) );
                                                            ans.set(dof5, ans.get(dof5)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof6) );
                                                            ans.set(dof6, ans.get(dof6)+coeffs[count]*x.get(dof)*x.get(dof2)*x.get(dof3)*x.get(dof4)*x.get(dof5) );

                                                            count++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return ans;
    }
    
    
    
    static DoubleMatrix2D evalSeriesHessian(double[] coeffs, DoubleMatrix1D x, int nd, boolean includeConst, 
            int order, int PCOrder, boolean isPC[]){

        //like evalSeriesGradient but Hessian
        
        int count = nd;//we can skip linear terms
        if(includeConst)//skip one more if there's a constant
            count++;
        
        DoubleMatrix2D ans = DoubleFactory2D.dense.make(nd,nd);

        
        for(int dof=0; dof<nd; dof++){
            for(int dof2=0; dof2<dof; dof2++){
                ans.set( dof, dof2, coeffs[count] );
                ans.set( dof2, dof, coeffs[count] );
                count++;
            }
            ans.set(dof, dof, 2*coeffs[count]);
            count++;
        }


        if(order>=3){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        
                        updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3);
                        count++;
                    }
                }
            }
        }
        else if(PCOrder>=3){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3);
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        }



        if(order>=4){

            for(int dof=0; dof<nd; dof++){
                for(int dof2=0; dof2<=dof; dof2++){
                    for(int dof3=0; dof3<=dof2; dof3++){
                        for(int dof4=0; dof4<=dof3; dof4++){
                            
                            updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3,dof4);
                            count++;
                        }
                    }
                }
            }
        }
        else if(PCOrder>=4){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){
                                            
                                            updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3,dof4);
                                            count++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        
        if(order>=5){
            
            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                                                                        
                                                    updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3,dof4,dof5);
                                                    count++;
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=5){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){
                                                                                                        
                                                    updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3,dof4,dof5);
                                                    count++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        if(order>=6){

            for(int dof=0; dof<nd; dof++){
                    for(int dof2=0; dof2<=dof; dof2++){
                            for(int dof3=0; dof3<=dof2; dof3++){
                                    for(int dof4=0; dof4<=dof3; dof4++){
                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        
                                                            updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3,dof4,dof5,dof6);
                                                            count++;
                                                    }
                                            }
                                    }
                            }
                    }
            }
        }
        else if(PCOrder>=6){

            for(int dof=0; dof<nd; dof++){
                if(isPC[dof]){

                    for(int dof2=0; dof2<=dof; dof2++){
                        if(isPC[dof2]){

                            for(int dof3=0; dof3<=dof2; dof3++){
                                if(isPC[dof3]){

                                    for(int dof4=0; dof4<=dof3; dof4++){
                                        if(isPC[dof4]){

                                            for(int dof5=0; dof5<=dof4; dof5++){
                                                if(isPC[dof5]){

                                                    for(int dof6=0; dof6<=dof5; dof6++){
                                                        if(isPC[dof6]){
                                                            
                                                            updateSeriesHessian(ans,x,coeffs[count],dof,dof2,dof3,dof4,dof5,dof6);
                                                            count++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return ans;
    }
    
    
    static void updateSeriesHessian(DoubleMatrix2D hess, DoubleMatrix1D x, double coeff, int...dofs ){
        //update the Hessian hess of a series at DOF values x
        //by adding in the Hessian of the term
        //coeff*prod_i x.get(dofs[i])
        int deg = dofs.length;
        
        for(int d1=0; d1<deg; d1++){
            for(int d2=0; d2<d1; d2++){//non diagonal Hessian contributions
                double dhess = coeff;
                
                for(int d3=0; d3<deg; d3++){
                    if(d3!=d1 && d3!=d2)
                        dhess *= x.get(dofs[d3]);
                }
                
                hess.set(dofs[d1],dofs[d2], hess.get(dofs[d1],dofs[d2])+dhess);
                hess.set(dofs[d2],dofs[d1], hess.get(dofs[d2],dofs[d1])+dhess);
            }
        }
    }
    
    

}