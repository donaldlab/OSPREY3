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
//	SubThreshSampler.java
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
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;


public class SubThreshSampler {
    //This class uses Metropolis to sample points in a voxel whose energies,
    //given by the objective function, are below the given threshold
    //used to generate samples for EPIC fits when sampling uniformly yields too many points
    //above the threshold
    
    double thresh;
    ObjectiveFunction of;
    DoubleMatrix1D DOFmin, DOFmax;
    //we want our samples x to lie in the bounds given by DOFmin, DOFmax
    //and to have of.get(x)<=thresh
    int numDOFs;
            

    DoubleMatrix1D samplingScale;//each step we take will perturb DOF #dof
    //by samplingScale.get(dof)+standard normal variable
    
    
    //we'll tune the samplingScale using constant multipliers
    //we'll try numCandScaleTuning candidates for each scaling, and try to get the desired acceptRatioTarget
    static int numCandScaleTuning = 25;
    static double[] acceptRatioTarget = new double[] {0.2,0.4};
    
    int useFrequency = 1;//out of samples accepted, frequency of samples to actually use
    
    
    DoubleMatrix1D x;//current sample value
    Random random = new Random();//a random
    
    static boolean adaptiveScale = true;//if true, then we adaptively choose our samplingScale at each step
    //it will be a function of x, and samplingScale will be specific to x

    public SubThreshSampler(double thresh, ObjectiveFunction of, DoubleMatrix1D DOFmin, DoubleMatrix1D DOFmax) {
        this.thresh = thresh;
        this.of = of;
        this.DOFmin = DOFmin;
        this.DOFmax = DOFmax;
        numDOFs = DOFmin.size();
    }
    
    
    
    void initScale(){
        //get an initial guess for samplingScale
        samplingScale = getScaleAdaptive(x);
    }
        
    DoubleMatrix1D getScaleAdaptive(DoubleMatrix1D sp){
        //Get an appropriate scale for sampling around some starting point sp
        /*
        //initially, we'll go about a tenth of the way across the voxel
        //in each dimension
        samplingScale = DOFmax.copy();
        samplingScale.assign(DOFmin,Functions.minus);
        samplingScale.assign( Functions.mult(0.1) );*/
        
        //Let's aim for a scaling based on what dimensions allow more flexibility
        //we'll let the scale in each dimension be ~ the distance we can go from our starting point (current x)
        //in that dimension (sum of farthest we can go in positive and negative directions)
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(numDOFs);
        
        for(int dim=0; dim<numDOFs; dim++){
            
            DoubleMatrix1D y = sp.copy();
            //we'll just try and estimate distances within a factor of two
            double upDist = 1e-6;//how far up we can go in this dimension
            double downDist = 1e-6;//how far down
            
            do {
                upDist*=2;
                y.set(dim, sp.get(dim)+upDist);
            } while(checkValidPt(y));
            
            do {
                downDist*=2;
                y.set(dim, sp.get(dim)-downDist);
            } while(checkValidPt(y));
            
            
            ans.set(dim,(upDist+downDist)/6);//divide by 6 so we can make moves within the allowed region
            //though tuneScale should handle uniform scaling of samplingScale better
            //if not adapativeScale..
        }
        
        return ans;
    }

    
    
    
    void burnIn(DoubleMatrix1D startingPoint){
        //burn in the sampler and figure out scale for sampling and how often to use a sample
        
        System.out.println("Starting burn-in for SubThreshSampler.  x="+startingPoint);
        
        x = startingPoint;

        initScale();
        
        if(!adaptiveScale)
            tuneScale();//want to start out with a good scale
                
        //ok now let's burn in until the first and second half of our burn-in periods
        //have nicely overlapping distributions (say, the difference in means is less than 
        //half the st. dev. for each variable (using the less of 1st, 2nd half st. dev.))
        ArrayList<DoubleMatrix1D> burnInSamp = new ArrayList<>();
        DoubleMatrix1D firstHalfSum = DoubleFactory1D.dense.make(numDOFs);
        DoubleMatrix1D secondHalfSum = DoubleFactory1D.dense.make(numDOFs);
        DoubleMatrix1D firstHalfSumSq = DoubleFactory1D.dense.make(numDOFs);
        DoubleMatrix1D secondHalfSumSq = DoubleFactory1D.dense.make(numDOFs);
        
        boolean done = false;

        int nhalf = 0;
        
        for(int b=0; !done; b++){//burn-in samples
            
            while(true) {//draw until acceptable candidate found
                DoubleMatrix1D y = nextCandidate();
                if(checkCandidate(y))//accepting
                    break;
            }
            
            burnInSamp.add(x);
            secondHalfSum.assign(x,Functions.plus);
            DoubleMatrix1D xsq = x.copy().assign(Functions.square);
            secondHalfSumSq.assign(xsq,Functions.plus);
            
            if(b%2==1){//stuff to be done after every pair of samples
                
                //update first and second halves
                DoubleMatrix1D y = burnInSamp.get(b/2);
                firstHalfSum.assign(y,Functions.plus);
                secondHalfSum.assign(y,Functions.minus);
                DoubleMatrix1D ysq = y.copy().assign(Functions.square);
                firstHalfSumSq.assign(ysq,Functions.plus);
                secondHalfSumSq.assign(ysq,Functions.minus);
                
                //now sums are up to date
                //let's try the overlap condition (except right at beginning, hence b>10)
                if(b>10){
                    nhalf = b/2+1;//number of samples in each half
                    DoubleMatrix1D meanDiff = secondHalfSum.copy().assign(firstHalfSum,Functions.minus);
                    meanDiff.assign(Functions.mult(1./nhalf));

                    DoubleMatrix1D std1 = getStDVec(firstHalfSum,firstHalfSumSq,nhalf);
                    DoubleMatrix1D std2 = getStDVec(secondHalfSum,secondHalfSumSq,nhalf);


                    done = true;//overlap condition met
                    for(int dof=0; dof<numDOFs; dof++){
                        if( meanDiff.get(dof) > 0.5*Math.min(std1.get(dof),std2.get(dof)) ){
                            done = false;
                            break;
                        }
                    }
                }
            }
            
            if(done)
                System.out.println("Burn-in complete at sample "+b);
            else if( (b+1)%1000000 == 0 )
                System.out.println("Burn-in sample "+b+" done.  x: "+x);
        }
        
        //we can get a useFrequency from the second half
        DoubleMatrix1D var = getStDVec(secondHalfSum,secondHalfSumSq,nhalf).assign(Functions.square);//full-sequence variance
        DoubleMatrix1D mean = secondHalfSum.copy().assign(Functions.mult(1./nhalf));
        
        //autocorrelation can get a little messy for high shifts
        //but if each dimension's autocorrelation has dropped below 0.1*var for some useFrequency < uf
        //then uf should be a good thinning (minimally correlated samples, for each dimension)
        //we'll use dimGood to keep track of which dimensions have had this
        BitSet dimGood = new BitSet();
        double bestAutorat[] = new double[numDOFs];//to keep track of progress in each dimension (best so far)
        
        
        for(int uf=1; uf<nhalf; uf++){
            
            //calculate autocorrelation (in the sense of covariance of sequence with shifted sequence)
            DoubleMatrix1D autocorr = DoubleFactory1D.dense.make(numDOFs);
            for(int s=nhalf; s<2*nhalf-uf; s++){//start with cross-sum
                DoubleMatrix1D crossTerm = burnInSamp.get(s).copy().assign(mean,Functions.minus);
                DoubleMatrix1D relSamp2 = burnInSamp.get(s+uf).copy().assign(mean,Functions.minus);;
                crossTerm.assign( relSamp2, Functions.mult );
                autocorr.assign( crossTerm, Functions.plus );
            }
            //normalize
            autocorr.assign(Functions.mult(1./(nhalf-uf)));
            
            //for each DOF,
            //autocorr is now the average product of sequence*(sequence shifted by uf)
            //minus the square of the mean of the sequence
            
            DoubleMatrix1D autorat = autocorr.copy().assign( var, Functions.div );
            
            for(int dim=0; dim<numDOFs; dim++){
                if(autorat.get(dim)<0.1)
                    dimGood.set(dim);
                bestAutorat[dim] = Math.min(bestAutorat[dim],autorat.get(dim));
            }
            
            if(dimGood.cardinality()==numDOFs){//all dimensions good
            //if( autorat.zSum()/numDOFs < 0.1 ){//on average, autocorrelation not more than 0.1x variance
                //autocorrelation should ultimately drop to 0, and it may have negative values
                //but it may be comparable to variance for very small uf
                useFrequency = uf;
                System.out.println("Setting useFrequency="+useFrequency+".  Variance: "+var+" Autocorr: "+autocorr);
                break;
            }
                        
            if(uf==nhalf-1){
                useFrequency = uf;
                System.out.println("Warning: high autocorrelation detected at all useFrequencies!  "
                        + "Setting useFrequency="+uf);
            }
            else if( (uf+1)%1000000 == 0 ){
                System.out.print("Trying useFrequency "+uf+
                        ".  Best autocorrelation/variance ratios so far: ");
                for(int dof=0; dof<numDOFs; dof++)
                    System.out.print(bestAutorat[dof]+" ");
                System.out.println();
            }
        }
    }
    
    
    
    DoubleMatrix1D getStDVec(DoubleMatrix1D sum, DoubleMatrix1D sumsq, int nsamp){
        //for n samples with given sum and sum of squares (component-wise)
        //get standard deviations for each component
        DoubleMatrix1D ans = sumsq.copy().assign(Functions.mult(1./nsamp));
        ans.assign( sum.copy().assign(Functions.mult(1./nsamp)).assign(Functions.square), Functions.minus  );
        ans.assign(Functions.sqrt);
        return ans;
    }
    
    
    void tuneScale(){
        int numAccepted = tryScale();
        
        if(numAccepted<numCandScaleTuning*acceptRatioTarget[0]){//too low acceptance ratio 
            while(numAccepted<numCandScaleTuning*acceptRatioTarget[0]){
                samplingScale.assign(Functions.mult(0.7));//scale down sampling scale so we get more acceptances
                numAccepted = tryScale();
            }
        }
        else if(numAccepted>numCandScaleTuning*acceptRatioTarget[1]){//too high
            while(numAccepted>numCandScaleTuning*acceptRatioTarget[1]){
                samplingScale.assign(Functions.mult(1.5));//scale up
                numAccepted = tryScale();
            }
        }
        
        System.out.println("Tuned scale for SubThreshSampler: "+samplingScale);
    }
    
    
    int tryScale(){
        //generate numCandScaleTuning candidates, accept as appropriate, count number accepted
        //this is to try out the current samplingScale
        int numAccepted = 0;
        
        for(int c=0; c<numCandScaleTuning; c++){
            DoubleMatrix1D y = nextCandidate();
            if(checkCandidate(y))//accept
                numAccepted++;
        }
        
        return numAccepted;
    }
    
    
    //once burnt in, we can get our next output sample as many times as needed
    DoubleMatrix1D nextSample(){
        for(int u=0; u<useFrequency; u++){
            while(true) {//draw until acceptable candidate found
                DoubleMatrix1D y = nextCandidate();
                if(checkCandidate(y))//accept
                    break;
            }
        }
        
        return x.copy();
    }
    
    
    
    
    boolean checkCandidate(DoubleMatrix1D y){
        //Metropolis acceptance condition
        //distribution is 1 for checkValidPt, else 0, so only accept valid pts
        //without adaptive scale, samplign distribution is symmetric, this is deterministic
        //else accept according to ratio of sampling distributions
        if(checkValidPt(y)){
            if(adaptiveScale){
                double qforward = Q(x,y,samplingScale);
                DoubleMatrix1D ySamplingScale = getScaleAdaptive(y);
                double qback = Q(y,x,ySamplingScale);
                if(qback/qforward > Math.random())//going to accept
                    samplingScale = ySamplingScale;
                else
                    return false;
            }
            x=y;
            return true;
        }
        else
            return false;
    }
    
    
    double Q(DoubleMatrix1D pt1, DoubleMatrix1D pt2, DoubleMatrix1D scale){
        //(Unnormalized) probability of getting to pt2 from pt1, sampling from a normal distribution
        //with the given standard deviations in each direction (the "sampling scale")
        DoubleMatrix1D diff = pt2.copy().assign(pt1,Functions.minus);

        //ok now we just need a standard normal probability
        //total probability is product of single-variable probabiliies
        double prob = 1;
        for(int dof=0; dof<numDOFs; dof++){
            double v = diff.get(dof)/scale.get(dof);
            prob *= Math.exp(-v*v/2);//standard normal probability, unnormalized
        }
        
        return prob;
    }
    
        
    boolean checkValidPt(DoubleMatrix1D y){
        //is y in the voxel and below the threshold?  (The distribution we're sampling is
        //uniform over the space of valid pts and 0 elsewhere)
        
        //first check voxel bounds...this is cheap
        for(int dof=0; dof<numDOFs; dof++){
            if(y.get(dof)<DOFmin.get(dof) || y.get(dof)>DOFmax.get(dof))
                return false;
        }
        
        if(of.getValue(y)>thresh)
            return false;
        
        return true;
    }
    
    
    DoubleMatrix1D nextCandidate(){
        DoubleMatrix1D y = x.copy();
        for(int dof=0; dof<numDOFs; dof++){
            y.set(dof, x.get(dof)+random.nextGaussian()*samplingScale.get(dof) );
        }
        
        return y;
    }
    
    
}
