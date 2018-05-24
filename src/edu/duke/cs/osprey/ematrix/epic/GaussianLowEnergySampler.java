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

package edu.duke.cs.osprey.ematrix.epic;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Random;
import org.apache.commons.math3.special.Erf;
import org.ojalgo.random.RandomUtils;
import static org.ojalgo.random.RandomUtils.erf;

/**
 *
 * This sampler samples preferentially from the low-energy region about the center
 * Used when this region is too small to be efficiently reached for this purpose
 * It samples Gaussianly about the center,
 * with variance determined by initial explorations
 * 
 * @author mhall44
 */
public class GaussianLowEnergySampler {
    
    double EPICThresh1;
    ObjectiveFunction of;
    DoubleMatrix1D DOFmin, DOFmax;
    //we want our samples x to lie in the bounds given by DOFmin, DOFmax
    //and to have of.get(x)<=thresh
    int numDOFs;
        
    //We'll draw from a Gaussian distribution about center with given sigmas
    //chosen to favor good region
    DoubleMatrix1D center;
    DoubleMatrix1D sigmas;
    
    Random random = new Random();//a random
            

    public GaussianLowEnergySampler(double thresh, ObjectiveFunction of, DoubleMatrix1D DOFmin, 
            DoubleMatrix1D DOFmax, DoubleMatrix1D center) {
        this.EPICThresh1 = thresh;
        this.of = of;
        this.DOFmin = DOFmin;
        this.DOFmax = DOFmax;
        this.center = center;
        numDOFs = DOFmin.size();
        
        //OK we need to figure out sigmas based on thresh
        
        
        double sigmaVal = chooseNumSigmas();
        
        //now go along axes to find sigma
        sigmas = DoubleFactory1D.dense.make(numDOFs);
        for(int dofNum=0; dofNum<numDOFs; dofNum++){
                        
            //bound the good region along dofNum axis
            double goodRegionLB = bisectForGoodRegionEnd(dofNum, DOFmin.get(dofNum));
            double goodRegionUB = bisectForGoodRegionEnd(dofNum, DOFmax.get(dofNum));
            
            //sigma will based on the farther of these from center
            double sigma;
            if( goodRegionUB - center.get(dofNum) >= center.get(dofNum) - goodRegionLB ){
                sigma = (goodRegionUB - center.get(dofNum)) / sigmaVal;
            }
            else
                sigma = (center.get(dofNum) - goodRegionLB) / sigmaVal;
            
            sigmas.set(dofNum, sigma);
        }
    }
    
    
    private double chooseNumSigmas(){
        //The good region will extend to roughly
        //this many sigmas in each dimension
        
        double insideTarget = 0.25;//how much of the multivariate distr should be in the "good zone"
        //underestimate since we can always reject to make up for uncertainty
        
        
        //OK zt will be the value of z = sum_i x_i^2 / 2sigma_i^2
        //such that integrating the probability of a Gaussian over the z>zt zone = insideTarget
        //This will tell us at what zt to place our threshold points
        //WLOG we do this integration with a standard normal distribution in numDOFs dimensions
        //we need consider only the radial integral: 
        //\int_0^a exp(-r^2/2) r^(numDOFs-1) dr
        //then insideTarget = radialIntegral(0 to sqrt(2*zt)) / radialIntegral(0 to infty)
        
        
        double fullInteg = 1;//radial integral from 0 to infty
        if(numDOFs%2==1)
            fullInteg = Math.sqrt(Math.PI/2);
        for(int k=numDOFs-2; k>0; k-=2)
            fullInteg *= k;
        
        
        //we find the right sqrt(2*zt) by bisection
        double target = fullInteg * insideTarget;//value we want radial integral to achieve
        double epsilon = 0.02;
        
        double lo = 0;//initial lower bound on sqrt(2*zt)
        double hi = 1;//upper bound
        while(radialIntegral(hi,numDOFs-1)<target)
            hi *= 2;
        
        do {
            double mid = (hi+lo)/2;
            double midVal = radialIntegral(mid,numDOFs-1);
            if( midVal < target )
                lo = mid;
            else
                hi = mid;
        } while( Math.abs(hi-lo) > epsilon*Math.abs(hi) );
        
        
        //OK so (hi+lo)/2 is now a good estimate of sqrt(2*zt)
        double numSigmas = (hi+lo)/2;
        //Iff the point where the i-axis crosses thresh is y_i, then 
        //sigma_i = y_i/sqrt(2*zt), so we set numSigmas = est(sqrt(2*zt)) = (hi+lo)/2
        return numSigmas;
    }
    
    
    private double radialIntegral(double a, int k){
        //\int_0^a exp(-r^2/2) r^k dr
        if(k>1)
            return (k-1)*radialIntegral(a,k-2) - Math.pow(a,k-1)*Math.exp(-a*a/2);
        else if(k==1)
            return 1-Math.exp(-a*a/2);
        else if(k==0)
            return Math.sqrt(Math.PI/2) * Erf.erf(a/Math.sqrt(2));
        else
            throw new RuntimeException("ERROR: Negative k not supported here: "+k);
    }
    
   
    private double bisectForGoodRegionEnd(int dofNum, double voxEnd){
        //Starting at center, see how far we can go towards voxEnd
        //in dimension dofNum
        //without surpassing the threshold
        of.setDOFs(center);
        double centerE = of.getValForDOF(dofNum, center.get(dofNum));
        if(of.getValForDOF(dofNum, voxEnd) <= centerE + EPICThresh1){
            //voxEnd still within thresh
            return voxEnd;
        }
        else {//must bisect
            //current bounds
            double hi = voxEnd;
            double lo = center.get(dofNum);
            double epsilon = 0.1;//relative error thresh
            
            do {
                double mid = (hi+lo)/2;
                double midVal = of.getValForDOF(dofNum,mid);
                if( midVal <= centerE+EPICThresh1 )
                    lo = mid;
                else
                    hi = mid;
            } while( Math.abs(hi-lo) > epsilon*Math.abs(hi-center.get(dofNum)) );
            
            return hi;
        }
    }
    
    
    DoubleMatrix1D nextSample(){//Sample from the distribution
        DoubleMatrix1D ans = DoubleFactory1D.dense.make(numDOFs);
        for(int dofNum=0; dofNum<numDOFs; dofNum++){
            do {
                double x = center.get(dofNum) + sigmas.get(dofNum) * random.nextGaussian();
                ans.set(dofNum, x);
            } while( ans.get(dofNum) < DOFmin.get(dofNum) 
                    || ans.get(dofNum) > DOFmax.get(dofNum) );
        }
        
        return ans;
    }    
    
}
