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

package edu.duke.cs.osprey.voxq;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.math.Functions;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * Calculate the energy difference between two voxels, using BAR-type method
 * 
 * @author mhall44
 */
public class VoxelsDeltaG {
    
    static int sampleBatchSize = 100;//samples will be drawn in batches of this size
    
    int numDOFs;
    //samples from each voxel
    ArrayList<Sample> samples1, samples2;
    //the samplers that generate them
    IntraVoxelSampler sampler1, sampler2;
    
    SampleNormalization sn1, sn2;//normalizations for the two voxels
    
    //current estimates of delta G (and relative error in integrals)
    double estDeltaG = 0;//0 is initial guess
    double integRelErr1 = Double.POSITIVE_INFINITY;
    double integRelErr2 = Double.POSITIVE_INFINITY;
    
    
    
    public VoxelsDeltaG(MoleculeModifierAndScorer mms1, MoleculeModifierAndScorer mms2, boolean alignByEnergy){
        //if alignByEnergy we try to align low-energy regions
        //else we assume aligning voxel bounds suffices
        //just need to set sampler1, sampler2, numDOFs
        //System.out.println("INITIALIZING IVS1");
        sampler1 = new IntraVoxelSampler(mms1);
        //System.out.println("INITIALIZING IVS2");
        sampler2 = new IntraVoxelSampler(mms2);
        //System.out.println("DONE INITIALIZING IVS");
        numDOFs = sampler1.numDOFs;
        if(sampler2.numDOFs!=numDOFs)
            throw new RuntimeException("ERROR: Not supporting delta G for voxels w/ different # DOFs currently...");
        
        
        //SETTING UP NORMALIZATIONS
        
        
        if(alignByEnergy){
                //first, draw some initial samples
            ArrayList<DoubleMatrix1D> fullSamples1 = new ArrayList<>();//will need full samples for alignment
            ArrayList<DoubleMatrix1D> fullSamples2 = new ArrayList<>();
            for(int n=0; n<sampleBatchSize; n++){
                //System.out.println("DRAWING SAMP 1");
                DoubleMatrix1D samp1 = sampler1.nextSample();
                fullSamples1.add(samp1);
                //System.out.println("DRAWING SAMP 2");
                DoubleMatrix1D samp2 = sampler2.nextSample();
                fullSamples2.add(samp2);
            }
            //and use these initial samples to figure out what alignment we want
            sn1 = new SampleNormalization(fullSamples1);
            sn2 = new SampleNormalization(fullSamples2);
        }
        else {
            sn1 = new SampleNormalization(mms1.getConstraints());
            sn2 = new SampleNormalization(mms2.getConstraints());
        }
    }
    
    
    
    
    private class Sample {
        double Ediff;//difference between voxel 2 and voxel 1 energies
        //at the point in DOF space corresponding to this sample
        //DEBUG!!!  should also have jacRatio (dz2/dy)/(dz1/dy)
        
        Sample(DoubleMatrix1D DOFVals, boolean isVox1){
            //generate sample given DOF values and whether they're drawn from voxel 1 or 2
            DoubleMatrix1D z1, z2;//corresponding points in voxels 1 and 2
            if(isVox1){
                z1 = DOFVals;
                z2 = sn2.unnormalize(sn1.normalize(DOFVals));
                
                if(sampler2.mms.isOutOfRange(z2)){//energies outside voxel considered infinite
                    Ediff = Double.POSITIVE_INFINITY;
                    return;
                }
            }
            else {
                z2 = DOFVals;
                z1 = sn1.unnormalize(sn2.normalize(DOFVals));
                
                if(sampler1.mms.isOutOfRange(z1)){//energies outside voxel considered infinite
                    Ediff = Double.NEGATIVE_INFINITY;
                    return;
                }
            }
            
            Ediff = sampler2.mms.getValue(z2) - sampler1.mms.getValue(z1);
        }
    }
    
    
    private class SampleNormalization {
        DoubleMatrix1D center;
        DoubleMatrix1D scaling;
        double jacDet;//determinant of Jacobian of (linear) mapping from normalized to actual DOFs
        
        SampleNormalization(DoubleMatrix1D[] constr){
            //just normalize so that each normalized DOF runs from -.5 to .5 between the given constraints
            center = constr[0].copy();
            center.assign(constr[1],Functions.plus);
            center.assign(Functions.mult(0.5));
            
            scaling = constr[1].copy();
            scaling.assign(constr[0],Functions.minus);
            
            jacDet = 1;
            for(double el : scaling.toArray())
                jacDet *= el;
        }
        
        SampleNormalization(ArrayList<DoubleMatrix1D> fullSamples){
            center = DoubleFactory1D.dense.make(numDOFs);
            scaling = DoubleFactory1D.dense.make(numDOFs);
            jacDet = 1;
            
            for(int dofNum=0; dofNum<numDOFs; dofNum++){
                ArrayList<Double> vals = new ArrayList<>();
                int numSamples = fullSamples.size();
                double cen = 0;
                for(DoubleMatrix1D samp : fullSamples){
                    vals.add(samp.get(dofNum));
                    cen += samp.get(dofNum);
                }
                center.set(dofNum, cen/numSamples);
                Collections.sort(vals);
                //estimate spread using interquartile range
                double sc = vals.get(3*numSamples/4) - vals.get(numSamples/4);
                jacDet *= sc;
                scaling.set(dofNum, sc);
            }
        }
        
        DoubleMatrix1D unnormalize(DoubleMatrix1D y){//normalized to unnormalized DOF vals
            return y.copy().assign(scaling, Functions.mult).assign(center, Functions.plus);
        }
        DoubleMatrix1D normalize(DoubleMatrix1D z){//unnormalized to normalized DOF vals
            return z.copy().assign(center, Functions.minus).assign(scaling, Functions.div);
        }
    }
    
    
    public double estDeltaG(double stdErr){
        
        double integRelErrTarget = stdErr / IntraVoxelSampler.RT;//desired relative error for each integral
                
        //now can add these samples to our lists
        samples1 = new ArrayList<>();
        samples2 = new ArrayList<>();
        for(int n=0; n<sampleBatchSize; n++){
            //samples1.add(new Sample(fullSamples1.get(n), true));
            //samples2.add(new Sample(fullSamples2.get(n), false));
            //WATCH OUT USING SAME SAMPLES TO DO NORMALIZATION & ENERGY CAUSES BIAS
            //ANY NORMALIZATION IS FINE BUT MUST BE INDEPENDENT OF ENERGY SAMPLES
            samples1.add(new Sample(sampler1.nextSample(), true));
            samples2.add(new Sample(sampler2.nextSample(), false));
        }
        
        
        //Next, estimate the integral based on said alignment
        //will need to iterate to converge self-consistently on a delta-E estimate
        while(true){
            //estimate energy from samples so far
            double newDeltaG = curDeltaGEstimate();
            
            //see if converged
            if(Math.abs(newDeltaG-estDeltaG)<stdErr && totIntegRelErr()<integRelErrTarget){
                return newDeltaG;//this estimate is good
            }
            else{
                estDeltaG = newDeltaG;
                //draw a new batch of samples from each voxel
                for(int n=0; n<sampleBatchSize; n++){
                    samples1.add(new Sample(sampler1.nextSample(), true));
                    samples2.add(new Sample(sampler2.nextSample(), false));
                }
            }
        }
    }
    
    
    private double totIntegRelErr(){
        return Math.sqrt(integRelErr1*integRelErr1 + integRelErr2*integRelErr2);
    }
    
    
    private static double fd(double E){//Fermi-dirac distribution
        return 1./(1+Math.exp(E/IntraVoxelSampler.RT));
    }
    
    static double mean(ArrayList<Double> arr){
        double ans = 0;
        for(double a : arr)
            ans += a;
        return ans / arr.size();
    }
    
    static double relStdDev(ArrayList<Double> arr, double mean){
        //Given arr and its mean, return standard deviation / mean
        double ans = 0;
        for(double a : arr)
            ans += (a-mean)*(a-mean);
        ans /= (arr.size()-1);
        ans = Math.sqrt(ans);
        ans /= mean;
        return ans;
    }
    
    double curDeltaGEstimate(){
        //Estimate delta G using BAR and current samples
        //return estimate, set relative errors for the two integrals we compute (integRelErrs)
        
        //compute the integrals
        //integrals are averages of f1, f2
        ArrayList<Double> f1 = new ArrayList<Double>();
        ArrayList<Double> f2 = new ArrayList<Double>();
        for(Sample s : samples1)//DEBUG!!!  fd should be multiplied by sqrt(s.jacRatio)
            f1.add(fd(s.Ediff-estDeltaG));
        for(Sample s : samples2)
            f2.add(fd(estDeltaG-s.Ediff));
        
        double integ1 = mean(f1);
        double integ2 = mean(f2);
        integRelErr1 = relStdDev(f1,integ1) / Math.sqrt(f1.size());
        integRelErr2 = relStdDev(f2,integ2) / Math.sqrt(f2.size());
        
        return estDeltaG - IntraVoxelSampler.RT * ( Math.log(integ1*sn2.jacDet) - Math.log(integ2*sn1.jacDet) );
    }
    
    public int numSamplesNeeded(){
        //if we just called estDeltaG, this will let us see how many samples were needed
        //should be same for samples1 and samples2
        return samples1.size();
    }
    
    
}
