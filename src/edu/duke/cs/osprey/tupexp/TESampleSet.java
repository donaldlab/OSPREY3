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

import edu.duke.cs.osprey.confspace.RCTuple;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author mhall44
 */
public class TESampleSet implements Serializable {
    
    
    ArrayList<int[]> samples = new ArrayList<>();//sample assignments
    ArrayList<Double> trueVals = new ArrayList<>();//corresponding true values
    ArrayList<Double> curFitVals = new ArrayList<>();//corresponding current tuple-expansion values
    
    ArrayList<Integer> tupleNumSamples = new ArrayList<>();//number of samples for each tuple
    
    ArrayList<Double> tupleResids = new ArrayList<>();//mean-square residual for samples involving this term 
    double totalResid = 0;//residual over all samples
    
    double worstResid = 0;//the worst residual
    
    ArrayList<Double> sampleResids = new ArrayList<>();

    
    //to keep track of terms that still need drawing
    //handled within constructor
    int numTuplesDone;
    boolean tuplesDone[];
    
    
    TupleExpander te;
    
        
    public TESampleSet(TupleExpander te) {
        
        this.te = te;
        
        //allocate sample arrays for each term...
        for( RCTuple tup : te.tuples ){
            tupleNumSamples.add(0);//no samples for any tuple yet
            tupleResids.add(0.);
        }
        
        
        for(int tup=0; tup<te.tuples.size(); tup++){
            
            if(tup>0 && tup%te.printedUpdateNumTuples==0)
                System.out.println(tup+" tuples done");
            
            updateSamples(tup);
        }
    }
    
       
    
    void updateSamples(int tup){
        //We just raised created or raised the order of approximation for the given tuple terms
        //we need to create more samples for it
        
        int numSampsNeeded = te.numSamplesNeeded(tup) - tupleNumSamples.get(tup);//tupleSamples.get(tup).size();//numDistinctSamples(tup);
        //We may need to have <numSamplesNeeded distinct samples because some tuples have a limited number of possible confs;
        //however, we would overfit if we counted a repeated sample from some other tuple with a small conf space
        // as satisfying the requirement for this tuple
        //(this tuple might have a larger conf space, so this could cause overfitting)
        
        for(int s=0; s<numSampsNeeded; s++){
            int[] sample = drawUnprunedSample(te.tuples.get(tup), true);
            //we know sample won't be null because we checked if tupleFeasible
            
            //don't add if not distinct.  Indicates small conf space
            if( ! isNewSampleDistinct(sample) )
                continue;
            
            double trueVal = te.scoreAssignmentList(sample);
            
            
            ArrayList<Integer> sampTuples = calcSampleTuples(sample);
            sampTuples.trimToSize();

            samples.add(sample);
            
            
            //DEBUG!!!!
            /*for(int sa : sample)
                System.out.print(sa+" ");
            System.out.println();*/

            //update tupleSamples
            for(int term : sampTuples){
                tupleNumSamples.set( term, tupleNumSamples.get(term)+1 );
            }

            trueVals.add(trueVal);
            curFitVals.add(0.);//to be replaced by updateFitVal
        }
    }
    
    
    
    boolean isNewSampleDistinct(int[] sample){
        //Is the new sample distinct from those listed?  (i.e. those that
        //share a tuple and thus could be the same)
       
        for(int sampNum2=0; sampNum2<samples.size(); sampNum2++){
            
            int[] sample2 = samples.get(sampNum2);
            boolean samplesIdentical = true;
                
            for(int pos=0; pos<te.numPos; pos++){
                if(sample[pos]!=sample2[pos]){
                    samplesIdentical = false;
                    break;
                }
            }

            if(samplesIdentical)
                return false;
        }
        
        //if we get here, it's distinct
        return true;
    }
    
            
    boolean tupleFeasible(RCTuple tup){
        int sample[] = drawUnprunedSample(tup, false);
        return (sample!=null);//were we able to draw a sample?
    }
    
    
    boolean finishSample(int[] sample){
        //Fill in the -1's in the sample
        //unlike in regular TESampleSet, these assignments are randomly selected
        //just avoiding pruned assignments or pairs
        //we return whether we did this without running out of possibilities for some res position
        //(this kind of failure requires a restart)
        
        
        while(true){
            
            int nextPos = getRandomUnfilledPos(sample);
            
            if(nextPos==-1)//all positions filled in
                return true;
            
            //now draw a random unpruned assignment for nextPos
            int distr[] = new int[te.numAllowed[nextPos]];
            boolean haveOptions = false;//have at least one option
            
            for(int a=0; a<te.numAllowed[nextPos]; a++){
                distr[a] = 1;

                if(te.canCheckPartialPruning) {
                    if (te.isPruned(new RCTuple(nextPos, a)))
                        distr[a] = 0;

                    //look for pruned pairs between a and assignments already in sample
                    for (int pos = 0; pos < te.numPos; pos++) {
                        if (sample[pos] != -1) {

                            if (isPairPrunedInSample(pos, sample[pos], nextPos, a, sample)) {
                                distr[a] = 0;
                                break;
                            }
                        }
                    }
                }
                else {
                    int augSamp[] = sample.clone();
                    augSamp[nextPos] = a;
                    if(te.isPruned(new RCTuple(augSamp)))
                        distr[a] = 0;
                }
                
                
                if(distr[a]==1)
                    haveOptions = true;
            }
            
            if(haveOptions)
                sample[nextPos] = discreteDraw(distr);
            else//no unpruned options
                return false;
        }
        
    }
    
    
    int getRandomUnfilledPos(int sample[]){
        //random position set to -1 in sample
        //otherwise uniformly distributed
        int distr[] = new int[te.numPos];
        boolean haveOptions = false;//some position is unfilled
        
        for(int pos=0; pos<te.numPos; pos++){
            if(sample[pos]==-1){
                distr[pos] = 1;
                haveOptions = true;
            }
        }
        
        if(haveOptions)
            return discreteDraw(distr);
        else
            return -1;//signals all positions filled
    }
    
    
    
    static int discreteDraw(int[] distr){
        //draw an integer m with probability proportional to distr[m]
        
        int tot = 0;
        for(int a : distr)
            tot += a;
        
        int randVal = new java.util.Random().nextInt(tot);
        int cumSum = 0;
        
        for(int ans=0; ans<distr.length; ans++){
            cumSum += distr[ans];
            if(cumSum>randVal)
                return ans;
        }
        
        throw new RuntimeException("ERROR: shouldn't get to end of discreteDraw!!!");
    }
    
    
    
    ArrayList<Integer> calcSampleTuples(int[] sample){
        //calculate the tuples involved in a given sample
        ArrayList<Integer> sampTuples = new ArrayList<>();
        
        for(int tup=0; tup<te.tuples.size(); tup++){
            
            if(te.sampleMatchesTuple(sample, te.tuples.get(tup))){
                sampTuples.add(tup);
            }
        }
        
        return sampTuples;
    }
    
    
    
    
    //update all the fit values
    void updateFitVals(FittingObjFcn fof){
        for(int s=0; s<samples.size(); s++)
            //curFitVals.set( s, te.fitValueForTuples(sampleTuples.get(s)) );
            curFitVals.set( s, te.fitValueForTuples(calcSampleTuples(samples.get(s))) );
        
        //and update residuals accordingly
        updateAllResids(fof);
    }
    
    
    
    
    void updateAllResids(FittingObjFcn fof){
        
        totalResid = 0;
        for(int t=0; t<te.tuples.size(); t++)
            tupleResids.set(t,0.);
        
        sampleResids = fof.computeAllResids(trueVals, curFitVals, te.constTerm);
        
        for(int s=0; s<samples.size(); s++){
            //double targetVal = trueVals.get(s) - te.constTerm;
            //double sampResid = fof.computeResid(curFitVals.get(s)-te.constTerm, targetVal);
            double sampResid = sampleResids.get(s);
            totalResid += sampResid;
            for(int tup : calcSampleTuples(samples.get(s)) )
                tupleResids.set( tup, tupleResids.get(tup)+sampResid );
            
            //sampleResids.add(sampResid);
        }
        
        totalResid /= samples.size();
        worstResid = 0;
        for(int t=0; t<te.tuples.size(); t++){
            tupleResids.set(t,tupleResids.get(t)/tupleNumSamples.get(t));//tupleSamples.get(t).size());
            worstResid = Math.max(worstResid,tupleResids.get(t));
        }
    }
    
    
    
    void addTuple(int tup){
        //add a tuple that has just been put into TupleExpander, updating samples accordingly
        //tup will be tuples.size()-1
        
        //record relationship to current samples
        tupleNumSamples.add(0);
        
        for(int s=0; s<samples.size(); s++){
            int sample[] = samples.get(s);
            
            if(te.sampleMatchesTuple(sample,te.tuples.get(tup))){
                tupleNumSamples.set( tup, tupleNumSamples.get(tup)+1 );
            }
        }
        
        
        //preallocate residuals
        tupleResids.add(0.);
        
        //draw any new samples needed
        updateSamples(tup);
    }
    
    
    void printResids(){
        System.out.println(samples.size()+" samples");
        System.out.println("Total resid: "+totalResid+" Worst resid: "+worstResid);
    }
    
    
    
    
    int[] drawUnprunedSample(RCTuple tup, boolean errorIfImpossible){
        //draw an unpruned sample containing the specified tuple
        //errorIfImpossible means we expect this to be possible, and should throw an error if it isn't
        //(return null if impossible and no error)
        
        //first we'll try drawing randomly...less overhead
        int[] sample = new int[te.numPos];
        
        
        
        //DEBUG!!!
        long startTime = System.currentTimeMillis();
        
        
        for(int tryNum=0; tryNum<5; tryNum++){
            Arrays.fill(sample, -1);
            te.assignTupleInSample(sample, tup);
            boolean finishSampleSuccessful = finishSample(sample);
            
            
            //DEBUG!!!
            long sampTime = System.currentTimeMillis() - startTime;
            //taking over 10 s is going to be an issue
            if(sampTime > 10000){
                System.out.println();
                System.out.println("finishSample sampling took over 10 s (ms shown): "+sampTime);
                System.out.println("Target tuple: "+tup.stringListing());
                System.out.println("tryNum: "+tryNum);
                if(finishSampleSuccessful)
                    System.out.println("Sample: "+new RCTuple(sample).stringListing());
                else
                    System.out.println("No sample found yet.");
                System.out.println();
            }
            //DEBUG!!
        
        
        
            
            if( finishSampleSuccessful )//finished successfully
                return sample;
        }
        
        //ok let's try something more systematic.  Depth-first search.  
        Arrays.fill(sample,-1);
        te.assignTupleInSample(sample, tup);
        
        //DEBUG!!  Want to see what's taking so long on second-round 30.5, 40.cont
        long sampStartTime = System.currentTimeMillis();
        
        
        sample = finishSampleDFS(sample);
        
        
        long sampTime = System.currentTimeMillis() - sampStartTime;
        //taking over 10 s is going to be an issue
        if(sampTime > 10000){
            System.out.println();
            System.out.println("Sampling took over 10 s (ms shown): "+sampTime);
            System.out.println("Target tuple: "+tup.stringListing());
            if(sample != null)
                System.out.println("Sample: "+new RCTuple(sample).stringListing());
            else
                System.out.println("No sample found.");
            System.out.println();
        }
        
        
        if(errorIfImpossible && sample==null){
            throw new RuntimeException("ERROR: No unpruned samples available for tuple "
                    + tup.stringListing());   
        }
        
        return sample;
    }
    
        
    int[] finishSampleDFS(int[] sample){
        //Complete the sample by DFS
        //For this we set up the list of RCs available at each position.  
        //Each time we'll search the position with the most eliminated RCs (by pairs+ pruning)
        //this way hopefully can clear out tuples that need pruning relatively quickly
        //random ordering of assignments at each position
        
        ArrayList<ArrayList<int[]>> allowedRCs = new ArrayList<>();
        ArrayList<Integer> numElim = new ArrayList<>();
        listCompatibleRCs(allowedRCs,numElim,sample);
        
        boolean tuplePossible = true;
        
        //see if any positions had no compatible RCs...indicates our tuple is impossible
        for(int pos=0; pos<te.numPos; pos++){
            if(sample[pos]==-1 && allowedRCs.get(pos).isEmpty()){
                tuplePossible = false;
                break;
            }
        }
        
        if(tuplePossible)
            tuplePossible = unprunedSampDFS(sample, allowedRCs, numElim);
               
        if(tuplePossible)
            return sample;
        else
            return null;
    }
    
    
    
    boolean unprunedSampDFS(int sample[], ArrayList<ArrayList<int[]>> allowedRCs, ArrayList<Integer> numElim){
        
        //first pick the level to assign next
        int nextPos = -1;
        double nextPosElimRatio = -1;//ratio of eliminated to non-eliminated RCs at positions
        //higher ratio --> more interactions --> good to expand next
        
        for( int pos=0; pos<te.numPos; pos++ ){
            if( sample[pos] == -1 ){//sample still open to be assigned
                //see how many unpruned RCs are eliminated by this...
                double elimRatio = 1.0*numElim.get(pos)/allowedRCs.get(pos).size();
                if(elimRatio > nextPosElimRatio){
                    nextPos = pos;
                    nextPosElimRatio = elimRatio;
                }
            }
        }
        
        if(nextPos==-1){//everything assigned already!
            return true;
        }
        else {//need to assign nextPos
            ArrayList<Integer> newOptions = shuffleOptions(allowedRCs.get(nextPos));
            
            for(int opt : newOptions){
                
                sample[nextPos] = opt;
                
                //determine what new options are eliminated at each unassigned position
                ArrayList<int[]> elimHere = new ArrayList<>();
                for(int pos=0; pos<te.numPos; pos++){
                    if(sample[pos]==-1){
                        
                        for(int q=allowedRCs.get(pos).size()-1; q>=0; q--){//reverse iteration to allow deletion
                            int[] uaOpt = allowedRCs.get(pos).get(q);

                            boolean pruned;
                            if(te.canCheckPartialPruning)
                                pruned = isPairPrunedInSample(nextPos,opt,uaOpt[0],uaOpt[1],sample);
                            else {
                                int augSamp[] = sample.clone();
                                augSamp[uaOpt[0]] = uaOpt[1];
                                pruned = te.isPruned(new RCTuple(augSamp));
                            }

                            if(pruned){
                                elimHere.add(uaOpt);
                                allowedRCs.get(pos).remove(q);
                                numElim.set( pos, numElim.get(pos)+1 );
                            }
                        }
                    }
                }
                                
                
                boolean success = unprunedSampDFS(sample,allowedRCs,numElim);
                //DFS, so if we were successful then we need not search any more
                if(success)
                    return true;
                
                
                //else restore options eliminated here and go on to the next option...
                for(int[] uaOpt : elimHere){
                    allowedRCs.get(uaOpt[0]).add(uaOpt);
                    numElim.set( uaOpt[0], numElim.get(uaOpt[0])-1 );
                }
                
                //allowedRCs and numElim are kept up to date until the position is assigned, then kept as is
                //as we go deeper
            }
            
            //if we get here, nothing worked in this branch, so return false
            sample[nextPos] = -1;//need to indicate that position is unassigned
            return false;
        }
    }
    
    
    
    private void listCompatibleRCs(ArrayList<ArrayList<int[]>> allowedRCs, 
            ArrayList<Integer> numElim, int[] sample){
        //list RCs at all unassigned positions compatible with the assigned positions in the sample

        for(int pos=0; pos<te.numPos; pos++){
                        
            ArrayList<int[]> forPos = new ArrayList<>();
            int numElimAtPos = 0;
            
            if(sample[pos]==-1){//need to identify options
                
                for(int rc=0; rc<te.numAllowed[pos]; rc++){
                    if(!te.isPruned(new RCTuple(pos,rc))){

                        boolean incompatible = false;
                        if(te.canCheckPartialPruning){
                            for(int pos2=0; pos2<te.numPos; pos2++) {
                                if (sample[pos2] != -1) {//pos2 assigned
                                    if (isPairPrunedInSample(pos2, sample[pos2], pos, rc, sample)) {
                                        incompatible = true;
                                        break;
                                    }
                                }
                            }
                        }
                        else {
                            int[] augSamp = sample.clone();
                            augSamp[pos] = rc;
                            incompatible = te.isPruned(new RCTuple(augSamp));
                        }


                        if(incompatible)
                            numElimAtPos++;
                        else
                            forPos.add(new int[] {pos,rc});
                    }
                }
            }
            
            allowedRCs.add(forPos);
            numElim.add(numElimAtPos);
        }
    }
    
    
    static ArrayList<Integer> shuffleOptions(ArrayList<int[]> optList){
        //given a bunch of (pos,rc) pairs, return the rc's in random order
        ArrayList<Integer> ans = new ArrayList<>();
        int numRCs = optList.size();
        
        ArrayList<Integer> ordering = randomOrdering(numRCs);
        
        for(int a=0; a<numRCs; a++){
            ans.add( optList.get(ordering.get(a))[1] );
        }
        
        return ans;
    }
    
    static ArrayList<Integer> randomOrdering(int n){
        //random order of integers [0,n)
        ArrayList<Integer> a = new ArrayList<>();
        ArrayList<Integer> ans = new ArrayList<>();
        
        Random rand = new Random();
        
        for(int b=0; b<n; b++)
            a.add(b);
        
        for(int b=0; b<n; b++){
            ans.add( a.remove(rand.nextInt(n-b)) );
        }
        
        return ans;
    }
    
    boolean isPairPrunedInSample(int pos1, int rc1, int pos2, int rc2, int[] sample){
        //check if the specified pair is pruned, or is part of a pruned higher tuple involving
        //the assigned options in the given sample
        //(pos1,rc1) is in the sample; (pos2,rc2) is not
        RCTuple pair = new RCTuple(pos1, rc1, pos2, rc2);
            
        if(te.isPruned(pair))
            return true;

        //and any higher-order pruned stuff...
        for(RCTuple prunedTup : te.higherOrderPrunedTuples(pair)){
            
            sample[pos2] = rc2;
            boolean pruned = te.sampleMatchesTuple(sample,prunedTup);
            sample[pos2] = -1;//revert to sample as before
            
            if(pruned)
                return true;
        }
        
        return false;
    }
    
    
}
