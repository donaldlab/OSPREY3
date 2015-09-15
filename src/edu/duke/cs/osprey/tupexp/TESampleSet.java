/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.RCTuple;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.PriorityQueue;
import java.util.Random;

/**
 *
 * @author mhall44
 */
public class TESampleSet implements Serializable {
    
    
    ArrayList<int[]> samples = new ArrayList<>();//sample assignments
    ArrayList<Double> trueVals = new ArrayList<>();//corresponding true values
    ArrayList<Double> curFitVals = new ArrayList<>();//corresponding current tuple-expansion values
    
    //ArrayList<ArrayList<Integer>> tupleSamples = new ArrayList<>();//for each tuple, lists samples involving this tuple
    ArrayList<Integer> tupleNumSamples = new ArrayList<>();//number of samples for each tuple
    
    ArrayList<BitSet> tupleSamplesAboveCutoff = new ArrayList<>();//which of the tupleSamples are above the bCutoff
    //ArrayList<ArrayList<Integer>> sampleTuples = new ArrayList<>();//for each sample, list tuples it involves
    ArrayList<Double> tupleResids = new ArrayList<>();//mean-square residual for samples involving this term 
    double totalResid = 0;//residual over all samples
    
    double worstResid = 0;//the worst residual
    
    ArrayList<Double> sampleResids = new ArrayList<>();

    
    //to keep track of terms that still need drawing
    //handled within constructor
    int numTuplesDone;
    boolean tuplesDone[];
    
    
    TupleExpander te;
    
    
    //we may need to repeat some drawing to get samples under bCutoff
    //we'll redraw sampleRCs up to maxDiscReps times
    //if that fails, then we probably can't get all the terms we need with the current bCutoff,
    //so we raise bCutoff to clear the minimum energy we got in all these reps
    double minEOverReps;
    
    static int maxDiscrReps = 100;
    
    
    int tupleFailCounts[];//if using tuple bCutoffs, then use this to decide when to raise
    //it keeps track of how many samples with the given tuple have failed to yield a sub-bCutoff
    //true value since the last success

    
    
    ArrayList<Double> sampleBCutoffs = new ArrayList<>();
    
    //DEBUG!!!!!
    //trying to find the bottleneck in sampling
    /*int numNoIter = 0;
    int numNoLoop = 0;
    int numNeededSamples = 0;*///one sample per tuple, maxing out at 10 per
    //should go up at every draw
    
    
    
    public TESampleSet(TupleExpander te) {
        
        this.te = te;
        
        //allocate sample arrays for each term...
        for( RCTuple tup : te.tuples ){
            //tupleSamples.add( new ArrayList<Integer>() );
            tupleNumSamples.add(0);//no samples for any tuple yet
            
            tupleSamplesAboveCutoff.add(new BitSet());
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
        
        
        //DEBUG!!!
        //This simpler version should be more overfitting-resistant,
        //probably will end up ditching bcutoff (pruning fills that function) and just using this
        if(te.bCutoff == Double.POSITIVE_INFINITY){
            updateSamples2(tup);
            //updateSamplesMultithread(tup);
            return;
        }
        
        int numSampsNeeded = te.numSamplesNeeded(tup) - tupleNumSamples.get(tup);//tupleSamples.get(tup).size();
        int numSubCutoffSampsNeeded = te.numSamplesNeeded(tup)/2 - 
                //(tupleSamples.get(tup).size()-tupleSamplesAboveCutoff.get(tup).cardinality());
                (tupleNumSamples.get(tup)-tupleSamplesAboveCutoff.get(tup).cardinality());
        
        
        if(te.subCutoffSamplesOnly)
            numSubCutoffSampsNeeded = numSampsNeeded;
        
        
        int maxSampsDrawn = maxDiscrReps*Math.max(numSampsNeeded,numSubCutoffSampsNeeded);
        //the most samples we're willing to draw
        
        PriorityQueue<ConfPair> bestAboveCutoff = new PriorityQueue<>();
        //best samples above cutoff so far
        //note the head of this queue is its highest energy member!  This ensures that we only keep 
        //at most numSubCutoffSampsNeeded confs in there
        
        for(int s=0; numSampsNeeded>0 || numSubCutoffSampsNeeded>0; s++){
                        
            
            int[] sample = drawUnprunedSample(te.tuples.get(tup), true);
            /*int[] sample = new int[te.numPos];

            //fill the current tuple into RCSample
            //for(int[] op : te.tuples.get(tup))
            //    sample[op[0]] = op[1];

            //if term isn't pruned there must be some way to finishRCSample
            boolean success;
            do{
                Arrays.fill(sample, -1);
                te.assignTupleInSample(sample, te.tuples.get(tup));
                success = finishSample(sample);
            }
            while(!success);*/

            //figure out which terms the sample applies to
            double trueVal = te.scoreAssignmentList(sample);
            
            if(trueVal<=te.bCutoff || numSampsNeeded>numSubCutoffSampsNeeded){
                //we'll take this sample now if it is below the standard bCutoff
                //or we still can use over-bcutoff samples
                ArrayList<Integer> sampTuples = calcSampleTuples(sample);


                //sampleTuples.add(sampTuples);
                samples.add(sample);
                int sampleIndex = samples.size()-1;

                //update tupleSamples
                for(int term : sampTuples){
                    //tupleSamples.get(term).add(sampleIndex);
                    tupleNumSamples.set( term , tupleNumSamples.get(term)+1 );
                    
                    if(trueVal>te.bCutoff)
                        //tupleSamplesAboveCutoff.get(term).set(tupleSamples.get(term).size()-1);
                        tupleSamplesAboveCutoff.get(term).set(tupleNumSamples.get(term)-1);
                }
                

                trueVals.add(trueVal);
                curFitVals.add(0.);//to be replaced by updateFitVal
                
                sampleBCutoffs.add(te.bCutoff);
            }
            else {
                bestAboveCutoff.add(new ConfPair(sample,new double[] {trueVal}));
                //we only need the best numSubCutoffSampsNeeded over-cutoff samples to be stored here
                while(bestAboveCutoff.size()>numSubCutoffSampsNeeded)
                    bestAboveCutoff.remove();
            }
            
            numSampsNeeded = te.numSamplesNeeded(tup) - tupleNumSamples.get(tup);//tupleSamples.get(tup).size();
            numSubCutoffSampsNeeded = te.numSamplesNeeded(tup)/2 - 
                //(tupleSamples.get(tup).size()-tupleSamplesAboveCutoff.get(tup).cardinality());
               (tupleNumSamples.get(tup)-tupleSamplesAboveCutoff.get(tup).cardinality());
            
            if(te.subCutoffSamplesOnly)
                numSubCutoffSampsNeeded = numSampsNeeded;
            
            if(s==maxSampsDrawn){
                //maxed out on samples to draw...indicates we're probably not getting
                //enough samples below standard bcutoff for this tuple
                //so introduce raised-bcutoff samples to avoid overfitting
                
                if(numSampsNeeded>numSubCutoffSampsNeeded)//should have drawn way more than enough over-cutoff samples by now
                    //since that's what drives us to max out...
                    throw new RuntimeException("ERROR: Maxed out on sampling in TESampleSet2 but still need over-cutoff samples...");
                
                //remove any but the best above-cutoff samples
                while(bestAboveCutoff.size()>numSubCutoffSampsNeeded)
                    bestAboveCutoff.remove();
                
                for(int q=0; q<numSubCutoffSampsNeeded; q++){
                    //add the best samples from the queue
                    //bCutoffs put at the true value, meaning these samples get
                    //the regular least-squares penalties for fitting purposes
                    ConfPair sv = bestAboveCutoff.poll();

                    sample = sv.conf;
                    trueVal = sv.energy[0];
                    
                    ArrayList<Integer> sampTuples = calcSampleTuples(sample);

                    //sampleTuples.add(sampTuples);
                    samples.add(sample);
                    int sampleIndex = samples.size()-1;

                    //update tupleSamples
                    for(int term : sampTuples){
                        //tupleSamples.get(term).add(sampleIndex);
                        tupleNumSamples.set( term, tupleNumSamples.get(term)+1 );
                        //these don't count as tupleSamplesAboveCutoff because they aren't fit as such
                    }

                    trueVals.add(trueVal);
                    curFitVals.add(0.);//to be replaced by updateFitVal

                    sampleBCutoffs.add(trueVal);
                }
                
                if(!bestAboveCutoff.isEmpty())//DEBUG!!!
                    throw new RuntimeException("ERROR: bestAboveCutoff too big!");
                
                break;
            }
        }
    }
    
    
    
    
    void updateSamples2(int tup){
        //We just raised created or raised the order of approximation for the given tuple terms
        //we need to create more samples for it
        
        //VERSION ASSUMING INFINITE BCUTOFF
        //AND AVOIDING SAMPLE DUPLICATION
        
        
        //ASSUMING LACK OF SAMPLE DUPLICATION
        //ENFORCING THIS BELOW BY ADDING ONLY DISTINCT ONES
        int numSampsNeeded = te.numSamplesNeeded(tup) - tupleNumSamples.get(tup);//tupleSamples.get(tup).size();//numDistinctSamples(tup);
        //We may need to have <numSamplesNeeded distinct samples because some tuples have a limited number of possible confs;
        //however, we would overfit if we counted a repeated sample from some other tuple with a small conf space
        // as satisfying the requirement for this tuple
        //(this tuple might have a larger conf space, so this could cause overfitting)
        //Alternative approach is to require 10 distinct samples for each tuple
        //(use DFS to confirm that there are really less)
        //But this may be expensive and unnecessary
        
        
        for(int s=0; s<numSampsNeeded; s++){
            int[] sample = drawUnprunedSample(te.tuples.get(tup), true);
            //we know sample won't be null because we checked if tupleFeasible
            
            
            //don't add if not distinct.  Indicates small conf space
            if( ! isNewSampleDistinct(sample/*,tupleSamples.get(tup)*/))
                continue;
            
            double trueVal = te.scoreAssignmentList(sample);
            
            
            ArrayList<Integer> sampTuples = calcSampleTuples(sample);
            sampTuples.trimToSize();

            //sampleTuples.add(sampTuples);
            samples.add(sample);
            int sampleIndex = samples.size()-1;

            //update tupleSamples
            for(int term : sampTuples){
                //tupleSamples.get(term).add(sampleIndex);
                tupleNumSamples.set( term, tupleNumSamples.get(term)+1 );
                if(trueVal>te.bCutoff)
                    //tupleSamplesAboveCutoff.get(term).set(tupleSamples.get(term).size()-1);
                    tupleSamplesAboveCutoff.get(term).set(tupleNumSamples.get(term)-1);
            }


            trueVals.add(trueVal);
            curFitVals.add(0.);//to be replaced by updateFitVal

            sampleBCutoffs.add(te.bCutoff);
        }
    }
    
    
    
    void updateSamplesMultithread(int tup){
        //Multithreaded version of updateSamples2
        //Requires energy terms that don't share objects they modify!
        
        //ASSUMING LACK OF SAMPLE DUPLICATION
        //ENFORCING THIS BELOW BY ADDING ONLY DISTINCT ONES
        int numSampsNeeded = te.numSamplesNeeded(tup) - tupleNumSamples.get(tup);//tupleSamples.get(tup).size();//numDistinctSamples(tup);
        //We may need to have <numSamplesNeeded distinct samples because some tuples have a limited number of possible confs;
        //however, we would overfit if we counted a repeated sample from some other tuple with a small conf space
        // as satisfying the requirement for this tuple
        //(this tuple might have a larger conf space, so this could cause overfitting)
        //Alternative approach is to require 10 distinct samples for each tuple
        //(use DFS to confirm that there are really less)
        //But this may be expensive and unnecessary
        
        //OK let's draw some samples and then score them in parallel
        ArrayList<int[]> sampleList = new ArrayList<>();
        
        for(int s=0; s<numSampsNeeded; s++){
            int[] sample = drawUnprunedSample(te.tuples.get(tup), true);
            //we know sample won't be null because we checked if tupleFeasible
            
            //don't add if not distinct.  Indicates small conf space
            if( isNewSampleDistinct(sample/*,tupleSamples.get(tup)*/))
                sampleList.add(sample);
        }
        
        double trueValList[] = scoreSamplesMultithread(sampleList);
        
        
        //now record the samples we created
        for(int s=0; s<sampleList.size(); s++){ 
            int sample[] = sampleList.get(s);
            double trueVal = trueValList[s];
            
            ArrayList<Integer> sampTuples = calcSampleTuples(sample);
            sampTuples.trimToSize();

            //sampleTuples.add(sampTuples);
            samples.add(sample);
            //int sampleIndex = samples.size()-1;

            //update tupleSamples
            for(int term : sampTuples){
                //tupleSamples.get(term).add(sampleIndex);
                tupleNumSamples.set( term, tupleNumSamples.get(term)+1 );
                if(trueVal>te.bCutoff)
                    //tupleSamplesAboveCutoff.get(term).set(tupleSamples.get(term).size()-1);
                    tupleSamplesAboveCutoff.get(term).set(tupleNumSamples.get(term)-1);
            }


            trueVals.add(trueVal);
            curFitVals.add(0.);//to be replaced by updateFitVal

            sampleBCutoffs.add(te.bCutoff);
        }
    }
    
    
    double[] scoreSamplesMultithread(final ArrayList<int[]> sampleList){
        //Given a list of samples, compute their scores in different threads
        
        int numSamp = sampleList.size();
        final double trueValList[] = new double[numSamp];
        
        Thread[] threads = new Thread[numSamp];
            
        for(int s=0; s<numSamp; s++){

            final int sampCount = s;
            threads[sampCount] = new Thread(){

                @Override
                public void run(){
                    double trueVal = te.scoreAssignmentList(sampleList.get(sampCount));

                    synchronized(trueValList){
                        trueValList[sampCount] = trueVal;
                    }
                }
            };
        }

        for(int s=0; s<numSamp; s++)//start all the threads
            threads[s].start();

        try {
            for(int s=0; s<numSamp; s++)//wait for all the threads to finish
                threads[s].join();
        }
        catch(InterruptedException e){
            throw new RuntimeException("ERROR: Threaded calc interrupted!");
        }

        return trueValList;
    }
    
    
    boolean isNewSampleDistinct(int[] sample/*, ArrayList<Integer> otherSampList*/){
        //Is the new sample distinct from those listed?  (i.e. those that
        //share a tuple and thus could be the same)
       
        //for(int sampNum2 : otherSampList){//not storing tupleSamples now, so just check them all
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
    
    
    /*int numDistinctSamples(int tup){
        //how many distinct samples are there containing tuple # tup?
        int numDistinct = 0;
        ArrayList<Integer> tupSamp = tupleSamples.get(tup);
        
        for(int s=0; s<tupSamp.size(); s++){
            
            int sample[] = samples.get(tupSamp.get(s));
            //compare to each of the previous samples containing tup, see if it's distinct from them
            boolean sampleDistinct = true;
            
            for(int s2=0; s2<s; s2++){
                int sample2[] = samples.get(tupSamp.get(s2));
                
                boolean samplesIdentical = true;
                
                for(int pos=0; pos<te.numPos; pos++){
                    if(sample[pos]!=sample2[pos]){
                        samplesIdentical = false;
                        break;
                    }
                }
                
                if(samplesIdentical){
                    sampleDistinct = false;
                    break;
                }
            }
            
            if(sampleDistinct)
                numDistinct++;
        }
        
        return numDistinct;
    }*/
    
    
        
    boolean tupleFeasible(RCTuple tup){
        //see if this tup yields sub-bcutoff samples in normal drawing
        //meant to predict (conservatively) if updateSamples will be OK getting samples for it
        
        
        if(te.bCutoff == Double.POSITIVE_INFINITY){
            //samples not restricted to be below some bcutoff,
            //so we can just do DFS to see if there is any unpruned conf containing tup
            int sample[] = drawUnprunedSample(tup, false);
            return (sample!=null);//were we able to draw a sample?
        }
        
        int goodSampCount = 0;//we'll demand this to be a bit higher than necessary in updateSamples,
        //for maxDiscrReps (max number of tries in updateSamples per sample needed)
        
        for(int s=0; s<maxDiscrReps; s++){
                        

            int sample[] = drawUnprunedSample(tup, false);
            boolean success = (sample!=null);//we were able to draw a sample 
            
            //DEBUG!!  Now with DFS drawing can break if no success, because won't ever have success
            if(!success)
                break;
                    
                    
            /*
            int[] sample = new int[te.numPos];

            //if term isn't pruned there must be some way to finishRCSample
            boolean success;
            
            //DEBUG!!
            //do{
                Arrays.fill(sample, -1);

                //fill the current tuple into RCSample
                //for(int[] op : te.tuples.get(tup))
                //    sample[op[0]] = op[1];
                te.assignTupleInSample(sample, tup);
                
                success = finishSample(sample);
            //}
            //while(!success);
            */
            
            if(success){

                //figure out which terms the sample applies to
                double trueVal = te.scoreAssignmentList(sample);

                if(trueVal<=te.bCutoff)
                    goodSampCount++;

                if(goodSampCount==3)//can be pretty sure of good updateSamples
                    return true;
            }
        }
        
        return false;//if we get here the drawing failed!
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
                
                if(te.isPruned(new RCTuple(nextPos, a)))
                    distr[a] = 0;
                
                //look for pruned pairs between a and assignments already in sample
                for(int pos=0; pos<te.numPos; pos++){
                    if(sample[pos]!=-1){
                        
                        if(isPairPrunedInSample(pos, sample[pos], nextPos, a, sample)){
                            distr[a] = 0;
                            break;
                        }
                        /*
                        RCTuple pair = new RCTuple(nextPos, a, pos, sample[pos]);
                        
                        
                        if(te.isPruned(pair)){
                            distr[a] = 0;
                            break;
                        }
                        
                        //and any higher-order pruned stuff...
                        for(RCTuple prunedTup : te.higherOrderPrunedTuples(pair)){
                            if(te.sampleMatchesTuple(sample,prunedTup)){
                                distr[a] = 0;
                                break;
                            }
                        }*/
                    }
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
    void updateFitVals(){
        for(int s=0; s<samples.size(); s++)
            //curFitVals.set( s, te.fitValueForTuples(sampleTuples.get(s)) );
            curFitVals.set( s, te.fitValueForTuples(calcSampleTuples(samples.get(s))) );
        
        //and update residuals accordingly
        updateAllResids();
    }
    
    
    
    
    void updateAllResids(){
        
        totalResid = 0;
        for(int t=0; t<te.tuples.size(); t++)
            tupleResids.set(t,0.);
        
        sampleResids = new ArrayList<>();
        
        for(int s=0; s<samples.size(); s++){
            double targetVal = trueVals.get(s);
            double sampBCutoff = sampleBCutoffs.get(s);
            double sampBCutoff2 = sampBCutoff + te.bCutoff2 - te.bCutoff;
            
            if(targetVal>sampBCutoff){
                if(targetVal>sampBCutoff2){
                    if(curFitVals.get(s)>sampBCutoff){
                        sampleResids.add(0.);
                        continue;//sample doesn't count towards residual
                    }
                    else
                        targetVal = sampBCutoff;
                }
                else {
                    if(curFitVals.get(s)<sampBCutoff)
                        targetVal = sampBCutoff;
                    else if(curFitVals.get(s)<=targetVal){
                        sampleResids.add(0.);
                        continue;
                    }
                }
            }
            
            double sampResid = (curFitVals.get(s)-targetVal)*(curFitVals.get(s)-targetVal);
            totalResid += sampResid;
            //for(int tup : sampleTuples.get(s))
            for(int tup : calcSampleTuples(samples.get(s)) )
                tupleResids.set( tup, tupleResids.get(tup)+sampResid );
            
            sampleResids.add(sampResid);
        }
        
        totalResid /= samples.size();
        worstResid = 0;
        for(int t=0; t<te.tuples.size(); t++){
            tupleResids.set(t,tupleResids.get(t)/tupleNumSamples.get(t));//tupleSamples.get(t).size());
            worstResid = Math.max(worstResid,tupleResids.get(t));
        }
        
        
        
        /*System.out.print("Current worst residual: "+worstResid+" for");
        for(int ind : worstResidIndex)
            System.out.print(" "+nwpf.termString(ind));
        System.out.println();*/
        
        //System.out.println("Current biggest term residual gradient norm: "+biggestGradNorm
        //        +" for "+biggestGradNormIndex);
    }
    
    
    
    void addTuple(int tup){
        //add a tuple that has just been put into TupleExpander, updating samples accordingly
        //tup will be tuples.size()-1
        
        //record relationship to current samples
        //tupleSamples.add( new ArrayList<Integer>() );
        tupleNumSamples.add(0);
        
        tupleSamplesAboveCutoff.add(new BitSet());
        
        for(int s=0; s<samples.size(); s++){
            int sample[] = samples.get(s);
            
            if(te.sampleMatchesTuple(sample,te.tuples.get(tup))){
                //sampleTuples.get(s).add(tup);
                
                //tupleSamples.get(tup).add(s);
                tupleNumSamples.set( tup, tupleNumSamples.get(tup)+1 );
                
                if(trueVals.get(s) > sampleBCutoffs.get(s))
                    //tupleSamplesAboveCutoff.get(tup).set(tupleSamples.get(tup).size()-1);
                    tupleSamplesAboveCutoff.get(tup).set(tupleNumSamples.get(tup)-1);
            }
        }
        
        
        //preallocate residuals
        tupleResids.add(0.);
        
        //draw any new samples needed
        updateSamples(tup);
    }
    
    
    void updateSampBelowCutoff(){
        //accommodate raising of bcutoff by updating the tupleSamplesAboveCutoff bitsets
        
        if(te.bCutoff == Double.POSITIVE_INFINITY)//expected case now...
            return;
        
        throw new RuntimeException("ERROR: Not expecting finite bcutoff here...");//would need the stuff below
        
        /*
        for(int term=0; term<te.tuples.size(); term++){
            for(int samp=0; samp<tupleSamples.get(term).size(); samp++){
                if(trueVals.get(tupleSamples.get(term).get(samp))>sampleBCutoffs.get(samp))
                    tupleSamplesAboveCutoff.get(term).set(samp,true);
                else
                    tupleSamplesAboveCutoff.get(term).set(samp,false);
            }
        }*/
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
        //For this we set up the list of RCs available at each position.  
        //Each time we'll search the position with the most eliminated RCs (by pairs+ pruning)
        //this way hopefully can clear out tuples that need pruning relatively quickly
        //random ordering of assignments at each position
        Arrays.fill(sample,-1);
        te.assignTupleInSample(sample, tup);
        
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

        //DEBUG!!  Want to see what's taking so long on second-round 30.5, 40.cont
        long sampStartTime = System.currentTimeMillis();
        
        if(tuplePossible)
            tuplePossible = unprunedSampDFS(sample, allowedRCs, numElim);
        
        
        //DEBUG!!!
        long sampTime = System.currentTimeMillis() - sampStartTime;
        //taking over 10 s is going to be an issue
        if(sampTime > 10000){
            System.out.println();
            System.out.println("Sampling took over 10 s (ms shown): "+sampTime);
            System.out.println("Target tuple: "+tup.stringListing());
            if(tuplePossible)
                System.out.println("Sample: "+new RCTuple(sample).stringListing());
            else
                System.out.println("No sample found.");
            System.out.println();
        }
        
        
        
        
        if(tuplePossible)
            return sample;
        else {
            if(errorIfImpossible)
                throw new RuntimeException("ERROR: No unpruned samples available for tuple "+tup.stringListing());
            else
                return null;
        }
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
                            if(isPairPrunedInSample(nextPos,opt,uaOpt[0],uaOpt[1],sample)){
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
                        for(int pos2=0; pos2<te.numPos; pos2++){
                            if(sample[pos2]!=-1){//pos2 assigned
                                if(isPairPrunedInSample(pos2,sample[pos2],pos,rc,sample)){
                                    incompatible = true;
                                    break;
                                }
                            }
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
    
    
    ArrayList<Integer> shuffleOptions(ArrayList<int[]> optList){
        //given a bunch of (pos,rc) pairs, return the rc's in random order
        ArrayList<Integer> ans = new ArrayList<>();
        int numRCs = optList.size();
        
        ArrayList<Integer> ordering = randomOrdering(numRCs);
        
        for(int a=0; a<numRCs; a++){
            ans.add( optList.get(ordering.get(a))[1] );
        }
        
        return ans;
    }
    
    ArrayList<Integer> randomOrdering(int n){
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
