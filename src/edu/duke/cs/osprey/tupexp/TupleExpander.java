/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tupexp;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;
import java.util.TreeSet;

/**
 *
 * @author mhall44
 */
/**
 *
 * Create an energy matrix that approximates a function of many integers (e.g., list of RCs for all
 * flexible residues) as a sum of terms dependent on low-order tuples (e.g., tuple
 * of RCs at a few positions)
 * 
 * 
 * @author mhall44
 */

public abstract class TupleExpander implements Serializable {
    
    int numPos;//number of positions
    int numAllowed[];//number of possible assignments at each position
    
    ArrayList<RCTuple> tuples = new ArrayList<>();//the tuples we're trying to expand in
    
    double[] tupleTerms;//terms for the tuples
    double constTerm = Double.NaN;
    //so we approximate our function as constTerm + sum_{tup in current tuples} tupleTerms[tup]
    
    TESampleSet trainingSamples=null, CVSamples=null;
    
    FittingObjFcn fof = new FittingObjFcn();//the objective function for fitting (e.g. basic or modified least squares)
    
    
    double pruningInterval;//what pruning interval is this expansion valid up to?
       
    
    static int printedUpdateNumTuples = 5;//100 w/o PB.  We print updates as we add tuples...this is how often

    boolean canCheckPartialPruning = true;//can evaluate new tuples
    //by using isPairPrunedInSample rather than checking isPruned on whole conf
    //PLUG can break this

    public TupleExpander (int numPos, int[] numAllowed/*, double constTerm*/, double pruningInterval, LUTESettings luteSettings) {
        this.numPos = numPos;
        this.numAllowed = numAllowed;
        //this.constTerm = constTerm;
        this.pruningInterval = pruningInterval;
        
        for(int pos=0; pos<numPos; pos++){
            /*tuplesForAssignments.add(new ArrayList<ArrayList<Integer>>());
            for(int a=0; a<numAllowed[pos]; a++)
                tuplesForAssignments.get(pos).add(new ArrayList<Integer>());*/
            
            assignmentSets.add(new ArrayList<ArrayList<Integer>>());
        }
        
        fof = new FittingObjFcn(pruningInterval,0.5,luteSettings.useRelWt,luteSettings.useThreshWt);
    }
    
    
    
    
    
    double computeInitGMECEst(){//let's find this by some random iterations...
        double ans = Double.POSITIVE_INFINITY;
        TESampleSet tss = new TESampleSet(this);
        
        System.out.println("Computing initial GMEC estimate...");
        
        
        //Let's make sure, by DFS, that there is at least one unpruned conf
        //If not then our init GMEC is infinity
        int testSamp[] = new int[numPos];
        Arrays.fill(testSamp, -1);
        if( tss.finishSampleDFS(testSamp) == null )
            return Double.POSITIVE_INFINITY;
        
        
        for(int iter=0; iter<50/*00*/; iter++){
            int sample[] = new int[numPos];
            boolean success;
            do {
                Arrays.fill(sample,-1);
                success = tss.finishSample(sample);
            } while (!success);
            
            double score = scoreAssignmentList(sample);
            ans = Math.min(ans,score);
        }

        System.out.println("Initial GMEC estimate: "+ans);
        
        return ans;
    }
    
    
    public double calcExpansion(ArrayList<RCTuple> tuplesToFit){
        //calculate the tuple coefficients
        //return cross-validation total residual        
        
        if(Double.isNaN(constTerm))//constTerm not computed yet
            constTerm = computeInitGMECEst();
        
        if( constTerm == Double.POSITIVE_INFINITY ){
            System.out.println("No conformations found for tuple expansion.  ");
            tupleTerms = new double[0];//no coefficients needed besides constTerm
            return 0;//all confs are pruned, and will be correctly assigned infinite energy
        }
        
        //this can be used for the first tuple expansion for this object, or to add tuples later (not take away though)
        
        System.out.println("About to calculate tuple expansion with "+tuplesToFit.size()+" tuples.");
        //System.out.println("constTerm: "+constTerm+" bCutoff: "+bCutoff+" bCutoff2: "+bCutoff2);

        
        setupSamples(tuplesToFit);//set up the training set (in the process, prune tuples that don't provide reasonable energies)
                
        fitLeastSquares();
        
        trainingSamples.updateFitVals(fof);
        CVSamples.updateFitVals(fof);
        
        System.out.println("TRAINING SAMPLES: ");
        trainingSamples.printResids();
        
        System.out.println("CV SAMPLES: ");
        CVSamples.printResids();
        
        return CVSamples.totalResid;
    }
    
    
    
    private void checkHighEnergyConfs(){
        //This is a function for debugging inadequate fit quality
        //Checks what kinds of high-energy confs each tuple has in the training set
        //Then, tries to redo the fitting with problem tuples removed--
        //first, tuples whose confs are all high-energy, then who have some high-energy confs
        //(high-energy here defined as lowest-energy sample + 100)
        
        System.out.println("Checking high-energy confs");
        
        double lowestSampleE = Collections.min( trainingSamples.trueVals );
        
        System.out.println("Lowest sample E: "+lowestSampleE);
        
        //figure out the lowest sample energy for each tuple
        double tupleBestE[] = new double[tuples.size()];
        Arrays.fill(tupleBestE, Double.POSITIVE_INFINITY);
        
        //also lowest LB and lowest contTerm.  Since pruning can be based on either of these
        //in different combinations.  
        //these here go with the bestE
        double tupleBestELB[] = new double[tuples.size()];
        double tupleBestEContE[] = new double[tuples.size()];
        
        //and then, for pruning purposes, I'd like to know the best LB and cont E overerall
        double tupleBestLB[] = new double[tuples.size()];
        double tupleBestContE[] = new double[tuples.size()];
        Arrays.fill(tupleBestLB, Double.POSITIVE_INFINITY);
        Arrays.fill(tupleBestContE, Double.POSITIVE_INFINITY);
        
        
        
        
        for(int s=0; s<trainingSamples.samples.size(); s++){
            double E = trainingSamples.trueVals.get(s);
            ArrayList<Integer> sampTups = trainingSamples.calcSampleTuples(trainingSamples.samples.get(s));
            
            double LB = ((ConfETupleExpander)this).sp.lowerBound(trainingSamples.samples.get(s));
            double contE = E - LB;
            
            for(int t : sampTups){
                if(E < tupleBestE[t]){
                    tupleBestE[t] = E;
                    tupleBestELB[t] = LB;
                    tupleBestEContE[t] = contE;
                }
                tupleBestLB[t] = Math.min(tupleBestLB[t], LB);
                tupleBestContE[t] = Math.min(tupleBestContE[t], contE);
            }
        }
        
        TreeSet<Integer> badTuples100 = new TreeSet<>();//tuples whose samples are all > lowestSampleE + 100
        TreeSet<Integer> badTuples50 = new TreeSet<>();//tuples whose samples are all

        
        
        //Now print these lowest energies and also figure out which tuples we might prune on this basis
        //This will tell us about what tuples might really be "bad"
        //(i.e., have no confs w/ reasonable energy)
        System.out.println("Tuple Best E's (then bestE LB and cont E; best LB, contE): ");
        for(int t=0; t<tuples.size(); t++){
            System.out.println(tupleBestE[t] + " " + 
                    tupleBestELB[t] + " " + tupleBestEContE[t] + " " + 
                    tupleBestLB[t] + " " + tupleBestContE[t] + " " + tuples.get(t).stringListing());
           
            if(tupleBestE[t] > lowestSampleE + 50)
                badTuples50.add(t);
            if(tupleBestE[t] > lowestSampleE + 100)
                badTuples100.add(t);
        }
        System.out.println("End Tuple Best E's");
        
        
        
        //OK now see if there are confs over 100 that have no "bad" tuples (in the 50 sense)
        //for each of these, the tuple with the worst tupleBestE is suspect
        TreeSet<Integer> suspectTuples = new TreeSet<>();
        suspectTuples.addAll(badTuples50);
        
        System.out.println("Bad samples w/o bad tuples (E, LB, worstTupE, worstTup): ");
        for(int s=0; s<trainingSamples.samples.size(); s++){
            if(trainingSamples.trueVals.get(s) > lowestSampleE+100){
                ArrayList<Integer> sampTups = trainingSamples.calcSampleTuples(trainingSamples.samples.get(s));
                
                double worstTupE = Double.NEGATIVE_INFINITY;
                int worstTup = -1;
                for(int t : sampTups){
                    if(tupleBestE[t]>worstTupE){
                        worstTup = t;
                        worstTupE = tupleBestE[t];
                        if(worstTupE > lowestSampleE + 50)
                            break;
                    }
                }
                
                if(worstTupE <= lowestSampleE + 50){//seriously bad conf E but tuple ok
                    //hmm...
                    double E = trainingSamples.trueVals.get(s);
                    double LB = ((ConfETupleExpander)this).sp.lowerBound(trainingSamples.samples.get(s));
                    System.out.println(E+" "+LB+" "+worstTupE+" "+tuples.get(worstTup).stringListing());
                    suspectTuples.add(worstTup);
                }
            }
        }
        System.out.println("End Bad samples w/o bad tuples");
        

        
        //OK the bad tuple sets will be wrong indices once we start deleting tuples
        //But we can fix this by converting to lists of tuples
        ArrayList<RCTuple> badList100 = new ArrayList<>();
        for(int t : badTuples100)
            badList100.add(tuples.get(t));
        
        ArrayList<RCTuple> badList50 = new ArrayList<>();
        for(int t : badTuples50)
            badList50.add(tuples.get(t));
        
        ArrayList<RCTuple> suspectList = new ArrayList<>();
        for(int t : suspectTuples)
            suspectList.add(tuples.get(t));
        
        
        
        System.out.println("REDOING FITTING W/O BAD TUPLES100");
        redoFittingWithoutTuples(badList100);
        System.out.println("REDOING FITTING W/O BAD TUPLES50");
        redoFittingWithoutTuples(badList50);
        System.out.println("REDOING FITTING W/O SUSPECT TUPLES");
        redoFittingWithoutTuples(suspectList);
        System.out.println("HIGH CONF-E CHECK DONE");
    }
    
    private void redoFittingWithoutTuples(ArrayList<RCTuple> badTuples){
        //OK now redo the fitting 
        trainingSamples = null;
        CVSamples = null;
        
        for(RCTuple t : badTuples){
            tuples.remove(t);
            pruneTuple(t);
        }
        
        ArrayList<RCTuple> tuplesToFit = tuples;
        tuples = new ArrayList<>();
        numSampsPerTuple = 10;
        
        setupSamples(tuplesToFit);//set up the training set (in the process, prune tuples that don't provide reasonable energies)
                
        fitLeastSquares();
        
        trainingSamples.updateFitVals(fof);
        CVSamples.updateFitVals(fof);
        
        System.out.println("TRAINING SAMPLES: ");
        trainingSamples.printResids();
        
        System.out.println("CV SAMPLES: ");
        CVSamples.printResids();
    }
    
    
    
    
    
    
    
    
    //Counts how many RC tuples were pruned by inability of DFS to find samples for them
    //(i.e., it follows from the existing pruning matrix that they can be pruned)
    int ezPruneCount = 0;
    
    void setupSamples(ArrayList<RCTuple> tuplesToFit){
        //generate the training sample set, and for each tuple we want to fit, either prune it or get samples for it
        
        numSampsPerTuple = 10;
                
        if(trainingSamples==null)//first tuple expansion for this object: initialize training samples
            trainingSamples = new TESampleSet(this);
        else {
            //figure out which tuples are already included in this expansion, so we don't try to add them again
            //note: assuming any such redundant tuples have the same ordering of positions!
            for(int tupNum=tuplesToFit.size()-1; tupNum>=0; tupNum--){
                RCTuple tup = tuplesToFit.get(tupNum);
                boolean removeTuple = false;
                
                for(RCTuple tupHere : tuples){
                    if(tupHere.isSameTuple(tup)){
                        tuplesToFit.remove(tupNum);
                        removeTuple = true;
                        break;
                    }
                }
                
                if(!removeTuple) {
                    //also want to remove any redundant tuples among the new ones...  (same assumption about ordering)
                    for(int tupNum2=0; tupNum2<tupNum; tupNum2++){
                        if( tup.isSameTuple( tuplesToFit.get(tupNum2) ) ){
                            tuplesToFit.remove(tupNum);
                            break;
                        }
                    }
                }
            }
        }
        
        
        System.out.println("Adding tuples to expansion and drawing training samples...");
        
        for(int tupNum=0; tupNum<tuplesToFit.size(); tupNum++){
            if(tupNum>0 && (tupNum%printedUpdateNumTuples==0))
                System.out.println(tupNum+" tuples added");
  
            tryAddingTuple(tuplesToFit.get(tupNum));//prune or get samples
        }
        
        
        System.out.println("EZPRUNE COUNT: "+ezPruneCount);
        
        System.out.println("Updating samples to finish training set...");
        
        for(int t=0; t<tuples.size(); t++)
            //replace any samples that were removed during pruning
            trainingSamples.updateSamples(t);
        
        
        //let's try to get about >=2x as many samples as tuples, to avoid overfitting
        //increase number of samples per tuple to achieve this
        if(trainingSamples.samples.size() < 2*tuples.size()){
            numSampsPerTuple *= 2*tuples.size()/trainingSamples.samples.size()+1;//+1 to round up
            
            for(int t=0; t<tuples.size(); t++)//get the additional samples we need
                trainingSamples.updateSamples(t);
        }
        
        System.out.println("Training set done.");
        System.out.println("Drawing CV samples.");
        
        //now make sure the CV samples are updated too
        if(CVSamples==null){
            CVSamples = new TESampleSet(this);
            for(int t=0; t<tuples.size(); t++)
                CVSamples.updateSamples(t);
        }
            
        System.out.println("CV set done.");
    }
    
    
    
    
    
    void fitLeastSquares(){
        //set up to call fitSeriesIterative
        int numTrainingSamples = trainingSamples.samples.size();
        
        double[] trueVals = new double[numTrainingSamples];
        
        for(int s=0; s<numTrainingSamples; s++){
            trueVals[s] = trainingSamples.trueVals.get(s) - constTerm;
        }
        
        //weights, bcutoffs (even bcutoffs2) may be needed but hopefully not (if so modify CG: full iterative quite expensive)
        /*double fitTerms[] = SeriesFitter.fitSeriesIterative(samp, trueVals, weights, lambda, false, 1,
                    bCutoffs, bCutoffs2, 1, null);*/
        
        TupleIndexMatrix tim = getTupleIndexMatrix();
        CGTupleFitter fitter = fof.makeTupleFitter(tim, trainingSamples.samples, tuples.size(), trueVals);
        //CGTupleFitter fitter = new CGTupleFitter(tim, trainingSamples.samples, tuples.size(), trueVals);
        
        double fitTerms[] = fitter.doFit();
        tupleTerms = fitTerms;
    }
 
    
    
    int numSampsPerTuple = 10;
    
    int numSamplesNeeded(int tup){
        //for a tuple (index in tuples), how many samples are needed?
        return numSampsPerTuple;//one param per tuple, so 10 samples should securely avoid overfitting
    }

    public void tryAddingTuple(RCTuple tup){
        boolean tupFeas = trainingSamples.tupleFeasible(tup);
        if(!tupFeas){
            //make sure we didn't miss a tuple we already know from previous samples
            //to be feasible.  This could result in pruning other tuples and general mess
            for(TESampleSet tss : new TESampleSet[]{trainingSamples,CVSamples}){
                if(tss!=null){
                    for(int[] sample : tss.samples){
                        if(sampleMatchesTuple(sample,tup)){
                            tupFeas = true;
                            break;
                        }
                    }
                }
            }
        }

        if(tupFeas){
            tuples.add(tup);
            int newTupleIndex = tuples.size()-1;//tup index in tuples

            trainingSamples.addTuple(newTupleIndex);
            if(CVSamples!=null)
                CVSamples.addTuple(newTupleIndex);
        }
        else {
            //mark as pruned
            pruneTuple(tup);
            ezPruneCount++;
        }
    }
       
    
    double fitValueForTuples(ArrayList<Integer> tuples){
        //given a list of tuples to which a term belongs
        //return the fit value
        double ans = constTerm;
        for(int tup : tuples)
            ans += tupleTerms[tup];
        return ans;
    }
    
    
    
    
    
    //HANDLING OF ASSIGNMENT SETS
    //These will be specified by negative values in a tuple array
    //Value -i for residue res denotes the set of assignments indexed by i in assignmentSets.get(res)
    ArrayList<ArrayList<ArrayList<Integer>>> assignmentSets = new ArrayList<>();
    
    
    int getAssignmentSet(int res, ArrayList<Integer> assignmentList){
        Collections.sort(assignmentList);
        
        ArrayList<ArrayList<Integer>> resASets = assignmentSets.get(res);
        
        for(int i=0; i<resASets.size(); i++){
            if(resASets.get(i).size() == assignmentList.size()){
                boolean listsEqual = true;
                for(int j=0; j<assignmentList.size(); j++){
                    if(resASets.get(i).get(j) != assignmentList.get(j))
                        listsEqual = false;
                }
            
                if(listsEqual)
                    return i;//found a match for this assignment set
            }
        }
        
        //if we get here, the assignment set is not currently listed, so we must do so
        resASets.add(assignmentList);
        return resASets.size() - 1;
    }
    
    
    void assignTupleInSample(int sample[], RCTuple tuple){
        //assign the sample to have the assignments specified by tuple
        //if there are assignment sets, pick randomly, though avoid pruned pairs
        for(int posCount=0; posCount<tuple.pos.size(); posCount++){
            
            int pos = tuple.pos.get(posCount);
            int rc = tuple.RCs.get(posCount);
            
            if(rc>=0)//specific assignments at position op[0]
                sample[pos] = rc;
            else{//assignment set...draw randomly
                ArrayList<Integer> aSet = assignmentSets.get(pos).get(-rc);
                int assignment = aSet.get( new Random().nextInt(aSet.size()) );
                
                while(!checkAssignmentUnpruned(sample,pos,assignment)){
                    //looks like this assignment is incompatible with rest of sample assigned so far
                    //redraw from rest of possible assignments...
                    ArrayList<Integer> aSetRed = new ArrayList<>();
                    for(int a : aSet){
                        if(a!=assignment)
                            aSetRed.add(a);
                    }
                    
                    aSet = aSetRed;
                    
                    if(aSet.isEmpty())
                        throw new RuntimeException("ERROR: Can't find compatible assignment for tuple in sample...");
                        //let's currently address this as an error...
                    
                    assignment = aSet.get( new Random().nextInt(aSet.size()) );
                }
                
                sample[pos] = assignment;
            }
        }
    }
    
    
    boolean checkAssignmentUnpruned(int sample[], int pos, int assignment){
        //given a partially defined sample (-1s for unassigned positions)
        //see if the given assignment at position pos is compatible pairs pruning-wise
        for(int pos2=0; pos2<sample.length; pos2++){
            
            RCTuple pair = new RCTuple(pos,assignment,pos2,sample[pos2]);
            
            if(sample[pos2]!=-1){
                if(isPruned(pair))
                    return false;
            }
        
            for(RCTuple prunedTup : higherOrderPrunedTuples(pair)){
                if(sampleMatchesTuple(sample,prunedTup))
                    return false;
            }
        }
                    
        //if we get here, we're unpruned
        return true;
    }
    
    
    
    boolean sampleMatchesTuple(int sample[], RCTuple tup){
        boolean termApplies = true;

        for(int posNum=0; posNum<tup.pos.size(); posNum++){
            
            int pos = tup.pos.get(posNum);
            int rc = tup.RCs.get(posNum);
            
            if(sample[pos]==-1)//undefined considered not to match
                return false;
            
            if(rc>=0){//sample must match specific assignment in tuple
                if( sample[pos] != rc ){
                    termApplies = false;
                    break;
                }
            }
            else {//sample must be in assignment set in tuple
                ArrayList<Integer> aSet = assignmentSets.get(pos).get(-rc);
                if( ! aSet.contains(sample[pos]) ){
                    termApplies = false;
                    break;
                }
            }
        }

        return termApplies;
    }
    
    
    public EnergyMatrix getEnergyMatrix(){
        //put the tuples into an energy matrix and return it
        EnergyMatrix ans = new EnergyMatrix(numPos,numAllowed,pruningInterval);
        
        ans.setConstTerm(constTerm);
        

        //first, let's put in values not explicitly in the tuple expansion
        //the tuple-expansion energy is the sum of values for all tuples in a conf
        //but an energy matrix expects values for all one-body and pairwise energies
        //so we'll specify these to be 0 for RCs and pairs not in the expansion.
        //HOWEVER, pruned tuples must be marked as impossible, i.e. infinite.  
        for(int pos=0; pos<numPos; pos++){
            for(int rc=0; rc<numAllowed[pos]; rc++){
                
                if(isPruned(new RCTuple(pos,rc)))
                    ans.setOneBody(pos, rc, Double.POSITIVE_INFINITY);
                else
                    ans.setOneBody(pos, rc, 0.);
                
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<numAllowed[pos2]; rc2++){
                        RCTuple pair = new RCTuple(pos,rc,pos2,rc2);
                        if(isPruned(pair))
                            ans.setPairwise(pos, rc, pos2, rc2, Double.POSITIVE_INFINITY);
                        else
                            ans.setPairwise(pos, rc, pos2, rc2, 0.);
                        
                        //higher tuples 0 by default...just mark pruned ones
                        for(RCTuple prunedTup : higherOrderPrunedTuples(pair))
                            ans.setTupleValue(prunedTup, Double.POSITIVE_INFINITY);
                    }
                }
                
            }
        }
        
        
        for(int tupNum=0; tupNum<tuples.size(); tupNum++){
            ans.setTupleValue( tuples.get(tupNum), tupleTerms[tupNum] );
        }
        
        return ans;
    }
    
    
    
    public TupleIndexMatrix getTupleIndexMatrix(){
        //Make a matrix of the indices in tuples of each tuple
        TupleIndexMatrix ans = new TupleIndexMatrix(numPos, numAllowed, pruningInterval);
        
        //first, fill in -1's for one-body and tuple terms, since values are expected at these matrix entries
        //-1 indicates absence of a tuple in the expansion, and will be overwritten as needed
        for(int pos=0; pos<numPos; pos++){
            for(int rc=0; rc<numAllowed[pos]; rc++){
                
                ans.setOneBody(pos, rc, -1);
                
                for(int pos2=0; pos2<pos; pos2++){
                    for(int rc2=0; rc2<numAllowed[pos2]; rc2++){
                        ans.setPairwise(pos, rc, pos2, rc2, -1);
                    }
                }
            }
        }
        
        
        for(int tupNum=0; tupNum<tuples.size(); tupNum++)
            ans.setTupleValue( tuples.get(tupNum), tupNum );
        
        return ans;
    }
    
    
    //functions dependent on what the assignments mean, etc.
    
    //score a list of assignments for each position
    abstract double scoreAssignmentList(int[] assignmentList);
    
    //prune, or check pruning of, a tuple
    abstract boolean isPruned(RCTuple tup);
    abstract void pruneTuple(RCTuple tup);
    
    
    //list higher-order pruned tuples that include the specified pair
    abstract ArrayList<RCTuple> higherOrderPrunedTuples(RCTuple tup);
    
    

    

    
}
