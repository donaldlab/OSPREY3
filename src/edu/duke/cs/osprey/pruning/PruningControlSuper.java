/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.pruning;

import edu.duke.cs.osprey.confspace.SearchProblemSuper;

/**
 *
 * @author hmn5
 */
//Implements PruningControl with superRCs
public class PruningControlSuper {

    SearchProblemSuper searchSpace;

    double pruningInterval;//relative energy threshold for pruning (includes I for iMinDEE)
    boolean typeDep;//type-dependent pruning
    double boundsThresh;//absolute threshold (for Bounds pruning)
    int algOption;//1-4, for different levels of pruning
    boolean useFlags;//Use pair pruning
    boolean useTriples;//Use triple pruning
    boolean preDACS;//do pruning appropriate for all-conf-space pruning before DACS
    boolean useEPIC;//Use EPIC in pruning (to get lower bound of continuous E for pruning candidate)
    boolean useTupExp;//prune based on tup-exp energy matrix

    double stericThresh;//Steric pruning threshold
    boolean onlyGoldstein = false;//for competitor pruning

    public PruningControlSuper(SearchProblemSuper searchSpace, double pruningInterval, boolean typeDep,
            double boundsThresh, int algOption, boolean useFlags, boolean useTriples,
            boolean preDACS, boolean useEPIC, boolean useTupExp, double stericThresh) {
        this.searchSpace = searchSpace;
        this.pruningInterval = pruningInterval;
        this.typeDep = typeDep;
        this.boundsThresh = boundsThresh;
        this.algOption = algOption;
        this.useFlags = useFlags;
        this.useTriples = useTriples;
        this.preDACS = preDACS;
        this.useEPIC = useEPIC;
        this.useTupExp = useTupExp;
        this.stericThresh = stericThresh;
    }

    public void prune() {

        System.out.println();
        System.out.println("BEGINNING PRUNING.  PRUNING INTERVAL: " + pruningInterval);
        System.out.println();

        long startTime = System.currentTimeMillis();

        PrunerSuper dee = new PrunerSuper(searchSpace, typeDep, boundsThresh, pruningInterval, useEPIC, useTupExp);

        //now go through the various types of pruning that we support
        //see KSParser
        //possibly start with steric pruning?  
        if (Double.isFinite(stericThresh) && !onlyGoldstein) {
            dee.pruneSteric(stericThresh);
        }

        boolean done = false;

        //numbers pruned so far
        int numPrunedRot = searchSpace.pruneMat.countPrunedRCs();
        int numPrunedPairs = 0;
        if ((useFlags) || (algOption >= 3)) //pairs pruning is performed
        {
            numPrunedPairs = searchSpace.pruneMat.countPrunedPairs();
        }

        for (int numRuns = 0; !done; numRuns++) { //repeat the pruning cycle until no more rotamers are pruned	

            System.out.println("Starting DEE cycle run: " + numRuns);

            //		if (doMinimize && !localUseMinDEEPruningEw) //precompute the interval terms in the MinDEE criterion
            //			rs.doCompMinDEEIntervals(mp.numberMutable, mp.strandMut, prunedRotAtRes, scaleInt, maxIntScale);
            //Depending on the chosen algorithm option, apply the corresponding pruning criteria;			
            dee.prune("GOLDSTEIN");

            /*
             if ((algOption>=3)) //simple Goldstein pairs
             dee.prune("GOLDSTEIN PAIRS MB");

             if ((useFlags)||(algOption>=3))
             dee.prune("BOUNDING FLAGS");

             dee.prune("CONFSPLIT1");
             //note: conf splitting is equivalent to pruning pairs (say (i_r,j_s)) with overlapping
             //competitors (say (i_t,j_s)), and then seeing what singles are pruned as a result
             //we already do this
            
             dee.prune("BOUNDS");
             */
            //check how many rotamers/pairs are pruned now
            int newNumPrunedRot = searchSpace.pruneMat.countPrunedRCs();
            int newNumPrunedPairs = 0;
            if ((useFlags) || (algOption >= 3)) //pairs pruning is performed
            {
                newNumPrunedPairs = searchSpace.pruneMat.countPrunedPairs();
            }

            if ((newNumPrunedRot == numPrunedRot) && (newNumPrunedPairs == numPrunedPairs) && (!preDACS)) {
                //no more rotamers pruned, so perform the computationally-expensive 2-sp split-Pruner and pairs

                if (algOption >= 3) { //simple Goldstein pairs
                    dee.prune("GOLDSTEIN PAIRS FULL");

                    if (useTriples) {
                        dee.prune("GOLDSTEIN TRIPLES");
                    }
                }

                /*
                 if ((algOption>=2)){ //2-sp conf splitting
                 dee.prune("CONFSPLIT2");
                 }

                    

                 if(algOption >= 4){
                 dee.prune("INDIRECT PAIRS");
                 dee.prune("INDIRECT");
                 }
                 */
                //check if 2-sp split-Pruner and pairs pruned new rotamers
                newNumPrunedRot = searchSpace.pruneMat.countPrunedRCs();
                newNumPrunedPairs = 0;
                if ((useFlags) || (algOption >= 3)) //pairs pruning is performed
                {
                    newNumPrunedPairs = searchSpace.pruneMat.countPrunedPairs();
                }
            }

            int numPrunedRotThisRun = newNumPrunedRot - numPrunedRot;
            int numPrunedPairsThisRun = newNumPrunedPairs - numPrunedPairs;

            System.out.println("Num pruned rot this run: " + numPrunedRotThisRun);
            System.out.println("Num pruned pairs this run: " + numPrunedPairsThisRun);
            System.out.println();

            if (numPrunedRotThisRun == 0 && numPrunedPairsThisRun == 0) {
                done = true;
            }

            numPrunedRot = newNumPrunedRot;
            numPrunedPairs = newNumPrunedPairs;
        }

        long pruneTime = System.currentTimeMillis() - startTime;

        System.out.println("Pruning time: " + pruneTime + " ms");
    }

    public void setOnlyGoldstein(boolean onlyGoldstein) {
        this.onlyGoldstein = onlyGoldstein;
    }

}
