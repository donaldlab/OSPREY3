/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.pruning;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.tools.TimeFormatter;

/**
 *
 * @author mhall44
 */
public class NewPruningControl {
    
    
    public static enum ReportMode {
    	
        Long(true) {
        	
            @Override
            public void preamble(double pruningInterval) {
                System.out.println();
                System.out.println("BEGINNING PRUNING.  PRUNING INTERVAL: " + pruningInterval);
                System.out.println();
            }
            
            @Override
            public void preRun(int run) {
                System.out.println("Starting DEE cycle run: " + run);
            }
            
            @Override
            public void postRun(int run, int numSinglesPruned, int numPairsPruned) {
                System.out.println("Num pruned rot this run: " + numSinglesPruned);
                System.out.println("Num pruned pairs this run: " + numPairsPruned);
                System.out.println();
            }
            
            @Override
            public void conclusion(int numRuns, int numSinglesPruned, int numPairsPruned, long pruneTimeMs) {
                System.out.println("Pruning time: " + pruneTimeMs + " ms" );
                System.out.println();
            }
        },
        
        Short(false) {
        	
            @Override
            public void conclusion(int numRuns, int numSinglesPruned, int numPairsPruned, long pruneTimeMs) {
                System.out.println(String.format("DEE took %s to prune %d singles and %d pairs in %d runs",
                    TimeFormatter.format(pruneTimeMs*TimeFormatter.NSpMS, 1),
                    numSinglesPruned, numPairsPruned, numRuns
                ));
            }
        },
        
        None(false) {
            // nothing to do
        };
    	
    	private boolean isVerbose;
    	
    	private ReportMode(boolean isVerbose) {
    		this.isVerbose = isVerbose;
    	}
    	
    	public boolean isVerbose() {
    		return isVerbose;
    	}
        
    	// do nothing by default on all these methods, override to report things
        public void preamble(double pruningInterval) {}
        public void preRun(int run) {}
        public void postRun(int run, int numSinglesPruned, int numPairsPruned) {}
		public void conclusion(int numRuns, int numSinglesPruned, int numPairsPruned, long pruneTimeMs) {}
    }
    
    PrecomputedMatrices precompMat;
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
    
    private ReportMode reportMode;
    
    public NewPruningControl(PrecomputedMatrices precompMat, double pruningInterval, boolean typeDep, 
            double boundsThresh, int algOption, boolean useFlags, boolean useTriples, 
            boolean preDACS, boolean useEPIC, boolean useTupExp, double stericThresh) {
        this.precompMat = precompMat;
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
        
    	// be verbose by default, since that's the old behavior
        reportMode = ReportMode.Long;
    }
    
    public void setPruningInterval(double val) {
    	pruningInterval = val;
    }
    
    public void setUseEPIC(boolean val) {
    	useEPIC = val;
    }
    
    public void setUseTupExp(boolean val) {
    	useTupExp = val;
    }
    
    public void setReportMode(ReportMode val) {
    	reportMode = val;
        if (reportMode == null) {
            reportMode = ReportMode.None;
        }
    }
    
    public void prune() {
        
        reportMode.preamble(pruningInterval);
        
        long startTime = System.currentTimeMillis();
        
        NewPruner dee = new NewPruner(precompMat,typeDep,boundsThresh,pruningInterval,useEPIC,useTupExp);
        dee.setVerbose(reportMode.isVerbose());
        
        //now go through the various types of pruning that we support
        //see KSParser
        
        //possibly start with steric pruning?  
        if(Double.isFinite(stericThresh) && !onlyGoldstein)
            dee.pruneSteric(stericThresh);
        
        
        boolean done = false;
        
        //numbers pruned so far
        int numPrunedRot = precompMat.getPruneMat().countPrunedRCs();
        int numPrunedPairs = 0;
        if ((useFlags)||(algOption>=3)) //pairs pruning is performed
            numPrunedPairs = precompMat.getPruneMat().countPrunedPairs();
        
        
        int numRuns;
        for (numRuns=0; !done; numRuns++){ //repeat the pruning cycle until no more rotamers are pruned	

        	reportMode.preRun(numRuns);
				
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
            int newNumPrunedRot = precompMat.getPruneMat().countPrunedRCs();
            int newNumPrunedPairs = 0;
            if ((useFlags)||(algOption>=3)) //pairs pruning is performed
                    newNumPrunedPairs = precompMat.getPruneMat().countPrunedPairs();

            
            if( (newNumPrunedRot==numPrunedRot) && (newNumPrunedPairs==numPrunedPairs) && (!preDACS) ) { 
            //no more rotamers pruned, so perform the computationally-expensive 2-sp split-Pruner and pairs

                    if ( algOption>=3 ) { //simple Goldstein pairs
                        dee.prune("GOLDSTEIN PAIRS FULL");
                        
                        if( useTriples )
                            dee.prune("GOLDSTEIN TRIPLES");
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
                    newNumPrunedRot = precompMat.getPruneMat().countPrunedRCs();
                    newNumPrunedPairs = 0;
                    if ((useFlags)||(algOption>=3)) //pairs pruning is performed
                        newNumPrunedPairs = precompMat.getPruneMat().countPrunedPairs();
            }
            
            int numPrunedRotThisRun = newNumPrunedRot - numPrunedRot;
            int numPrunedPairsThisRun = newNumPrunedPairs - numPrunedPairs;
            
            reportMode.postRun(numRuns, numPrunedRotThisRun, numPrunedPairsThisRun);
            
            if(numPrunedRotThisRun==0 && numPrunedPairsThisRun==0)
                done = true;
            
            numPrunedRot = newNumPrunedRot;
            numPrunedPairs = newNumPrunedPairs;
        }

        long pruneTime = System.currentTimeMillis() - startTime;
        
        reportMode.conclusion(numRuns, numPrunedRot, numPrunedPairs, pruneTime);
    }

    public void setOnlyGoldstein(boolean onlyGoldstein) {
        this.onlyGoldstein = onlyGoldstein;
    }
    
}
