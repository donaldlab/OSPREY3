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

package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.astar.GMECMutSpace;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.externalMemory.Queue;
import edu.duke.cs.osprey.tupexp.LUTESettings;

/**
 *
 * A more complicated version of the SimpleGMECFinder, supporting DEE and related algorithms (EPIC, LUTE, etc.)
 *
 * 
 * @author mhall44
 */
public class DEEGMECFinder extends SimpleGMECFinder {
    
    public static class Builder extends SimpleGMECFinder.Builder {
        
        EnergyMatrix emat;
        SimpleConfSpace confSpace;
        EnergyCalculator ecalc;
        String name;
        
        public Builder(EnergyMatrix emat, SimpleConfSpace confSpace, EnergyCalculator ecalc, ConfEnergyCalculator confEcalc, String name){
            super(null,confEcalc);//DEEGMECFinder generates its confSearch later (need to make some choices)
            this.emat  = emat;//might just do this as part of precompMat
            this.confSpace = confSpace;
            this.ecalc = ecalc;
            this.name = name;
        }
		
		
	public DEEGMECFinder build() {
            return new DEEGMECFinder(
                    emat,
                    confSpace,
                    ecalc,
                    confEcalc,
                    pruner,
                    logPrinter,
                    consolePrinter,
                    printIntermediateConfsToConsole,
                    useExternalMemory,
                    name
            );
	}
    }
    
    public EnergyMatrix emat;
    public SimpleConfSpace confSpace;
    public EnergyCalculator ecalc;
    public double I0=5;
    public double Ew=0;
    public boolean doIMinDEE = true;//may infer from flexibility?
    String name;
    
    String GMECMutFile = null;//if not null, limits the space of sequences considered for the GME
    
    public EPICSettings epicSettings = new EPICSettings();//DEBUG!!  maybe have a func to turn on EPIC, in a default way??  others too
    public LUTESettings luteSettings = new LUTESettings();
    public PruningSettings pruningSettings = new PruningSettings();
    
    private PrecomputedMatrices precompMat = null;//maybe toss in emat with these actually

    private DEEGMECFinder(EnergyMatrix emat, SimpleConfSpace confSpace, EnergyCalculator ecalc, ConfEnergyCalculator confEcalc, GMECFinder.ConfPruner pruner, 
                ConfPrinter logPrinter, ConfPrinter consolePrinter, boolean printIntermediateConfsToConsole, boolean useExternalMemory,
                String name){
            super(null,confEcalc,pruner,logPrinter,consolePrinter,printIntermediateConfsToConsole,true,useExternalMemory,null);
            this.emat = emat;
            this.confSpace = confSpace;
            this.ecalc = ecalc;
            this.name = name;
    }
    
    
    
    public ConfSearch.EnergiedConf calcGMEC() {
        return calcGMEC(I0).poll();
    }
    
    
    public Queue.FIFO<ConfSearch.EnergiedConf> calcGMEC(double interval) {
        //important: the interval argument here is specifically the iMinDEE interval.  
        //Ew can be added separately
        //SimpleGMECFinder doesn't care becaues it doesn't prune (so feed it Ew)
        
        System.out.println("Calculating GMEC with interval = " + interval);
        
        boolean printEPICEnergy = epicSettings.shouldWeUseEPIC();
        
        // 11/11/2015 JJ: This logic belongs out here. A function that does nothing if a flag is false should 
        // have its flag promoted outside of the function, unless it's used multiple times. In that case
        // the function needs to be named accordingly.
        if (epicSettings.shouldWeUseEPIC()) {
            checkEPICThresh2(interval);//Make sure EPIC thresh 2 matches current interval
        }

        //precompute the pruning and maybe EPIC or LUTE matrices
        //must be done separately for each round of iMinDEE
        precompMat = new PrecomputedMatrices(Ew+interval, Ew, name, emat, confSpace, ecalc, confEcalc,
            epicSettings, luteSettings, pruningSettings);
        search = prepareConfSearch();
        
        
        Queue.FIFO<ConfSearch.EnergiedConf> goodConfs = super.find(Ew);//assuming pruning and confSearch are right, find the GMEC
        if(printEPICEnergy){
            System.out.println("GMEC EPIC energy: "+precompMat.epicMat.minimizeEnergy(goodConfs.peek().getAssignments()));
        }
        //and everything within Ew of it.  
        
        //things like minScoreConf can be made fields in SimpleGMECFinder
        
        if(doIMinDEE){//iMinDEE...figure out if a second round is needed
            double lowestBound = lowestPairwiseBound();
            double firstGMECEnergy = goodConfs.peek().getEnergy();
            
            // could the minGMEC have been pruned due to a pruning interval that's too small?
            if (firstGMECEnergy > lowestBound + interval) {

                // yeah, it could have been. we can't prove minEnergyConf is the minGMEC
                // we have to pick a new interval and try again
                System.out.println("Pruning interval is too small. minGMEC could have been pruned.");
                System.out.println("Will estimate new interval based on conformations evaluated so far and restart");


                double nextInterval = firstGMECEnergy - lowestBound;

                // pad the new interval a bit to avoid numerical instability
                double intervalPad = 0.001;
                //if(searchSpace.useVoxelG)//GMEC energy has statistical error (~0.1 kcal/mol maybe)
                //    intervalPad = 0.2;
                //Let's not expose continuous entropy to Python yet
                nextInterval += intervalPad;

                return calcGMEC(nextInterval);
            }
        }
        
        return goodConfs;//this should have everything within Ew of the GMEC
    }
    
    private double lowestPairwiseBound() {
        //Make a new search and find the lowest bound
        ConfSearch lbSearch;
        if(needHigherOrderTerms() || epicSettings.shouldWeUseEPIC() || GMECMutFile!=null){//need a separate search for the lower bound
            lbSearch = ConfTree.makeFull(precompMat, parseGMECMutFile(), 
                    false, false, new EPICSettings(), confSpace.getNumPos());
        }
        else
            lbSearch = new ConfAStarTree.Builder(emat, new RCs(precompMat.pruneMat)).build();
        System.out.println("Computing lowest pairwise-minimized bound");
        ScoredConf lbConf = lbSearch.nextConf();
        double bound = lbConf==null ? Double.POSITIVE_INFINITY : lbConf.getScore();
        System.out.println("Lowest bound: "+bound);
        return bound;
    }
    
    
    private ConfSearch prepareConfSearch(){
        //get the conf search ready
        //search refers to a SearchProblem here...change
        if (needHigherOrderTerms() || epicSettings.shouldWeUseEPIC() || GMECMutFile!=null) {
                // if we need higher-order or EPIC terms, use the old A* code
            return ConfTree.makeFull(precompMat, parseGMECMutFile(), 
                    luteSettings.shouldWeUseLUTE(), epicSettings.shouldWeUseEPIC(), epicSettings, confSpace.getNumPos());
        }
        else if(luteSettings.shouldWeUseLUTE())
            return new ConfAStarTree.Builder(precompMat.luteMat, new RCs(precompMat.pruneMat)).build();
        else
            return new ConfAStarTree.Builder(emat, new RCs(precompMat.pruneMat)).build();
    }
    
    private boolean needHigherOrderTerms(){
        //need higher-order energy terms
        if(luteSettings.shouldWeUseLUTE())
            return precompMat.luteMat.hasHigherOrderTerms();
        else
            return emat.hasHigherOrderTerms();
    }
    

    
    private void checkEPICThresh2(double curInterval){
        if(curInterval+Ew>epicSettings.EPICThresh2){//need to raise EPICThresh2 
            //to the point that we can guarantee no continuous component of the GMEC
            //or desired ensemble will need to reach it
            System.out.println("Raising EPICThresh2 to "+(curInterval+Ew)+" based on "
                    + "iMinDEE interval and energy window");
            epicSettings.EPICThresh2 = curInterval+Ew;
        }
    }
    
    public GMECMutSpace parseGMECMutFile(){
        if(GMECMutFile==null)//no GMEC mut file
            return null;
        else
            return new GMECMutSpace(GMECMutFile, confSpace);
    }
    
    //TODO:
    //For MMS/MoleculeObjectiveFunction: //handleBlocksTogetherMaybe();//DEBUG!!!
    //also must check that that moving strands are made entirely of flexible res -- makeShell maybe
    //also should the first round of imindee avoid claims like "Found GMEC"?
    //deal with everything loaded in loadData, and in defaults.cfg so people can customize their runs
    //maybe also make a function to print all the settings?
    //deal with Poisson-Boltzmann
    //EPIC
    //COMETS
    //make printEPICEnergy a SimpleGMECFinderField...false by default, DEEGMECFinder can turn it on, make sure everything is printed
        //that was printed in the real GMECFinder
    //Merge MoleculeModifierAndScorer with MoleculeObjectiveFunction
    //Make base conf e tuple expander so can reuse pruning code
    
}
