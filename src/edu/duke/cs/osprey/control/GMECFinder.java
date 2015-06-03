/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.util.ArrayList;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.pruning.PruningControl;

/**
 *
 * @author mhall44
 */
public class GMECFinder {
    
    //Many of the parameters that control the optimization process (like what DEE algs to use,
    //whether to use MPLP, etc. can be stored in this class
    //I propose that they be grouped into classes, similar to "EPICSettings," that handle settings for
    //particular aspects of functionality (energy function parameters, DEE/A* alg settings, etc.)
    //But we can also just store things as fields if people prefer
    
    //KStarCalculator will be set up similarly to this class
    
    ConfigFileParser cfp;
    
    SearchProblem searchSpace;
    
    double Ew;//energy window for enumerating conformations: 0 for just GMEC
    double I0=0;//initial value of iMinDEE pruning interval
    boolean doIMinDEE;//do iMinDEE
    
    
    boolean useContFlex;
    //boolean enumByPairwiseLowerBound;//are we using a search method that 
    //enumerates by pairwise lower bound (minDEE-style)
    //or by (possibly approximated) true energy (rigid, EPIC, or tuple expander)?
    
    
    boolean outputGMECStruct;//write the GMEC structure to a PDB file
    
    boolean useEPIC = false;
    boolean useTupExp = false;
    
    boolean checkApproxE = true;//Use the actual energy function to evaluate
    //each enumerated conformation rather than just relying on the EPIC or tup-exp approximation
    
    boolean useEllipses = false;
    
    
    public GMECFinder (ConfigFileParser cfgP){
        //fill in all the settings
        
        cfp = cfgP;
        
        Ew = cfgP.params.getDouble("Ew",0);
        doIMinDEE = cfgP.params.getBool("imindee",false);
        if(doIMinDEE){
            I0 = cfgP.params.getDouble("Ival",5);
        }
        
        useContFlex = cfgP.params.getBool("doMinimize",false);
        useTupExp = cfgP.params.getBool("UseTupExp",false);
        useEPIC = cfgP.params.getBool("UseEPIC",false);
        
        checkApproxE = cfgP.params.getBool("CheckApproxE",true);
        
        if(doIMinDEE && !useContFlex)
            throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");
        
        outputGMECStruct = cfgP.params.getBool("OUTPUTGMECSTRUCT", false);
        
        useEllipses = cfgP.params.getBool("useEllipses", false);
        
        //FOR NOW minimization-aware is by lower bound...
        enumByLowerBound = useContFlex;
    }
    
    
   
    void calcGMEC(){
        //Calculate the GMEC
        
        double curInterval = I0;//For iMinDEE.  curInterval will need to be an upper bound
        //on GMEC-lowestBound
        //but we can start with an estimate I0 and raise if needed
        boolean needToRepeat;
        
        //GMEC is the lowest-energy conformation enumerated in this whole search
        int GMECConf[] = null;
        double bestESoFar = Double.POSITIVE_INFINITY;
        
        searchSpace = cfp.getSearchProblem();
        
        do {
            needToRepeat = false;
            
            //initialize a search problem with current Ival
            checkEPICThresh2(curInterval);//Make sure EPIC thresh 2 matches current interval

            precomputeMatrices(Ew+curInterval);//precompute the energy, pruning, and maybe EPIC or tup-exp matrices
            //must be done separately for each round of iMinDEE
            
            //Finally, do A*, which will output the top conformations
            ConfSearch search = initSearch(searchSpace);//e.g. new AStarTree from searchSpace & params
            //can have options to instantiate other kinds of search here too...choose based on params
            
            double lowestBound  = Double.POSITIVE_INFINITY;
            if( (useEPIC||useTupExp) && doIMinDEE)//lowest bound must be calculated without EPIC/tup exp, to ensure valid iMinDEE interval
                lowestBound = lowestPairwiseBound(searchSpace);
            
            double lowerBound;
            int conformationCount=0;
            
            System.out.println();
            System.out.println("BEGINNING CONFORMATION ENUMERATION");
            System.out.println();
            
            do {
                int conf[] = search.nextConf();
                
                if(conf==null){//no confs left in tree. Effectively, infinite lower bound on remaining confs
                    lowerBound = Double.POSITIVE_INFINITY;
                }
                else {//tree not empty
                                        
                    double confE = getConfEnergy(conf);//MINIMIZED, EPIC, OR MATRIX E AS APPROPRIATE
                    //this is the true energy; if !checkApproxE may be an approximation to it
                    
                    
                    lowerBound = confE;//the lower bound and confE are the same in rigid
                    //or non-checkApproxE EPIC and tup-exp calculations
                    //lowerBound is always what we enumerate in order of
                    if(useContFlex){
                        if(useTupExp||useEPIC){
                            if(checkApproxE)//confE is the "true" energy function
                                lowerBound = searchSpace.approxMinimizedEnergy(conf);
                        }
                        else
                            lowerBound = searchSpace.lowerBound(conf);
                    }
                    
                    if(confE<bestESoFar){
                        bestESoFar = confE;
                        GMECConf = conf;
                        System.out.println("New best energy: "+confE);
                    }

                    lowestBound = Math.min(lowestBound,lowerBound);

                    System.out.println("");
                    System.out.println("Time taken: "+((System.currentTimeMillis()-startTime)/1000));
                    System.out.println("CONFORMATION "+(++conformationCount));
                    printConf(conf,confE,lowerBound,bestESoFar);
                }
                
                if(doIMinDEE){//if there are no conformations with minimized energy
                    //within curInterval of the lowestBound
                    //then we have to repeat with a higher curInterval
                    if( (lowerBound > lowestBound + curInterval)
                            && (bestESoFar > lowestBound + curInterval) ){
                        
                        curInterval = bestESoFar - lowestBound + 0.001;//give it a little buffer for numerical issues
                        System.out.println("Raising pruning interval to "+curInterval);
                        needToRepeat=true;
                        break;
                    }
                }
                
                if(bestESoFar==Double.POSITIVE_INFINITY){//no conformations found
                    System.out.println("A* returned no conformations.");
                    break;
                }

            } while( bestESoFar+Ew >= lowerBound );//lower bound above GMEC + Ew...can stop enumerating
            
        } while(needToRepeat);
        
        if(outputGMECStruct && GMECConf!=null)
            searchSpace.outputMinimizedStruct( GMECConf, searchSpace.name+".GMEC.pdb" );
        
        System.out.println("GMEC calculation complete.  ");
        double gmecE = searchSpace.minimizedEnergy(GMECConf);
        System.out.println("GMEC energy: "+gmecE);
    }
    
    
    private void precomputeMatrices(double pruningInterval){
            //Precalculate TupleMatrices needed for GMEC computation.  Some of these may already be computed.  
            //All of these matrices except the basic pairwise energy matrix are pruning-dependent:
            //we can prune conformations whose energies are within pruningInterval
            //of the lowest pairwise lower bound
        
        
            //First calculate the pairwise energy matrix, if not already present
            searchSpace.loadEnergyMatrix();
            
            //Next, do DEE, which will fill in the pruning matrix
            PruningControl pruning = cfp.setupPruning(searchSpace,pruningInterval);
            
            pruning.prune();//pass in DEE options, and run the specified types of DEE            
            
            //precomputing EPIC or tuple-expander matrices is much faster
            //if only done for unpruned RCs.  Less RCs to handle, and the fits are far simpler.  
            if(useEPIC)
                searchSpace.loadEPICMatrix();
            if(useTupExp)//preferably do this one EPIC loaded (much faster if can fit to EPIC)
                searchSpace.loadTupExpEMatrix();
    }
    
    
    double getConfEnergy(int[] conf){
        //MINIMIZED, EPIC, OR MATRIX E AS APPROPRIATE
        //whatever is being used as the "true" energy while enumerating
        if(useContFlex){
            if( (useEPIC||useTupExp) && (!checkApproxE) )
                return searchSpace.approxMinimizedEnergy(conf);
            else
                return searchSpace.minimizedEnergy(conf);
        }
        else
            return searchSpace.lowerBound(conf);
    }
    
    void printConf(int[] conf, double confE, double lowerBound, double bestESoFar){
        
        System.out.println("ENUMERATING CONFORMATION.  RCs (residue-based numbers):");
        for(int rc : conf)
            System.out.print(rc + " ");
        System.out.println();
        
        System.out.println("Residue types: ");
        for(int pos=0; pos<searchSpace.confSpace.numPos; pos++){
            String resType = searchSpace.confSpace.posFlex.get(pos).RCs.get(conf[pos]).AAType;
            System.out.print( resType + " " );
        }
        System.out.println();
        
        System.out.println("Rotamer numbers: ");
        for(int pos=0; pos<searchSpace.confSpace.numPos; pos++){
            int rotNum = searchSpace.confSpace.posFlex.get(pos).RCs.get(conf[pos]).rotNum;
            System.out.print( rotNum + " " );
        }
        System.out.println();
        
        
        String energyStatement = "Lower bound/enumeration energy: "+lowerBound+" Energy: "+confE+" Best so far: "+bestESoFar;
        //Lower bound/enumeration energy is what we enumerate in order of
        //(either a lower bound on the actual energy, or the same as Energy)
        
        System.out.println(energyStatement);
    }
                
    
    
    
    ConfSearch initSearch(SearchProblem searchSpace){
        //initialize some kind of combinatorial search, like A*
        //FOR NOW just using A*; may also want BWM*, WCSP, or something according to settings
        return new ConfTree(searchSpace);
    }
    
    
    private double lowestPairwiseBound(SearchProblem prob){
        //In an EPIC calculation, our enumeration will probably include much less conformations,
        //but for iMinDEE purposes we still need to know what our lowest bound would have been
        //if we enumerated w/o EPIC (i.e., what is the minimum energy calculated using the lower-bound energy matrix)
        
        System.out.println();
        System.out.println("Calculating no-EPIC lowest energy bound");
        
        SearchProblem probNoEPIC = new SearchProblem(prob);//shallow copy
        probNoEPIC.useEPIC = false;
        probNoEPIC.useTupExpForSearch = false;
        
        ConfSearch searchNoEPIC = initSearch(probNoEPIC);
        int LBConf[] = searchNoEPIC.nextConf();//lowest conf for this search
        double LB = probNoEPIC.lowerBound(LBConf);
        
        System.out.println("No-EPIC lowest energy bound: "+LB);
        System.out.println();
        
        return LB;
    }

    
    private void checkEPICThresh2(double curInterval){
                         
        if(useEPIC){
            if(curInterval+Ew>searchSpace.epicSettings.EPICThresh2){//need to raise EPICThresh2 
                //to the point that we can guarantee no continuous component of the GMEC
                //or desired ensemble will need to reach it
                System.out.println("Raising EPICThresh2 to "+(curInterval+Ew)+" based on "
                        + "iMinDEE interval and energy window");
                searchSpace.epicSettings.EPICThresh2 = curInterval+Ew;
            }
        }
    }
    
}
