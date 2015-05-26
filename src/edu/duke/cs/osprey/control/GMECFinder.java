/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.SearchSpace;
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
    
    SearchSpace searchSpace;
    
    double Ew;//energy window for enumerating conformations: 0 for just GMEC
    double I0=0;//initial value of iMinDEE pruning interval
    boolean doIMinDEE;//do iMinDEE
    
    
    boolean useContFlex;
    boolean enumByLowerBound;//are we using a search method that enumerates by lower bound (minDEE-style)
    //or by true energy (rigid, EPIC, or tuple expander)?
    //Note for EPIC or tuple expander we may still be doing iMinDEE 
    //(used I0 to prune a lot, need to check if the conformations we enumerate violate iMinDEE condition)
    
    
    boolean outputGMECStruct;//write the GMEC structure to a PDB file
    
    boolean useEPIC = false;//FOR NOW!
    boolean useTupExp = false;//FOR NOW!
    
    boolean useEllipses = false;
    
    
    public GMECFinder (ConfigFileParser cfgP){
        //fill in all the settings
        
        cfp = cfgP;
        searchSpace = cfgP.getSearchSpace();
        
        Ew = cfgP.params.getDouble("Ew",0);
        doIMinDEE = cfgP.params.getBool("imindee",false);
        if(doIMinDEE){
            I0 = cfgP.params.getDouble("Ival",5);
        }
        
        useContFlex = cfgP.params.getBool("doMinimize",false);
        
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

        
        do {
            needToRepeat = false;
            
            //First calculate the pairwise energy matrix, if not already present
            searchSpace.loadEnergyMatrix();
            if(useEPIC)
                throw new UnsupportedOperationException("ERROR: EPIC not yet supported");
                //searchSpace.loadEPICMatrix();
            //and also the EPIC matrix if appropriate

            //Next, do DEE, which will fill in prunedRot
            //we can repeat the stuff from here down for iMinDEE
            PruningControl pruning = cfp.setupPruning(searchSpace,Ew+curInterval);
            
            
            pruning.prune();//pass in DEE options, and run the specified types of DEE
            //These will be subclasses of the abstract class DEE

            //Finally, do A*, which will output the top conformations
            ConfSearch search = initSearch(searchSpace);//e.g. new AStarTree from searchSpace & params
            //can have options to instantiate other kinds of search here too...choose based on params

            double lowestBound  = Double.POSITIVE_INFINITY;
            double lowerBound;
            
            do {
                int conf[] = search.nextConf();
                
                if(conf==null){//no confs left in tree. Effectively, infinite lower bound on remaining confs
                    lowerBound = Double.POSITIVE_INFINITY;
                }
                else {//tree not empty
                
                    double confE = getConfEnergy(conf);//MINIMIZED, EPIC, OR MATRIX E AS APPROPRIATE

                    if(enumByLowerBound)
                        lowerBound = searchSpace.lowerBound(conf);
                    else//we have effectively a 
                        lowerBound = confE;
                    
                    if(confE<bestESoFar){
                        bestESoFar = confE;
                        GMECConf = conf;
                    }

                    lowestBound = Math.min(lowestBound,lowerBound);

                    printConf(conf,confE,lowerBound,bestESoFar);
                    //if doing a fit true energy (e.g. EPIC),
                    //this is the only place we may want the actual 
                    //searchSpace.minimizedEnergy(conf).  compute if needed.
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
    }
    
    
    
    double getConfEnergy(int[] conf){
        //MINIMIZED, EPIC, OR MATRIX E AS APPROPRIATE
        //whatever is being used as the "true" energy while enumerating
        if(useContFlex){
            if(useEPIC||useTupExp)
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
        
        System.out.println("Lower bound: "+lowerBound+" Energy: "+confE+" Best so far: "+bestESoFar);
        //if doing a fit true energy (e.g. EPIC),
        //this is the only place we may want the actual 
        //searchSpace.minimizedEnergy(conf).  compute if needed.
    }
                
    
    
    
    ConfSearch initSearch(SearchSpace searchSpace){
        //initialize some kind of combinatorial search, like A*
        //FOR NOW just using A*; may also want BWM*, WCSP, or something according to settings
        return new ConfTree(searchSpace);
    }
    

    
}
