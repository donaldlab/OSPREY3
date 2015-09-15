/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.util.ArrayList;

import edu.duke.cs.osprey.astar.ConfTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.confspace.PositionConfSpaceSuper;
import edu.duke.cs.osprey.confspace.ConfSpaceSuper;
/**
 *
 * @author hmn5
 */
public class KaDEEFinder {
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
    
    public KaDEEFinder(ConfigFileParser cfp){
        this.cfp = cfp;
        Ew = cfp.params.getDouble("Ew",0);
        doIMinDEE = cfp.params.getBool("imindee",false);
        if(doIMinDEE){
            I0 = cfp.params.getDouble("Ival",5);
        }
        
        useContFlex = cfp.params.getBool("doMinimize",false);
        useTupExp = cfp.params.getBool("UseTupExp",false);
        useEPIC = cfp.params.getBool("UseEPIC",false);
        
        checkApproxE = cfp.params.getBool("CheckApproxE",true);
        
        if(doIMinDEE && !useContFlex)
            throw new RuntimeException("ERROR: iMinDEE requires continuous flexibility");
        
        outputGMECStruct = cfp.params.getBool("OUTPUTGMECSTRUCT", false);
        
        useEllipses = cfp.params.getBool("useEllipses", false);
    }
    void doKaDEE(){
        //Calculate the GMEC
        
        double curInterval = I0;//For iMinDEE.  curInterval will need to be an upper bound

        searchSpace = cfp.getSearchProblem();
        ConfSpaceSuper confSpaceSuper = searchSpace.confSpaceSuper;

        ArrayList<ArrayList<Integer>> posToMerge = new ArrayList<>();
        for (int i=0; i<confSpaceSuper.posFlex.size();i++){
            ArrayList<Integer> newPos = new ArrayList<>();
            if (i==0){
                newPos.add(i);
                i++;
                newPos.add(i);
            }
            else{
                newPos.add(i);
            }
            posToMerge.add(newPos);
        }
        confSpaceSuper.mergePosition(posToMerge);

    }
}
