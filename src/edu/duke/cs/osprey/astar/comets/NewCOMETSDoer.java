/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar.comets;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import java.util.ArrayList;

/**
 *
 * This version of NewCOMETSDoer is intended to be usable via Python
 * 
 * 
 * @author mhall44
 */
public class NewCOMETSDoer {
    
    
    NewCOMETSTree tree;//The tree used for the COMETS search
    int numSeqsWanted;//How many sequences to enumerate
    
    LME objFcn;//objective function for the COMETS search
    LME[] constraints;//constraints for the COMETS search
    int numStates;//number of states considered
    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null;//AA types allowed at each mutable position
    
    ArrayList<ArrayList<Integer>> mutable2StatePosNums = new ArrayList<>();
    //For each state, a list of which flexible positions are mutable
    //these will be listed directly in Multistate.cfg under "STATEMUTRES0" etc.
    
        
    public NewCOMETSDoer (SimpleConfSpace[] confSpaces, PrecomputedMatrices[] precompMats, LME objFcn, LME[] constraints,
            ArrayList<ArrayList<Integer>> mutable2StatePosNums, ArrayList<ArrayList<String>> AATypeOptions,
            int numMaxMut, String wtSeq[], int numSeqsWanted,
            boolean outputGMECStructs, ConfEnergyCalculator[] confECalc) {
    	    	
        //fill in all the settings
        //each state will have its own config file parser
        
        this.objFcn = objFcn;
        this.constraints = constraints;
        this.mutable2StatePosNums = mutable2StatePosNums;
        this.AATypeOptions = AATypeOptions;
        numStates = confSpaces.length;
        numTreeLevels = AATypeOptions.size();
        this.numSeqsWanted = numSeqsWanted;

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)

        tree = new NewCOMETSTree(numTreeLevels, objFcn, constraints, 
            AATypeOptions, numMaxMut, wtSeq, numStates, confSpaces, precompMats, 
            mutable2StatePosNums, outputGMECStructs, confECalc);
    }
    
            
    
    public ArrayList<String> calcBestSequences(){
                    
        System.out.println("Performing multistate A*");
        
            
        //how many sequences to enumerate

        long startAStarTime = System.currentTimeMillis();
        
        ArrayList<String> bestSequences = new ArrayList<>();

        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
        	//this will find the best sequence and print it
        	ScoredConf conf = tree.nextConf();
        	if (conf == null) {
        		//empty sequence...indicates no more sequence possibilities
        		break;
        	} else {
                bestSequences.add(tree.seqAsString(conf.getAssignments()));
        	}
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));
        
        return bestSequences;
    }
        
        
    
}
