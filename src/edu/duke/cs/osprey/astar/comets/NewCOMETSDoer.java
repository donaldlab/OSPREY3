/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

