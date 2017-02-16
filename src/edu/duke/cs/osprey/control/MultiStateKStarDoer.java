package edu.duke.cs.osprey.control;

import java.util.ArrayList;

import edu.duke.cs.osprey.multistatekstar.LMKSS;
import edu.duke.cs.osprey.multistatekstar.MultiStateKStarTree;

public class MultiStateKStarDoer {

	MultiStateKStarTree tree;//tree used for the MultiStateKStar search
	int numSeqsWanted;//number of desired sequences

	LMKSS objFcn;//objective function for MultiStateKStar search, i.e. the f-score
	LMKSS[] constraints;//constraints for search. (partial)sequences that violate constraints are pruned
	int numStates;//number of states considered
	int numTreeLevels;//number of mutable positions

	ArrayList<ArrayList<String>> AATypeOptions;//AA types allowed at each mutable position

	ArrayList<ArrayList<Integer>> mutable2StatePosNums;
	//For each state, a list of which flexible positions are mutable
	//these will be listed directly in Multistate.cfg under "STATEMUTRES0" etc.
	
	String stateArgs[][];//search arguments for each state
	
	public MultiStateKStarDoer(String args[]) {
		//fill in all the settings
        //each state will have its own config file parser
        
        System.out.println();
        System.out.println("Performing multistate K*");
        System.out.println();
        
      //check format of args
        if(!args[0].equalsIgnoreCase("-c"))
            throw new RuntimeException("ERROR: bad arguments (should start with -c)");
        
        ParamSet sParams = new ParamSet();
        sParams.setVerbosity(false);
        sParams.addParamsFromFile(args[3]);//read multistate parameters
        sParams.addDefaultParams();
        
        String defaultCFGName = args[1];//Each state could have its own KStar.cfg, but this is the default
        
        numStates = sParams.getInt("NUMSTATES");
        numTreeLevels = sParams.getInt("NUMMUTRES");
        int numStateConstr = sParams.getInt("NUMSTATECONSTR");
        
        stateArgs = new String[numStates][];
	}
}
