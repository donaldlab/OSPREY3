package edu.duke.cs.osprey.control;

import java.util.ArrayList;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.multistatekstar.MSConfigFileParser;
import edu.duke.cs.osprey.pruning.PruningControl;

public class ContinuousMRFDoer {

	MSConfigFileParser cfp;//config file parser
	SearchProblem[] search;//searchproblems by strand
	ArrayList<ArrayList<Integer>> mutRes;//mutable residues by strand
	int numStates;

	public ContinuousMRFDoer(String args[]) {

		//check format of args
		if(!args[0].equalsIgnoreCase("-c"))
			throw new RuntimeException("ERROR: bad arguments (should start with -c)");

		//read params
		cfp = new MSConfigFileParser(args);
		cfp.getParams().setVerbosity(false);
		cfp.loadData();

		//make sure weare doing continuous minimization
		if(!cfp.getParams().getBool("DOMINIMIZE")||!cfp.getParams().getBool("IMINDEE"))
			throw new RuntimeException("ERROR: CMRF requires continuous minimization. "+
					"Set IMINDEE and DOMINIMIZE to true");

		numStates = cfp.getParams().getInt("NUMOFSTRANDS")+1;

		//create searchproblem for each (un)bound state
		mutRes = getMutableRes();
		
		System.out.println();
		System.out.println("Preparing search problems and energy matrices");
		System.out.println();
		
		search = makeStateSearchProblems();
		
		System.out.println();
		System.out.println("Finished preparing search problems and energy matrices");
		System.out.println();

	}

	private ArrayList<ArrayList<Integer>> getMutableRes() {
		int state;
		ArrayList<ArrayList<Integer>> ans = new ArrayList<>();
		
		//unbound states
		for(state=0;state<numStates-1;++state) {
			ArrayList<Integer> mutRes = new ArrayList<>();
			StringTokenizer st = new StringTokenizer(cfp.getParams().getValue("STRANDMUT"+state));
			while(st.hasMoreTokens()) mutRes.add(Integer.valueOf(st.nextToken()));
			mutRes.trimToSize();
			ans.add(mutRes);
		}
		
		//bound state
		ans.add(new ArrayList<>());
		for(int unbound=0;unbound<numStates-1;++unbound) 
			ans.get(state).addAll(ans.get(unbound));
		ans.trimToSize();
		return ans;
	}

	public SearchProblem[] makeStateSearchProblems() {
		SearchProblem[] ans = new SearchProblem[numStates];
		for(int state=0; state<numStates;++state) {
			//make search problem
			ans[state] = cfp.getSearchProblem(state, state, mutRes.get(state), true);
			//load emat
			ans[state].loadEnergyMatrix();
			//prune emat
			if(!cfp.getParams().getBool("UsePoissonBoltzmann")) {
				PruningControl pc = cfp.setupPruning(ans[state], 
						cfp.getParams().getDouble("Ival")+cfp.getParams().getDouble("Ew"), 
						cfp.getParams().getBool("UseEpic"), 
						cfp.getParams().getBool("UseTupExp"));
				//silence output
				pc.setReportMode(null);
				pc.prune();
			}
			
			System.out.println();
			System.out.println("State "+state+" matrices ready");
			System.out.println();
		}
		return ans;
	}

}
