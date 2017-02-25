package edu.duke.cs.osprey.control;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.multistatekstar.InputValidation;
import edu.duke.cs.osprey.multistatekstar.KStarRatio;
import edu.duke.cs.osprey.multistatekstar.MultiStateKStarSettings;
import edu.duke.cs.osprey.multistatekstar.LMV;
import edu.duke.cs.osprey.multistatekstar.MultiStateConfigFileParser;
import edu.duke.cs.osprey.multistatekstar.MultiStateKStarTree;
import edu.duke.cs.osprey.multistatekstar.MultiStateSearchProblem;
import edu.duke.cs.osprey.multistatekstar.SearchProblemSettings;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.tools.Stopwatch;
import edu.duke.cs.osprey.tools.StringParsing;

public class MultiStateKStarDoer {

	MultiStateKStarTree tree;//tree used for the MultiStateKStar search
	int numSeqsWanted;//number of desired sequences
	int numMaxMut;//max number of mutations from wt

	LMV objFcn;//objective function for MultiStateKStar search, i.e. the f-score
	LMV[] constraints;//constraints for search. (partial)sequences that violate constraints are pruned
	int numStates;//number of states considered
	int numTreeLevels;//number of mutable positions

	ArrayList<String[]> wtSeqs;//bound state wild type sequences for each state

	ArrayList<ArrayList<ArrayList<ArrayList<String>>>> AATypeOptions;//AA types allowed at 
	//each mutable position for each substate

	ArrayList<ArrayList<ArrayList<Integer>>> mutable2StateResNums;
	//For each state, a list of which flexible positions are mutable
	//these will be listed directly in system cfg files

	String stateArgs[][];//search arguments for each state

	MultiStateConfigFileParser[] cfps;//config file parsers for each state

	SearchProblem[][] searchDisc;//continuous search problems
	SearchProblem[][] searchCont;//discrete search problems
	
	public ConfEnergyCalculator.Async[][] ecalcs;//global energy calculator objects

	public MultiStateKStarDoer(String args[]) {
		//fill in all the settings
		//each state will have its own config file parser

		System.out.println();
		System.out.println("Performing multistate K*");
		System.out.println();

		//check format of args
		if(!args[0].equalsIgnoreCase("-c"))
			throw new RuntimeException("ERROR: bad arguments (should start with -c)");

		// multistate spec parameters
		ParamSet sParams = new ParamSet();
		sParams.setVerbosity(false);
		sParams.addParamsFromFile(args[4]);//read multistate parameters
		sParams.addDefaultParams();

		numStates = sParams.getInt("NUMSTATES");
		numTreeLevels = sParams.getInt("NUMMUTRES");
		numMaxMut = sParams.getInt("NUMMAXMUT");
		int numConstr = sParams.getInt("NUMSTATECONSTR");

		stateArgs = new String[numStates][];

		objFcn = new LMV(sParams.getValue("OBJFCN"), numStates);

		constraints = new LMV[numConstr];
		for(int constr=0; constr<numConstr; constr++)
			constraints[constr] = new LMV(sParams.getValue("STATECONSTR"+constr), numStates);
		// might need to adjust these with wt later

		cfps = new MultiStateConfigFileParser[numStates];

		searchDisc = new SearchProblem[numStates][];
		searchCont = new SearchProblem[numStates][];

		System.out.println();
		System.out.println("Checking multistate K* parameters for consistency");
		System.out.println();

		mutable2StateResNums = new ArrayList<>();
		AATypeOptions = new ArrayList<>();
		wtSeqs = new ArrayList<>();
		InputValidation inputValidation = new InputValidation(AATypeOptions, mutable2StateResNums);
		for(int state=0; state<numStates; state++) {

			System.out.println();
			System.out.println("Checking state "+state+" parameters");

			cfps[state] = makeStateCfp(state, sParams);
			inputValidation.handleStateParams(state, cfps[state].getParams(), sParams);
			mutable2StateResNums.add(stateMutableRes(state, cfps[state], numTreeLevels));

			for(int subState=0; subState<mutable2StateResNums.get(state).size(); ++subState){
				inputValidation.handleAATypeOptions(state, subState, cfps[state]);
				//get bound substate wt sequence
				if(subState==mutable2StateResNums.get(state).size()-1)
					wtSeqs.add(cfps[state].getWtSeq(mutable2StateResNums.get(state).get(subState)));
			}

			System.out.println("State "+state+" parameters checked");
			System.out.println();
		}

		System.out.println();
		System.out.println("Preparing search problems and matrices for multistate K*");
		System.out.println();

		for(int state=0; state<numStates; state++) {

			searchCont[state] = makeStateSearchProblems(state, true, cfps[state]);//continuous flex
			searchDisc[state] = makeStateSearchProblems(state, false, cfps[state]);//discrete flex

			System.out.println();
			System.out.println("State "+state+" matrices ready");
			System.out.println();
		}
		
		System.out.println();
		System.out.println("Preparing energy calculators for multistate K*");
		System.out.println();
		
		ecalcs = makeEnergyCalculators(true);
		
		System.out.println();
		System.out.println("Finished preparing energy calculators for multistate K*");
		System.out.println();
	}
	
	private ConfEnergyCalculator.Async[][] makeEnergyCalculators(boolean continuous) {
		SearchProblem[][] sps = continuous ? searchCont : searchDisc;
		ConfEnergyCalculator.Async[][] ans = new ConfEnergyCalculator.Async[sps.length][];
		
		for(int state=0;state<sps.length;++state) {
			SearchProblem[] statesps = sps[state];
			ans[state] = new ConfEnergyCalculator.Async[statesps.length];
			
			for(int substate=0;substate<ans[state].length;++substate) {
				ans[state][substate] = MultiStateKStarSettings.makeEnergyCalculator(cfps[state], sps[state][substate]);
			}
		}
		
		return ans;
	}
	
	private void cleanupEnergyCalculators() {
		if(ecalcs == null) return;
		
		for(int state=0;state<ecalcs.length;++state) {
			for(int substate=0;substate<ecalcs[state].length;++substate) {
				ConfEnergyCalculator.Async ecalc = ecalcs[state][substate];
				if(ecalc != null) {
					ecalcs[state][substate].cleanup();
					ecalcs[state][substate] = null;
				}
			}
		}
		
		ecalcs = null;
	}
	
	public void cleanup() {
		cleanupEnergyCalculators();
	}

	/**
	 * Verify algorithm results by doing exhaustive search
	 */
	void exhaustiveMultistateSearch(){

		System.out.println();
		System.out.println("Checking MultiStateKStar by exhaustive search");
		System.out.println();

		Stopwatch stopwatch = new Stopwatch().start();
		
		ArrayList<ArrayList<String[]>> seqList = listAllSeqs();
		int numSeqs = seqList.get(0).size();
		BigDecimal[][] stateKSRs = new BigDecimal[numStates][numSeqs];

		for(int state=0; state<numStates; state++){
			for(int seqNum=0; seqNum<numSeqs; seqNum++){
				stateKSRs[state][seqNum] = calcStateKSRatio(state, seqList.get(state).get(seqNum));
			}
		}

		System.out.println();
		System.out.println("Finished checking MultiStateKStar by exhaustive search in "+stopwatch.getTime(2));
		System.out.println();
	}

	private BigDecimal calcStateKSRatio(int state, String[] boundStateAATypes) {

		//get arraylist formatted sequence for each substate
		ArrayList<ArrayList<ArrayList<String>>> subStateAATypes = new ArrayList<>();
		for(int subState=0;subState<mutable2StateResNums.get(state).size();++subState){
			//mutable2StateResNums.get(state).get(subState)) contains substate flexible residues
			//last substate is the bound state
			subStateAATypes.add(new ArrayList<>());
			ArrayList<Integer> subStateResNums = mutable2StateResNums.get(state).get(subState);
			ArrayList<Integer> boundStateResNums = mutable2StateResNums.get(state).get(mutable2StateResNums.get(state).size()-1);
			for(int resNum : subStateResNums){
				int index = boundStateResNums.indexOf(resNum);
				ArrayList<String> aa = new ArrayList<>(); aa.add(boundStateAATypes[index]);
				subStateAATypes.get(subState).add(aa);
			}
		}

		ParamSet sParams = cfps[state].getParams();

		//make LMVs
		int numUbConstr = sParams.getInt("NUMUBCONSTR");
		int numPartFuncs = sParams.getInt("NUMUBSTATES")+1;
		LMV[] sConstraints = new LMV[numUbConstr];
		for(int constr=0;constr<numUbConstr;constr++)
			sConstraints[constr] = new LMV(sParams.getValue("UBCONSTR"+constr), numPartFuncs);

		//populate search problems
		MultiStateSearchProblem[] sps = new MultiStateSearchProblem[numPartFuncs];
		for(int subState=0;subState<numPartFuncs;++subState){

			SearchProblemSettings spSet = new SearchProblemSettings();
			spSet.AATypeOptions = subStateAATypes.get(subState);
			ArrayList<String> mutRes = new ArrayList<>();
			for(int i:mutable2StateResNums.get(state).get(subState)) mutRes.add(String.valueOf(i));
			spSet.mutRes = mutRes;

			sps[subState] = sParams.getBool("DOMINIMIZE") ? new MultiStateSearchProblem(searchCont[state][subState], spSet)
					: new MultiStateSearchProblem(searchDisc[state][subState], spSet);
		}

		//make k* settings
		MultiStateKStarSettings ksSettings = new MultiStateKStarSettings();
		ksSettings.pruningWindow = sParams.getDouble("IVAL") + sParams.getDouble("EW");
		ksSettings.targetEpsilon = sParams.getDouble("EPSILON");
		ksSettings.stericThreshold = sParams.getDouble("STERICTHRESH");
		ksSettings.cfp = cfps[state];
		ksSettings.search = sps;
		ksSettings.constraints = sConstraints;
		ksSettings.ecalcs = ecalcs[state];
		ksSettings.sequence = boundStateAATypes;
		ksSettings.isReportingProgress = true;

		KStarRatio ksRatio = new KStarRatio(ksSettings);
		ksRatio.compute(Integer.MAX_VALUE);
		return ksRatio.getKStarRatio();
	}

	/**
	 * creates and prunes all-seq energy matrices for each substate within a state
	 * @param state
	 * @param cont
	 * @param stateCfp
	 * @return
	 */
	private SearchProblem[] makeStateSearchProblems(int state, boolean cont, 
			MultiStateConfigFileParser stateCfp) {

		ParamSet sParams = stateCfp.getParams();
		int numUbStates = sParams.getInt("NUMUBSTATES");		

		SearchProblem[] subStateSps = new SearchProblem[numUbStates+1];

		for(int subState=0;subState<subStateSps.length;++subState) {
			subStateSps[subState] = stateCfp.getSearchProblem(state, subState, 
					mutable2StateResNums.get(state).get(subState), cont);

			//make emats
			subStateSps[subState].loadEnergyMatrix();

			//prune
			if(!sParams.getBool("UsePoissonBoltzmann")) {
				PruningControl pc = stateCfp.setupPruning(subStateSps[subState], 
						sParams.getDouble("Ival")+sParams.getDouble("Ew"), 
						sParams.getBool("UseEpic"), 
						sParams.getBool("UseTupExp"));
				//silence output
				pc.setReportMode(null);
				pc.prune();
			}

			System.out.println();
			System.out.println("State "+state+"."+subState+" matrices ready");
			System.out.println();
		}

		return subStateSps;
	}

	/**
	 * creates mutable residue list for all substates in all states
	 * @param state
	 * @param stateCfp
	 * @param numTreeLevels
	 * @return
	 */
	private ArrayList<ArrayList<Integer>> stateMutableRes(int state, MultiStateConfigFileParser stateCfp, int numTreeLevels){
		ParamSet sParams = stateCfp.getParams();
		int numUbStates = sParams.getInt("NUMUBSTATES");
		ArrayList<ArrayList<Integer>> m2s = new ArrayList<>();
		for(int ubState=0;ubState<=numUbStates;++ubState) m2s.add(new ArrayList<>());

		for(int ubState=0;ubState<numUbStates;++ubState) {
			StringTokenizer st = new StringTokenizer(sParams.getValue("UBSTATEMUT"+ubState));
			while(st.hasMoreTokens()) m2s.get(ubState).add(Integer.valueOf(st.nextToken()));
			//append to complex residues
			m2s.get(numUbStates).addAll(m2s.get(ubState));
		}

		if(m2s.get(numUbStates).size()!=numTreeLevels){
			throw new RuntimeException("ERROR: SeqTree has "+numTreeLevels+" mutable positions "
					+" but "+m2s.size()+" are listed for state "+state);
		}

		return m2s;
	}

	/**
	 * makes a config file parser for each state
	 * @param state
	 * @param sParams
	 * @return
	 */
	private MultiStateConfigFileParser makeStateCfp(int state, ParamSet sParams) {
		//We expect input of the form KStar0.cfg System0.cfg DEE0.cfg
		String stateConfigFiles = sParams.getValue("STATECFGFILES"+state);
		String stateKStFile = StringParsing.getToken(stateConfigFiles, 1);
		String stateSysFile = StringParsing.getToken(stateConfigFiles, 2);
		String stateDEEFile = StringParsing.getToken(stateConfigFiles, 3);

		stateArgs[state] = new String[] {"-c", stateKStFile, "n/a", stateSysFile, stateDEEFile};

		MultiStateConfigFileParser stateCfp = new MultiStateConfigFileParser(stateArgs[state], false);
		//verbosity is false by default; too many arguments clutter the display
		stateCfp.loadData();

		return stateCfp;
	}

	protected void printAllSeqs(){
		ArrayList<ArrayList<String[]>> stateSeqLists = listAllSeqs();
		for(int state=0;state<stateSeqLists.size();++state){

			int numSeqs=stateSeqLists.get(state).size();

			System.out.println();
			System.out.println("State"+state+": "+numSeqs+" sequences with <= "+numMaxMut+" mutation(s) from wild-type");
			System.out.println();

			for(String[] seq : stateSeqLists.get(state)){ 
				for(String aa : seq) System.out.print(aa+" ");
				System.out.println();
			}
		}
	}

	private ArrayList<ArrayList<String[]>> listAllSeqs(){
		//for all bound states, list all possible sequences for the mutable residues,
		//based on AATypeOptions
		ArrayList<ArrayList<String[]>> ans = new ArrayList<>();

		//pre-allocate buffer
		String[] buf = new String[numTreeLevels];

		for(int state=0;state<numStates;++state) {
			int subState = AATypeOptions.get(state).size()-1;
			ArrayList<ArrayList<String>> subStateAATypeOptions = AATypeOptions.get(state).get(subState);
			ArrayList<String[]> stateOutput = new ArrayList<>();

			//get allowed sequences for this state's bound complex
			listAllSeqsHelper(subStateAATypeOptions, stateOutput, wtSeqs.get(state), buf, 0, 0);
			stateOutput.trimToSize();
			ans.add(stateOutput);
		}

		ans.trimToSize();
		return ans;
	}

	private void listAllSeqsHelper(ArrayList<ArrayList<String>> subStateAATypeOptions, 
			ArrayList<String[]> stateOutput, String[] wt, String[] buf, int depth, int dist){
		//List all sequences for the subset of mutable positions with max distance
		//from wt starting at depth=0 and going to the last mutable position
		if(depth==numTreeLevels){
			String[] seq = new String[numTreeLevels];
			System.arraycopy(buf, 0, seq, 0, numTreeLevels);
			stateOutput.add(seq);
			return;
		}

		for(int aaIndex=0;aaIndex<subStateAATypeOptions.get(depth).size();++aaIndex){
			buf[depth]=subStateAATypeOptions.get(depth).get(aaIndex);
			int nDist=buf[depth].equalsIgnoreCase(wt[depth]) ? dist : dist+1;
			if(nDist>numMaxMut) continue;
			listAllSeqsHelper(subStateAATypeOptions, stateOutput, wt, buf, depth+1, nDist);
		}
	}

	public ArrayList<String> calcBestSequences() {

		System.out.println();
		System.out.println("Performing multistate K*");
		System.out.println();

		//how many sequences to enumerate

		long startKStarTime = System.currentTimeMillis();

		ArrayList<String> bestSequences = new ArrayList<>();

		for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
			int seq[] = tree.nextConf().getAssignments();//this will find the best sequence and print it
			if(seq == null)//empty sequence...indicates no more sequence possibilities
				break;
			else
				bestSequences.add(tree.seqAsString(seq));
		}

		long stopTime = System.currentTimeMillis();
		System.out.println("Sequence enumeration time: "+((stopTime-startKStarTime)/(60.0*1000.0)));

		//exhaustiveMultistateSearch();
		return bestSequences;
	}
}
