package edu.duke.cs.osprey.control;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.multistatekstar.LMV;
import edu.duke.cs.osprey.multistatekstar.MultiStateConfigFileParser;
import edu.duke.cs.osprey.multistatekstar.MultiStateKStarTree;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
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

		MultiStateConfigFileParser[] cfps = new MultiStateConfigFileParser[numStates];

		SearchProblem[][] spsDisc = new SearchProblem[numStates][];
		SearchProblem[][] spsCont = new SearchProblem[numStates][];

		System.out.println();
		System.out.println("Checking multistate K* parameters for consistency");
		System.out.println();

		mutable2StateResNums = new ArrayList<>();
		wtSeqs = new ArrayList<>();
		for(int state=0; state<numStates; state++) {

			System.out.println();
			System.out.println("Checking state "+state+" parameters");

			cfps[state] = makeStateCfp(state, sParams);
			checkParamConsistency(state, cfps[state].getParams(), sParams);
			mutable2StateResNums.add(stateMutableRes(state, cfps[state], numTreeLevels));

			for(int subState=0; subState<mutable2StateResNums.get(state).size(); ++subState){
				handleAATypeOptions(state, subState, cfps[state]);
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

			spsCont[state] = makeStateSearchProblems(state, true, cfps[state]);//continuous flex
			spsDisc[state] = makeStateSearchProblems(state, false, cfps[state]);//discrete flex

			System.out.println();
			System.out.println("State "+state+" matrices ready");
			System.out.println();
		}

		printAllSeqs();
	}

	private void handleAATypeOptions(int state, int subState, MultiStateConfigFileParser stateCfp) {
		//Given the config file parser for a state, make sure AATypeOptions
		//matches the allowed AA types for this state

		Molecule wtMolec = PDBFileReader.readPDBFile( stateCfp.params.getValue("PDBName"), null );

		ArrayList<Integer> mutRes = mutable2StateResNums.get(state).get(subState);
		ArrayList<ArrayList<String>> subStateAAOptions = stateCfp.getAllowedAAs(mutRes);

		if(AATypeOptions==null) 
			AATypeOptions = new ArrayList<>();

		for(int mutPos=0; mutPos<mutRes.size(); mutPos++) {
			ArrayList<String> subStateResOptions = subStateAAOptions.get(mutPos);
			if(stateCfp.params.getBool("AddWT")) {
				Residue res = wtMolec.getResByPDBResNumber(String.valueOf(mutRes.get(mutPos)));
				//always add wt in pos 0
				if(StringParsing.containsIgnoreCase(subStateResOptions, res.template.name))
					subStateResOptions.remove(res.template.name);
				subStateResOptions.add(0, res.template.name);
			}

			if(AATypeOptions.size()<=state)//add storage for state
				AATypeOptions.add(new ArrayList<>());

			if(AATypeOptions.get(state).size()<=subState)//add storage for substate
				AATypeOptions.get(state).add(new ArrayList<>());

			if(AATypeOptions.get(state).get(subState).size()<=mutPos)//need to fill in based on this substate
				AATypeOptions.get(state).get(subState).add(subStateResOptions);

			//check for correspondence in AAs among states
			ArrayList<String> defaultSubStateResOptions = AATypeOptions.get(0).get(subState).get(mutPos);

			if(defaultSubStateResOptions.size()!=subStateResOptions.size()){
				throw new RuntimeException("ERROR: Current state has "+
						subStateResOptions.size()+" AA types allowed for position "+mutPos
						+" compared to "+defaultSubStateResOptions.size()+" for previous states");
			}

			for(int a=0;a<defaultSubStateResOptions.size();++a){
				String aa1 = defaultSubStateResOptions.get(a);
				String aa2 = subStateResOptions.get(a);

				//only amino acids must correspond between states
				if(!DAminoAcidHandler.isStandardLAminoAcid(aa1) || 
						!DAminoAcidHandler.isStandardLAminoAcid(aa2)) continue;

				if(!aa1.equalsIgnoreCase(aa2))
					throw new RuntimeException("ERROR: Current state has AA type "+
							aa2+" where previous states have AA type "+aa1+", at position "+mutPos);
			}
		}
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

	/**
	 * performs a sanity check on config files
	 * @param state
	 * @param sParams
	 * @param msParams
	 */
	private void checkParamConsistency(int state, ParamSet sParams, ParamSet msParams) {
		//parameter sanity check
		//check number of constraints
		int numUbConstr = sParams.getInt("NUMUBCONSTR");
		ArrayList<String> ubConstr = sParams.searchParams("UBCONSTR");
		ubConstr.remove("NUMUBCONSTR");
		if(numUbConstr != ubConstr.size())
			throw new RuntimeException("ERROR: NUMUBCONSTR != number of listed constraints");

		int numUbStates = sParams.getInt("NUMUBSTATES");
		if(numUbStates<2) throw new RuntimeException("ERROR: NUMUBSTATES must be >=2");

		String ubStateMutNums = sParams.getValue("UBSTATEMUTNUMS");
		StringTokenizer st = new StringTokenizer(ubStateMutNums);
		if(st.countTokens() != numUbStates) throw new RuntimeException("ERROR: "
				+ "the number of tokens in UBSTATEMUTNUMS should be the same as "
				+ "NUMUBSTATES");

		if(numUbStates!=sParams.searchParams("UBSTATELIMITS").size())
			throw new RuntimeException("ERROR: need an UBSTATELIMITS line for each NUMUBSTATES");

		int numMutsRes = 0;
		while(st.hasMoreTokens()) numMutsRes += Integer.valueOf(st.nextToken());
		if(numMutsRes != msParams.getInt("NUMMUTRES")) throw new RuntimeException("ERROR: "
				+"UBSTATEMUTNUMS does not sum up to NUMMUTRES");

		ArrayList<Integer> globalMutList = new ArrayList<>();
		for(int ubState=0;ubState<numUbStates;++ubState) {
			//num unbound state residues must match number of listed residues
			int numUbMutRes = Integer.valueOf(StringParsing.getToken(ubStateMutNums, ubState+1));
			String ubMutRes = sParams.getValue("UBSTATEMUT"+ubState);
			st = new StringTokenizer(ubMutRes);
			ArrayList<Integer> ubStateMutList = new ArrayList<>();
			while(st.hasMoreTokens()) ubStateMutList.add(Integer.valueOf(st.nextToken()));
			ubStateMutList = new ArrayList<>(new HashSet<>(ubStateMutList));
			if(ubStateMutList.size()!=numUbMutRes) throw new RuntimeException("ERROR: the "
					+"number of distinct mutable residues in UBSTATEMUT"+ubState+
					" is not equal to the value specified in UBSTATEMUTNUMS");

			globalMutList.addAll(ubStateMutList);

			//listed unbound state residues must be within limits
			ArrayList<Integer> ubStateLims = new ArrayList<>();
			st = new StringTokenizer(sParams.getValue("UBSTATELIMITS"+state));
			if(st.countTokens()!=2) throw new RuntimeException("ERROR: UBSTATELIMITS"+state
					+" must have 2 tokens");
			while(st.hasMoreTokens()) ubStateLims.add(Integer.valueOf(st.nextToken()));
			Collections.sort(ubStateLims);
			for(int res : ubStateMutList){
				if(res<ubStateLims.get(0) && res>ubStateLims.get(1)) throw new RuntimeException("ERROR: "
						+"mutable residue "+res+" exceeds the boundaries of UBSTATELIMITS"+state);
			}

			//ResAllowed must exist for each mutable residue
			for(int res: ubStateMutList) {
				if(sParams.getValue("RESALLOWED"+res, "").length()==0)
					throw new RuntimeException("ERROR: RESALLOWED"+res+" must be delcared");
			}
		}

		//check that all RESALLOWED is in list of mutable residues
		ArrayList<String> raKeys = sParams.searchParams("RESALLOWED");
		if(raKeys.size() > numMutsRes) {
			for(String raVal : raKeys) {
				raVal = raVal.replaceAll("RESALLOWED", "").trim();
				if(!globalMutList.contains(Integer.valueOf(raVal)))
					throw new RuntimeException("ERROR: RESALLOWED"+raVal+" is not in the list of UBSTATEMUT");
			}
		}

		//check that ubState limits are mutually exclusive 
		ArrayList<ArrayList<Integer>> ubStateLimits = new ArrayList<>();
		for(int ubState=0;ubState<numUbStates;++ubState) {
			st = new StringTokenizer(sParams.getValue("UBSTATELIMITS"+ubState));
			ArrayList<Integer> tmp = new ArrayList<Integer>();
			while(st.hasMoreTokens()) tmp.add(Integer.valueOf(st.nextToken()));
			Collections.sort(tmp);
			ubStateLimits.add(tmp);
		}
		if(ubStateLimits.get(0).get(0) <= ubStateLimits.get(1).get(1) && 
				ubStateLimits.get(1).get(0) <= ubStateLimits.get(0).get(1))
			throw new RuntimeException("ERROR: UBSTATELIMITS are not disjoint");
	}

	private void printAllSeqs(){
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
