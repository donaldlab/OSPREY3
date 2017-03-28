package edu.duke.cs.osprey.control;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.function.ToDoubleFunction;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.multistatekstar.MSConfigFileParser;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.CMRF;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.CMRFNodeDomain;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.EnergyFunctionMap;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.Kernel;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.KernelGaussian;
import edu.duke.cs.osprey.pruning.PruningControl;

public class CMRFDoer {

	MSConfigFileParser cfp;//config file parser
	SearchProblem[] search;//searchproblems by strand
	ArrayList<ArrayList<Integer>> mutRes;//mutable residues by strand
	int numStates;

	public CMRFDoer(String args[]) {

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

	public void compute() {
		
		CMRF toy = new CMRF(2);
		toy.runToyCMRF(18, 0.15, 10);
//		toy.runToyCMRF(18, 0.15, 50);
//		toy.runToyCMRF(18, 0.15, 10000);
		System.exit(0);
		
		int numNodes = search[2].flexRes.size();
		CMRF cmrf = new CMRF(numNodes);
		EnergyFunctionMap efm = new EnergyFunctionMap(search[2], null);
		efm.populateOneBodyRCData();
		//efm.populateOneBody2Energy();
		//efm.populatePairWise2Energy();
		double[][] kDB0 = efm.getKernelDomainBounds(0);
		double[][] kDB1 = efm.getKernelDomainBounds(1);

		//test
		//unary rc
		Kernel k0 = new KernelGaussian(kDB0, 0.25);
		CMRFNodeDomain nd0 = new CMRFNodeDomain(
				efm.getNode(0).getDOFMin(),
				efm.getNode(0).getDOFMax(),
				k0,
				(point)->(efm.getNode(0).getEnergy())
				);

		//unary rc
		Kernel k1 = new KernelGaussian(kDB1, 0.25);
		CMRFNodeDomain nd1 = new CMRFNodeDomain(
				efm.getNode(1).getDOFMin(),
				efm.getNode(1).getDOFMax(),
				k1,
				(point)->(efm.getNode(1).getEnergy())
				);

		//pairwise rc
		HashMap<Integer, CMRFNodeDomain[]> h = new HashMap<>();
		h.put(0, new CMRFNodeDomain[]{nd0});
		h.put(1, new CMRFNodeDomain[]{nd1});
		
		//pairwise energy
		ToDoubleFunction<double[]>f = (point)->(efm.getPairWiseEnergy(efm.getNode(0), efm.getNode(1)));
		ToDoubleFunction<double[]>f1 = (point)->(efm.getPairWiseEnergy(efm.getNode(1), efm.getNode(0)));
		
		HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>> map1 = new HashMap<>();
		map1.put(nd1, f);
		HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> map2 = new HashMap<>();
		map2.put(nd0, map1);
		HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>> map3 = new HashMap<>();
		map3.put(1, map2);
		HashMap<Integer, HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>>> map4 = new HashMap<>();
		map4.put(0, map3);
		
		HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>> map5 = new HashMap<>();
		map5.put(nd0, f1);
		HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> map6 = new HashMap<>();
		map6.put(nd1, map5);
		HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>> map7 = new HashMap<>();
		map7.put(0, map6);
		map4.put(1, map7);

		

//		cmrf.addNodes(h, map4);
//		System.out.println();
//		System.out.println("Running SCMF");
//		cmrf.runSCMF();
//		System.out.println("Finished!");
//		
//		System.out.println();
//		System.out.println("Running TRBP");
//		cmrf.runTRBP(0);
//		System.out.println("Finished!");
	}

}
