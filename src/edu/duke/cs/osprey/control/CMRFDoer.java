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
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.RCData;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.RCDatum;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.SCMF;
import edu.duke.cs.osprey.partitionfunctionbounds.continuous.TRBP;
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
		/*
		if(!cfp.getParams().getBool("DOMINIMIZE")||!cfp.getParams().getBool("IMINDEE"))
			throw new RuntimeException("ERROR: CMRF requires continuous minimization. "+
					"Set IMINDEE and DOMINIMIZE to true");
	    */
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

		//		CMRF.runToyCMRF(1, 0.15, 1);
		//		toy.runToyCMRF(18, 0.15, 50);
		//		toy.runToyCMRF(18, 0.15, 10000);
		//		System.exit(0);

		EnergyFunctionMap efm = new EnergyFunctionMap(search[2], null);
		RCDatum.EFM = efm;
		RCData.EFM = efm;
		efm.populateOneBodyRCData();
		//efm.populateOneBody2Energy();
		//efm.populatePairWise2Energy();

		//kernel domain bounds
		double[][][] kdb = new double[efm.getNumResidues()][][];
		//kernels
		Kernel[] k = new Kernel[kdb.length];
		//cmrf node domains
		CMRFNodeDomain[] nd = new CMRFNodeDomain[k.length];
		//node domain hashmap
		HashMap<Integer, CMRFNodeDomain[]> ndMap = new HashMap<>();

		for(int residue=0; residue<kdb.length; ++residue) {

			kdb[residue] = efm.getKernelDomainBounds(residue);
			
			k[residue] = new KernelGaussian(kdb[residue], 1.0);

			final int fresidue = residue;
			nd[residue] = new CMRFNodeDomain(efm.getNode(residue).getDOFMin(), 
					efm.getNode(residue).getDOFMax(), 
					k[residue],
					(point)->(efm.getNode(fresidue).getEnergy(point)));

			ndMap.put(residue, new CMRFNodeDomain[]{nd[residue]});
		}

		HashMap<Integer, 
		HashMap<Integer, 
		HashMap<CMRFNodeDomain, 
		HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>>> edgeMap = 
		new HashMap<>();

		for (int i=0; i<nd.length; i++) {
			HashMap<Integer, HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>>> map1 = 
					new HashMap<>();
			for (int j=0; j<nd.length; j++) {
				if (i==j) { continue; } 
				HashMap<CMRFNodeDomain, HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>>> map2 = new HashMap<>();
				for (CMRFNodeDomain d1 : ndMap.get(i)) {
					HashMap<CMRFNodeDomain, ToDoubleFunction<double[]>> map3 = new HashMap<>();
					for (CMRFNodeDomain d2 : ndMap.get(j)) {
						final int fi = i, fj = j;
						map3.put(d2, (point)->efm.addPairWise(efm.getNode(fi), efm.getNode(fj)).getEnergy(point));
					}
					map2.put(d1, map3);
				}
				map1.put(j,  map2);
			}
			edgeMap.put(i, map1);
		}
		
		CMRF c = new CMRF(nd.length);
		c.addNodes(ndMap, edgeMap);

		SCMF s = new SCMF(c);
		double logZLB = s.runSCMF();
		
		TRBP t = new TRBP(c);
		double logZUB = t.runTRBP(2); // no iterations of LBP
		
		double[] ret = new double[2];
		ret[0] = logZLB; ret[1] = logZUB;
	}

}
