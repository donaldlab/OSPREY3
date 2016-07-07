package edu.duke.cs.osprey.partcr;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.order.AStarOrder;
import edu.duke.cs.osprey.astar.conf.order.StaticScoreHMeanAStarOrder;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.MPLPPairwiseHScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.mplp.NodeUpdater;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.ShellDistribution;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.partcr.pickers.ConfPicker;
import edu.duke.cs.osprey.partcr.pickers.WalkingConfPicker;
import edu.duke.cs.osprey.partcr.scorers.RCScorer;
import edu.duke.cs.osprey.partcr.scorers.VolumeRCScorer;
import edu.duke.cs.osprey.partcr.splitters.NAryRCSplitter;
import edu.duke.cs.osprey.partcr.splitters.RCSplitter;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;

public class Playground extends TestBase {
	
	public static void main(String[] args)
	throws Exception {
		
		// config
		initDefaultEnvironment();
		
		MultiTermEnergyFunction.setNumThreads(4);
		
		// make the search problem
		//String aaNames = "ALA VAL LEU ILE PHE TYR TRP CYS MET SER THR LYS ARG HIE HID ASP GLU ASN GLN GLY";
		//String aaNames = "ALA VAL LEU ILE GLU ASN GLN GLY";
		String aaNames = "ALA VAL LEU ILE";
		//String aaNames = "ALA";
		String flexRes = "38 39 40 41 42 43 44";
		//String flexRes = "38 39 40 41";
		//String flexRes = "41 42 43 44";
		ArrayList<String> flexResList = new ArrayList<>(Arrays.asList(flexRes.split(" ")));
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<flexResList.size(); i++) {
			allowedAAs.add(new ArrayList<>(Arrays.asList(aaNames.split(" "))));
		}
		boolean doMinimize = true;
		boolean addWt = false;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		SearchProblem search = new SearchProblem(
			"test", "test/1CC8/1CC8.ss.pdb", 
			flexResList, allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion,
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots
		);
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		// NOTE: AllOnSingles is much faster than the others, even though it gives a little bit looser bounds
		// the speed/tightness tradeoff probably favors AllOnSingles, so let's just use that
		ShellDistribution dist = ShellDistribution.AllOnSingles;
		//ShellDistribution dist = ShellDistribution.Even;
		//ShellDistribution dist = ShellDistribution.HalfSinglesHalfPairs;
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues, dist);
		
		// compute the energy and dof matrices
		File ematFile = new File(String.format("/tmp/emat.partcr.%s.dat", dist));
		if (ematFile.exists()) {
			System.out.println("\nReading matrices...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		} else {
			System.out.println("\nComputing matrices...");
			search.emat = ecalc.calcEnergyMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// don't do any pruning, it's pointless in the continuous case anyway
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		
		// get the min bound conformation
		RCs rcs = new RCs(search.pruneMat);
		System.out.println("Finding min bound GMEC among " + search.confSpace.getNumConformations().doubleValue() + " conformations...");
		AStarOrder order = new StaticScoreHMeanAStarOrder();
		AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
		ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
		tree.initProgress();
		int[] minBoundConf = tree.nextConf();
		System.out.println("min bound conf: " + Arrays.toString(minBoundConf));
		
		// keep track of improvements over time
		List<Integer> minBoundCounts = new ArrayList<>();
		List<Long> timesMs = new ArrayList<>();
		timesMs.add(0L);
		
		// get the min bound and minimized energies
		System.out.println("Minimizing...");
		double minBoundEnergy = search.emat.getInternalEnergy(new RCTuple(minBoundConf));
		double minimizedEnergy = search.minimizedEnergy(minBoundConf);
		double boundError = minimizedEnergy - minBoundEnergy;
		System.out.println("min bound energy: " + minBoundEnergy);
		System.out.println("minimized energy: " + minimizedEnergy);
		System.out.println("Bound error:      " + boundError);
			
		// estimate how many confs there are between the min bound and the min GMEC
		System.out.println("estimating conformations...");
		List<ConfAStarNode> nodes = enumerateConformations(tree, minimizedEnergy);
		int minBoundCount = nodes.size();
		minBoundCounts.add(minBoundCount);
		
		// TODO: some enterprising student could try to optimize the PartCR configuration
		// i.e., see what settings/heuristics work best on a wide variety of design problems
		
		// configure PartCR
		
		//ConfPicker picker = new FirstConfPicker();
		ConfPicker picker = new WalkingConfPicker();
		
		//RCScorer scorer = new NopRCErrorScorer();
		//RCScorer scorer = new SplitsRCScorer();
		RCScorer scorer = new VolumeRCScorer();
		
		RCSplitter splitter = new NAryRCSplitter();
		//RCSplitter splitter = new BinaryRCSplitter();
		
		double Ew = 0;
		PartCR pcr = new PartCR(search, Ew, ecalc, nodes);
		pcr.setPicker(picker);
		pcr.setScorer(scorer);
		pcr.setSplitter(splitter);
		
		// run it!!
		pcr.autoIterate();
	}
	
	private static List<ConfAStarNode> enumerateConformations(ConfAStarTree tree, double minimizedEnergy) {
		List<ConfAStarNode> nodes = new ArrayList<>();
		while (true) {
			ConfAStarNode node = tree.nextLeafNode();
			if (node == null) {
				break;
			}
			nodes.add(node);
			if (node.getGScore() >= minimizedEnergy) {
				break;
			}
		}
		return nodes;
	}
}
