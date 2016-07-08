package edu.duke.cs.osprey.partcr;

import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.partcr.pickers.ConfPicker;
import edu.duke.cs.osprey.partcr.pickers.WalkingConfPicker;
import edu.duke.cs.osprey.partcr.scorers.RCScorer;
import edu.duke.cs.osprey.partcr.scorers.VolumeRCScorer;
import edu.duke.cs.osprey.partcr.splitters.NAryRCSplitter;
import edu.duke.cs.osprey.partcr.splitters.RCSplitter;
import edu.duke.cs.osprey.tools.TimeFormatter;

public class PartCR {
	
	private SearchProblem search;
	private double Ew;
	private SimpleEnergyCalculator ecalc;
	private List<ConfAStarNode> nodes;
	private SplitWorld splitWorld;
	private ConfPicker picker;
	private RCScorer scorer;
	private RCSplitter splitter;
	private double bestMinimizedEnergy;
	private long minimizationNs;
	private long iterationNs;
	private int numIterations;
	
	public PartCR(SearchProblem search, double Ew, SimpleEnergyCalculator ecalc, List<ConfAStarNode> nodes) {
		
		this.search = search;
		this.Ew = Ew;
		this.ecalc = ecalc;
		this.nodes = nodes;
		
		// make a separate conf space for splitting RCs
		splitWorld = new SplitWorld(search, ecalc.getEnergyFunctionGenerator(), ecalc.getShellDistribution());
		
		// init default options
		// NOTE: these options worked best for me on very limited test cases
		// but more work could certainly be done to pick better defaults
		picker = new WalkingConfPicker();
		scorer = new VolumeRCScorer();
		splitter = new NAryRCSplitter();
		
		bestMinimizedEnergy = Double.POSITIVE_INFINITY;
		minimizationNs = 0;
		iterationNs = 0;
		numIterations = 0;
	}
	
	public void setPicker(ConfPicker val) {
		picker = val;
	}
	
	public void setScorer(RCScorer val) {
		scorer = val;
	}
	
	public void setSplitter(RCSplitter val) {
		splitter = val;
	}
	
	public List<ConfAStarNode> getNodes() {
		return nodes;
	}
	
	public long getAvgMinimizationTimeNs() {
		return minimizationNs/numIterations;
	}
	
	public long getAvgIterationTimeNs() {
		return iterationNs/numIterations;
	}
	
	public void autoIterate() {
		// three strikes seems like a good default
		autoIterate(3);
	}
	
	public void autoIterate(int maxNumStrikes) {
		
		// the basic idea here is to keep iterating as long as
		// we can prune conformations faster than it takes to enumerate/minimize them
		
		int initialNumConfs = nodes.size();
		int numConfs = initialNumConfs;
		long startTimeNs = System.nanoTime();
		int numStrikes = 0;
		
		while (true) {
			
			iterate();
			
			// based on timing info, how many nodes should we prune this iteration?
			int targetPruning = (int)(getAvgIterationTimeNs()/getAvgMinimizationTimeNs());
			
			// did we prune enough to keep iterating?
			int numPruned = numConfs - nodes.size();
			numConfs = nodes.size();
			if (numPruned < targetPruning) {
				
				numStrikes++;
				boolean shouldStop = numStrikes >= maxNumStrikes;
				
				System.out.println(String.format("Pruned %d/%d conformations. %d/%d strikes. %s",
					numPruned, targetPruning, numStrikes, maxNumStrikes,
					shouldStop ? " time to stop." : ""
				));
				
				if (shouldStop) {
					// YOU'RE OUTTA HERE!!
					break;
				}
				
			} else {
			
				System.out.println(String.format("Pruned %d/%d conformations. keep iterating", numPruned, targetPruning));
			}
		}
		
		long initialTimeNs = getAvgMinimizationTimeNs()*initialNumConfs;
		long afterTimeNs = getAvgMinimizationTimeNs()*numConfs;
		long diffTimeNs = System.nanoTime() - startTimeNs;
		System.out.println(String.format("\nPartCR took %s and saved you %s for a net savings of %s!",
			TimeFormatter.format(diffTimeNs, 1),
			TimeFormatter.format(initialTimeNs - afterTimeNs, 1),
			TimeFormatter.format(initialTimeNs - afterTimeNs - diffTimeNs, 1)
		));
		System.out.println(String.format("PartCR pruned %.1f%% of low-energy conformations",
			100.0*(initialNumConfs - numConfs)/initialNumConfs
		));
	}
	
	public void iterate() {
		
		numIterations++;
		System.out.println("\nPartCR iteration " + numIterations);
		
		long startIterNs = System.nanoTime();
		
		int numPos = splitWorld.getSearchProblem().confSpace.numPos;
		
		// pick a conformation to analyze
		int[] conf = new int[numPos];
		picker.pick(nodes).getConf(conf);
		
		// analyze the conf and put our protein in the minimized conformation
		// NOTE: it's very important that the protein be in the minimized conformation to calculate bound errors correctly
		double boundEnergy = calcBoundEnergy(conf);
		System.out.println("minimizing conformation...");
		long startMinNs = System.nanoTime();
		double minimizedEnergy = calcMinimizedEnergy(conf);
		long diffMinNs = System.nanoTime() - startMinNs;
		minimizationNs += diffMinNs;
		
		if (numIterations == 1) {
			System.out.println(String.format("initial conformations: %d, estimated time to enumerate: %s",
				nodes.size(), TimeFormatter.format(getAvgMinimizationTimeNs()*nodes.size(), 1)
			));
		}
		
		bestMinimizedEnergy = Math.min(bestMinimizedEnergy, minimizedEnergy);
		
		// score all the positions
		double boundEnergyCheck = 0;
		double minimizedEnergyCheck = 0;
		TreeMap<Double,Integer> positionsByScore = new TreeMap<>();
		for (int pos=0; pos<numPos; pos++) {
			RC rcObj = splitWorld.getSearchProblem().confSpace.posFlex.get(pos).RCs.get(conf[pos]);
			
			double posBoundEnergy = calcPosBoundEnergy(conf, pos);
			double posMinimizedEnergy = calcPosMinimizedEnergy(pos);
			
			boundEnergyCheck += posBoundEnergy;
			minimizedEnergyCheck += posMinimizedEnergy;
			
			double err = posMinimizedEnergy - posBoundEnergy;
			double score = scorer.calcScore(splitWorld, rcObj, err);
			
			positionsByScore.put(score, pos);
		}
		
		// just in case...
		checkEnergy(boundEnergyCheck, boundEnergy);
		checkEnergy(minimizedEnergyCheck, minimizedEnergy);
		
		// split the RC at the position with the highest score
		System.out.println("splitting residue conformation...");
		int splitPos = positionsByScore.lastEntry().getValue();
		RC rcObj = splitWorld.getRC(splitPos, conf[splitPos]);
		List<RC> splitRCs = splitter.split(splitPos, rcObj);
		splitWorld.replaceRc(splitPos, rcObj, splitRCs);
		
		// update energy matrix
		System.out.println("updating matrices...");
		splitWorld.updateMatrices(nodes.get(0).getGScore(), bestMinimizedEnergy, Ew);
		
		// prune nodes based on the new bounds
		System.out.println("pruning conformations...");
		Iterator<ConfAStarNode> iter = nodes.iterator();
		while (iter.hasNext()) {
			ConfAStarNode node = iter.next();
			
			// use the split world to get a tighter bound
			double improvedBoundEnergy = splitWorld.improveNode(node).getGScore();
			
			// if the new bound pushes the node over the threshold, prune it
			if (improvedBoundEnergy > bestMinimizedEnergy + Ew) {
				iter.remove();
			}
		}
		
		// update timing
		long diffIterNs = System.nanoTime() - startIterNs;
		iterationNs += diffIterNs;
		
		System.out.println(String.format("conformations remaining: %d, estimated time to enumerate: %s",
			nodes.size(), TimeFormatter.format(getAvgMinimizationTimeNs()*nodes.size(), 1)
		));
	}
	
	private void checkEnergy(double observed, double expected) {
		double absErr = Math.abs(observed - expected);
		double relErr = absErr/Math.abs(expected);
		final double Epsilon = 1e-10;
		if (relErr > Epsilon) {
			throw new Error(String.format("Energies don't match. This is a bug!\n\texpected: %f\n\tobserved: %f\n\tabs err:  %f\n\trel err:  %f",
				expected, observed, absErr, relErr
			));
		}
	}

	private double calcBoundEnergy(int[] conf) {
		
		double energy = 0;
		
		EnergyMatrix emat = splitWorld.getSearchProblem().emat;
		int numPos = emat.getNumPos();
		
		for (int pos1=0; pos1<numPos; pos1++) {
			int rc1 = conf[pos1];
			
			energy += emat.getOneBody(pos1, rc1);
			
			for (int pos2=0; pos2<pos1; pos2++) {
				int rc2 = conf[pos2];
				
				energy += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}
		
		return energy;
	}
	
	private double calcMinimizedEnergy(int[] conf) {
		return search.minimizedEnergy(conf);
	}
	
	private double calcPosBoundEnergy(int[] conf, int pos1) {
		
		// for position energies, distribute all the adjacent pairwise energies over this position
		// instead of just picking the ones where pos2 < pos1
		
		EnergyMatrix emat = splitWorld.getSearchProblem().emat;
		int numPos = emat.getNumPos();
		
		int rc1 = conf[pos1];
		double energy = emat.getOneBody(pos1, rc1);
		
		for (int pos2=0; pos2<numPos; pos2++) {
			if (pos1 != pos2) {
				int rc2 = conf[pos2];
				
				energy += emat.getPairwise(pos1, rc1, pos2, rc2)/2;
			}
		}
		
		return energy;
	}
	
	private double calcPosMinimizedEnergy(int pos1) {
		
		// for position energies, distribute all the adjacent pairwise energies over this position
		// instead of just picking the ones where pos2 < pos1
		
		double energy = 0;
		
		// single energy
		energy += ecalc.getSingleEfunc(pos1).getEnergy();
		
		// pairwise energy
		for (int pos2=0; pos2<ecalc.getNumPos(); pos2++) {
			if (pos2 != pos1) {
				energy += ecalc.getPairEfunc(pos1, pos2).getEnergy()/2;
			}
		}
		
		return energy;
	}
}
