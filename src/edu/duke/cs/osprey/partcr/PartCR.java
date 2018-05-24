/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.partcr;

import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
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
	private List<ScoredConf> confs;
	private SplitWorld splitWorld;
	private ConfPicker picker;
	private RCScorer scorer;
	private RCSplitter splitter;
	private double bestMinimizedEnergy;
	private long minimizationNs;
	private long iterationNs;
	private int numIterations;
	
	public PartCR(SearchProblem search, ForcefieldParams ffparams, double Ew, List<ScoredConf> confs) {
		
		this.search = search;
		this.Ew = Ew;
		this.ecalc = new SimpleEnergyCalculator.Cpu(ffparams, search.confSpace, search.shellResidues);
		this.confs = confs;
		
		// make a separate conf space for splitting RCs
		splitWorld = new SplitWorld(search, ffparams);
		
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
	
	public List<ScoredConf> getConfs() {
		return confs;
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
		
		int initialNumConfs = confs.size();
		int numConfs = initialNumConfs;
		long startTimeNs = System.nanoTime();
		int numStrikes = 0;
		
		while (true) {
			
			iterate();
			
			// based on timing info, how many nodes should we prune this iteration?
			int targetPruning = (int)(getAvgIterationTimeNs()/getAvgMinimizationTimeNs());
			int numPruned = numConfs - confs.size();
			numConfs = confs.size();
			
			// are there even that many conformations left?
			if (numConfs <= targetPruning) {
				System.out.println(String.format("Pruned %d/%d conformations. Only %d conformations left, time to stop.",
					numPruned, targetPruning, numConfs
				));
				break;
			}
			
			// did we prune enough to keep iterating?
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
		long savingsNs = initialTimeNs - afterTimeNs;
		long netSavingsNs = savingsNs - diffTimeNs;
		System.out.println(String.format("\nPartCR took %s and saved %s of minimization time for a net %s",
			TimeFormatter.format(diffTimeNs, 1),
			TimeFormatter.format(savingsNs, 1),
			String.format("%s of %s", (netSavingsNs > 0) ? "SAVINGS" : "LOSS", TimeFormatter.format(Math.abs(netSavingsNs), 1))
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
		
		// pick a conformation to analyze, and translate it into the split world
		ScoredConf pickedConf = picker.pick(confs);
		ScoredConf translatedPickedConf = splitWorld.translateConf(pickedConf);
		
		// analyze the conf and put our protein in the minimized conformation
		// NOTE: it's very important that the protein be in the minimized conformation to calculate bound errors correctly
		System.out.println("minimizing conformation...");
		long startMinNs = System.nanoTime();
		EnergiedConf analyzeConf = new EnergiedConf(translatedPickedConf, calcMinimizedEnergy(pickedConf.getAssignments()));
		long diffMinNs = System.nanoTime() - startMinNs;
		minimizationNs += diffMinNs;
		bestMinimizedEnergy = Math.min(bestMinimizedEnergy, analyzeConf.getEnergy());
		
		if (numIterations == 1) {
			System.out.println(String.format("initial conformations: %d, estimated time to enumerate: %s",
				confs.size(), TimeFormatter.format(getAvgMinimizationTimeNs()*confs.size(), 1)
			));
		}
		
		// score all the positions (using the split world)
		double boundEnergyCheck = search.emat.getConstTerm();
		double minimizedEnergyCheck = 0;
		TreeMap<Double,Integer> positionsByScore = new TreeMap<>();
		for (int pos=0; pos<numPos; pos++) {
			RC rcObj = splitWorld.getSearchProblem().confSpace.posFlex.get(pos).RCs.get(analyzeConf.getAssignments()[pos]);
			
			double posBoundEnergy = calcPosBoundEnergy(analyzeConf.getAssignments(), pos);
			double posMinimizedEnergy = calcPosMinimizedEnergy(pos);
			
			boundEnergyCheck += posBoundEnergy;
			minimizedEnergyCheck += posMinimizedEnergy;
			
			double err = posMinimizedEnergy - posBoundEnergy;
			double score = scorer.calcScore(splitWorld, rcObj, err);
			
			positionsByScore.put(score, pos);
		}
		
		// make sure the sums of our position energies match the conf energies
		// otherwise, our re-distribution of energy terms is wrong and it's a bug
		checkEnergy(boundEnergyCheck, analyzeConf.getScore());
		checkEnergy(minimizedEnergyCheck, analyzeConf.getEnergy());
		
		// split the RC at the position with the highest score
		System.out.println("splitting residue conformation...");
		int splitPos = positionsByScore.lastEntry().getValue();
		RC rcObj = splitWorld.getRC(splitPos, analyzeConf.getAssignments()[splitPos]);
		List<RC> splitRCs = splitter.split(splitPos, rcObj);
		splitWorld.replaceRc(splitPos, rcObj, splitRCs);
		
		// NOTE: since SplitWorld uses lazy evaluation for computing the energy matrix,
		// the actual energy calculation that gets split across later DEE/A* calls.
		// so it's actually counter-productive to separate resizeMatrices() (which calculates energies)
		// and improveBound() (which calls A*) in the log, or for timing purposes
		System.out.println("calculating energies and pruning conformations...");
		
		splitWorld.resizeMatrices();
		
		// prune nodes based on the new bounds
		Iterator<ScoredConf> iter = confs.iterator();
		while (iter.hasNext()) {
			ScoredConf conf = iter.next();
			
			// use the split world to get a tighter bound
			double improvedBoundEnergy = splitWorld.translateConf(conf).getScore();
			
			if (improvedBoundEnergy > bestMinimizedEnergy + Ew) {
				
				// new bound is above pruning threshold, so prune the conf
				iter.remove();
			}
		}
		
		// update timing
		long diffIterNs = System.nanoTime() - startIterNs;
		iterationNs += diffIterNs;
		
		System.out.println(String.format("finished iteration in %s", TimeFormatter.format(diffIterNs, 1)));
		System.out.println(String.format("conformations remaining: %d, estimated time to enumerate: %s",
			confs.size(), TimeFormatter.format(getAvgMinimizationTimeNs()*confs.size(), 1)
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
		energy += ecalc.makeSingleEfunc(pos1).getEnergy();
		
		// pairwise energy
		for (int pos2=0; pos2<ecalc.confSpace.numPos; pos2++) {
			if (pos2 != pos1) {
				energy += ecalc.makePairEfunc(pos1, pos2).getEnergy()/2;
			}
		}
		
		return energy;
	}
}
