package edu.duke.cs.osprey.partcr;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.TreeMap;

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
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCIndexMap;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.DofMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.Result;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.ShellDistribution;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;

public class Playground extends TestBase {
	
	private static class Confs {
		
		private List<int[]> confs;
		private List<Double> minBounds;
		private List<Double> minEnergies;
		
		public Confs() {
			confs = new ArrayList<>();
			minBounds = new ArrayList<>();
			minEnergies = new ArrayList<>();
		}
		
		public void add(int[] conf) {
			add(conf, null, null);
		}
		
		public void add(int[] conf, Double minBound, Double minEnergy) {
			confs.add(conf);
			minBounds.add(minBound);
			minEnergies.add(minEnergy);
		}
		
		public int size() {
			return confs.size();
		}
		
		public int[] getConf(int i) {
			return confs.get(i);
		}
		
		public double getMinBound(int i) {
			return minBounds.get(i);
		}
		
		public double getMinEnergy(int i) {
			return minEnergies.get(i);
		}
		
		public boolean hasEnergies(int i) {
			return minBounds.get(i) != null;
		}
		
		public void removeEnergies(int i) {
			minBounds.set(i, null);
			minEnergies.set(i, null);
		}
		
		public void setEnergies(int i, double minBound, double minEnergy) {
			minBounds.set(i, minBound);
			minEnergies.set(i, minEnergy);
		}
		
		public int getMinBoundConfIndex() {
			double bestEnergy = Double.POSITIVE_INFINITY;
			int bestConfi = -1;
			for (int i=0; i<confs.size(); i++) {
				if (minBounds.get(i) < bestEnergy) {
					bestEnergy = minBounds.get(i);
					bestConfi = i;
				}
			}
			return bestConfi;
		}
		
		public double getMinEnergy() {
			double minEnergy = Double.POSITIVE_INFINITY;
			for (double energy : minEnergies) {
				minEnergy = Math.min(minEnergy, energy);
			}
			return minEnergy;
		}
	}
	
	private static class SplitWorld {
	
		private ConfSpace confSpace;
		private RCSplits splits;
		private List<RCIndexMap> rcMaps;
		private EnergyMatrix emat;
		private DofMatrix dofmat;
		private SimpleEnergyCalculator ecalc;
		private EnergyFunction efunc;
		
		public SplitWorld(SearchProblem search, DofMatrix dofmat, EnergyFunctionGenerator egen, ShellDistribution dist) {
			
			// copy the confspace, the positions, and the rcs
			this.confSpace = new ConfSpace(search.confSpace);
			for (int pos=0; pos<this.confSpace.numPos; pos++) {
				this.confSpace.posFlex.set(pos, new PositionConfSpace(search.confSpace.posFlex.get(pos)));
			}
			
			this.splits = new RCSplits(confSpace);
			
			this.rcMaps = new ArrayList<>();
			for (int pos=0; pos<this.confSpace.numPos; pos++) {
				this.rcMaps.add(null);
			}
			
			this.emat = new EnergyMatrix(search.emat);
			this.dofmat = dofmat;
			this.ecalc = new SimpleEnergyCalculator(egen, confSpace, search.shellResidues, dist);
			this.efunc = egen.fullConfEnergy(confSpace, search.shellResidues);
		}
		
		public void replaceRc(int pos, RC rc, List<RC> splitRCs) {
			
			// keep track of which children go to which parents
			splits.split(rc, splitRCs);
			
			// split the rc and save the index map for updateMatrices()
			RCIndexMap map = confSpace.posFlex.get(pos).replaceRC(rc, splitRCs);
			rcMaps.set(pos, map);
		}
		
		public RCIndexMap getRcMap(int pos) {
			return rcMaps.get(pos);
		}
		
		public void updateMatrices() {
		
			int numPos = confSpace.numPos;
			
			assert (emat != null);
			assert (dofmat != null);
			
			// copy as many energies as possible from the old matrix
			// to a new matrix that is sized for the new conf space (with the added rotamers)
			EnergyMatrix newEmat = new EnergyMatrix(confSpace, emat.getPruningInterval());
			DofMatrix newDofmat = new DofMatrix(confSpace);
			
			for (int pos1=0; pos1<numPos; pos1++) {
				RCIndexMap rcMap1 = rcMaps.get(pos1);
				for (int rc1=0; rc1<newEmat.getNumConfAtPos(pos1); rc1++) {

					// get the old rc1, if any
					Integer oldRc1 = rc1;
					if (rcMap1 != null) {
						oldRc1 = rcMap1.newToOld(rc1);
					}
					
					// one-body
					if (oldRc1 != null) {
						
						// copy values from the old matrices
						newDofmat.setOneBody(pos1, rc1, dofmat.getOneBody(pos1, oldRc1));
						newEmat.setOneBody(pos1, rc1, emat.getOneBody(pos1, oldRc1));
						
					} else {
						
						// calculate new values
						Result result = ecalc.calcSingle(pos1, rc1);
						newDofmat.setOneBody(pos1, rc1, result.getDofValues());
						newEmat.setOneBody(pos1, rc1, result.getEnergy());
					}
					
					// pairwise
					for (int pos2=0; pos2<pos1; pos2++) {
						RCIndexMap rcMap2 = rcMaps.get(pos2);
						for (int rc2=0; rc2<newEmat.getNumConfAtPos(pos2); rc2++) {
							
							// get the old rc2, if any
							Integer oldRc2 = rc2;
							if (rcMap2 != null) {
								oldRc2 = rcMap2.newToOld(rc2);
							}
							
							if (oldRc1 != null && oldRc2 != null) {
								
								// copy values from the old matrices
								newDofmat.setPairwise(pos1, rc1, pos2, rc2, dofmat.getPairwise(pos1, oldRc1, pos2, oldRc2));
								newEmat.setPairwise(pos1, rc1, pos2, rc2, emat.getPairwise(pos1, oldRc1, pos2, oldRc2));
							
							} else {
								
								// calculate new values
								Result result = ecalc.calcPair(pos1, rc1, pos2, rc2);
								newDofmat.setPairwise(pos1, rc1, pos2, rc2, result.getDofValues());
								newEmat.setPairwise(pos1, rc1, pos2, rc2, result.getEnergy());
							}
						}
					}
				}
			}
			
			// clear the maps now that we're done with them
			for (int pos1=0; pos1<numPos; pos1++) {
				rcMaps.set(pos1, null);
			}
				
			emat = newEmat;
			dofmat = newDofmat;
		}
		
		public double minimize(int[] conf) {
			return confSpace.minimizeEnergy(conf, efunc, null);
		}
		
		public RCs makeRCs(ConfAStarNode leafNode) {
			
			int[] conf = new int[confSpace.numPos];
			leafNode.getConf(conf);
			
			// make a pruning matrix that only leaves the parent rc at unsplit positions
			// and the split rcs at the split positions
			PruningMatrix pruneMat = new PruningMatrix(confSpace, 0);
			for (int pos=0; pos<confSpace.numPos; pos++) {
				
				List<RC> rcsAtPos = confSpace.posFlex.get(pos).RCs;
				RCSplits.RCInfo info = splits.getRCInfo(pos, conf[pos]);
				
				if (info.isSplit()) {
					
					// prune all but the splits
					for (int rc=0; rc<rcsAtPos.size(); rc++) {
						if (!info.isChild(rc)) {
							pruneMat.setOneBody(pos, rc, true);
						}
					}
					
				} else {
					
					// prune all but the parent
					for (int rc=0; rc<rcsAtPos.size(); rc++) {
						if (!info.isParent(rc)) {
							pruneMat.setOneBody(pos, rc, true);
						}
					}
				}
			}
			
			return new RCs(pruneMat);
		}

		public ConfAStarNode improveNode(ConfAStarNode node) {
			RCs subRcs = makeRCs(node);
			AStarOrder order = new StaticScoreHMeanAStarOrder();
			AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), emat, 1, 0.0001);
			ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(emat), hscorer, subRcs);
			return tree.nextLeafNode();
		}
	}
	
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
		//String flexRes = "38 39 40 41 42 43 44";
		String flexRes = "38 39 40 41";
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
		int numPos = search.confSpace.numPos;
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		// NOTE: AllOnSingles is much faster than the others, even though it gives a little bit looser bounds
		// the speed/tightness tradeoff probably favors AllOnSingles, so let's just use that
		ShellDistribution dist = ShellDistribution.AllOnSingles;
		//ShellDistribution dist = ShellDistribution.Even;
		//ShellDistribution dist = ShellDistribution.HalfSinglesHalfPairs;
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(egen, search.confSpace, search.shellResidues, dist);
		
		// compute the energy and dof matrices
		File ematFile = new File(String.format("/tmp/emat.partcr.%s.dat", dist));
		File dofmatFile = new File(String.format("/tmp/dofmat.partcr.%s.dat", dist));
		DofMatrix dofmat = null;
		if (ematFile.exists() && dofmatFile.exists()) {
			System.out.println("\nReading matrices...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
			dofmat = (DofMatrix)ObjectIO.readObject(dofmatFile.getAbsolutePath(), true);
		} else {
			System.out.println("\nComputing matrices...");
			search.emat = new EnergyMatrix(search.confSpace, 0);
			dofmat = new DofMatrix(search.confSpace);
			ecalc.calcMatrices(search.emat, dofmat);
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
			ObjectIO.writeObject(dofmat, dofmatFile.getAbsolutePath());
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
		List<String> splitRotamers = new ArrayList<>();
		List<Double> minBoundEnergies = new ArrayList<>();
		List<Double> minimizedEnergies = new ArrayList<>();
		List<Integer> minBoundCounts = new ArrayList<>();
		
		// get the min bound and minimized energies
		System.out.println("Minimizing...");
		double minBoundEnergy = search.emat.getInternalEnergy(new RCTuple(minBoundConf));
		double minimizedEnergy = search.minimizedEnergy(minBoundConf);
		double boundError = minimizedEnergy - minBoundEnergy;
		double bestMinimizedEnergy = minimizedEnergy;
		System.out.println("min bound energy: " + minBoundEnergy);
		System.out.println("minimized energy: " + minimizedEnergy);
		System.out.println("Bound error:      " + boundError);
		minBoundEnergies.add(minBoundEnergy);
		minimizedEnergies.add(minimizedEnergy);
			
		// estimate how many confs there are between the min bound and the min GMEC
		System.out.println("estimating conformations...");
		int minBoundCount = countBounds(tree, minimizedEnergy);
		System.out.println("confs to enumerate: " + minBoundCount);
		minBoundCounts.add(minBoundCount);
		
		// make a separate conf space for splitting RCs
		SplitWorld splitWorld = new SplitWorld(search, dofmat, egen, dist);
		
		Confs confs = new Confs();
		confs.add(minBoundConf, minBoundEnergy, minimizedEnergy);
		
		//RCSplitter splitter = new NAryRCSplitter();
		RCSplitter splitter = new BinaryRCSplitter();
		
		int numPartCrIters = 20;
		for (int i=0; i<numPartCrIters; i++) {
			System.out.println("\nPartCR iteration " + (i+1));
			//try {
			
				String splitRotsDesc = "";
				
				// analyze RC DOFs
				boolean reportDofs = false;
				for (int pos1=0; pos1<numPos; pos1++) {
					int rc1 = minBoundConf[pos1];
					RC rc1Obj = splitWorld.confSpace.posFlex.get(pos1).RCs.get(rc1);
					int numDofs1 = rc1Obj.DOFs.size();
					assert (splitWorld.dofmat.getOneBody(pos1, rc1).length == numDofs1);
					
					boolean areDofsInRange = true;
					
					if (reportDofs) {
						System.out.println(String.format("Analysis: %d,%d, %d dofs", pos1, rc1, numDofs1));
						for (int j=0; j<numDofs1; j++) {
							System.out.println(String.format("\tdof %d: %f,%f", j, rc1Obj.DOFmin.get(j), rc1Obj.DOFmax.get(j)));
						}
					}
					double[] dofVals = new double[rc1Obj.DOFs.size()];
					for (int j=0; j<dofVals.length; j++) {
						dofVals[j] = rc1Obj.DOFs.get(j).getCurVal();
					}
					if (reportDofs) {
						System.out.println(String.format("\tmin dofs: %s  %b", Arrays.toString(dofVals), areDofsInRange(dofVals, rc1Obj)));
					}
					areDofsInRange &= areDofsInRange(dofVals, rc1Obj);
					dofVals = splitWorld.dofmat.getOneBody(pos1, rc1);
					if (reportDofs) {
						System.out.println(String.format("\tbound dofs single: %s  %b", Arrays.toString(dofVals), areDofsInRange(dofVals, rc1Obj)));
					}
					areDofsInRange &= areDofsInRange(dofVals, rc1Obj);
					
					for (int pos2=0; pos2<numPos; pos2++) {
						if (pos2 != pos1) {
							int rc2 = minBoundConf[pos2];
							int numDofs2 = splitWorld.confSpace.posFlex.get(pos2).RCs.get(rc2).DOFs.size();
							assert (splitWorld.dofmat.getOneBody(pos2, rc2).length == numDofs2);
							
							double[] pairDofs = splitWorld.dofmat.getPairwise(pos1, rc1, pos2, rc2);
							assert (numDofs1 + numDofs2 == pairDofs.length);
							if (pos2 < pos1) {
								pairDofs = Arrays.copyOfRange(pairDofs, 0, numDofs1);
								assert (pairDofs.length == rc1Obj.DOFs.size());
							} else {
								pairDofs = Arrays.copyOfRange(pairDofs, numDofs2, pairDofs.length);
								assert (pairDofs.length == rc1Obj.DOFs.size());
							}
							if (reportDofs) {
								System.out.println(String.format("\tbound dofs pair %d,%d: %s  %b", pos2, rc2, Arrays.toString(pairDofs), areDofsInRange(pairDofs, rc1Obj)));
							}
							areDofsInRange &= areDofsInRange(pairDofs, rc1Obj);
						}
					}
					
					assert (areDofsInRange);
				}
				
				// analyze bound error
				boolean reportBounds = false;
				double minBoundEnergyCheck = 0;
				double minimizedEnergyCheck = 0;
				for (int pos1=0; pos1<numPos; pos1++) {
					int rc1 = minBoundConf[pos1];
					
					// single energies
					double singleBoundEnergy = splitWorld.emat.getOneBody(pos1, rc1);
					double singleMinimizedEnergy = ecalc.getSingleEfunc(pos1).getEnergy();
					
					double singleErr = singleMinimizedEnergy - singleBoundEnergy;
					if (reportBounds) {
						System.out.println(String.format("%5d:  %10f  %10f  %10f", pos1, singleBoundEnergy, singleMinimizedEnergy, singleErr));
					}
					minBoundEnergyCheck += singleBoundEnergy;
					minimizedEnergyCheck += singleMinimizedEnergy;
					
					// pairwise energies
					for (int pos2=0; pos2<pos1; pos2++) {
						int rc2 = minBoundConf[pos2];
						
						double pairwiseBoundEnergy = splitWorld.emat.getPairwise(pos1, rc1, pos2, rc2);
						double pairwiseMinimizedEnergy = ecalc.getPairEfunc(pos1, pos2).getEnergy();
						
						double pairwiseErr = pairwiseMinimizedEnergy - pairwiseBoundEnergy;
						if (reportBounds) {
							System.out.println(String.format("%2d,%2d:  %10f  %10f  %10f", pos1, pos2, pairwiseBoundEnergy, pairwiseMinimizedEnergy, pairwiseErr));
						}
						minBoundEnergyCheck += pairwiseBoundEnergy;
						minimizedEnergyCheck += pairwiseMinimizedEnergy;
					}
				}
				assert (Math.abs(minBoundEnergy - minBoundEnergyCheck) < 1e-12);
				assert (Math.abs(minimizedEnergy - minimizedEnergyCheck) < 1e-12);
				
				
				// PER-POSITION SPLITTING
				
				// sort all the positions by errors
				TreeMap<Double,Integer> positionsByError = new TreeMap<>();
				for (int pos=0; pos<numPos; pos++) {
					int rot = splitWorld.confSpace.posFlex.get(pos).RCs.get(minBoundConf[pos]).rotNum;
					double posBoundEnergy = getPosBoundEnergy(splitWorld.emat, minBoundConf, pos);
					double posMinEnergy = getPosMinEnergy(ecalc, pos);
					double posErr = posMinEnergy - posBoundEnergy;
					System.out.println(String.format("pos=%2d, rot=%3d, bound=%10f, energy=%10f, err=%10f",
						pos, rot, posBoundEnergy, posMinEnergy, posErr
					));
					positionsByError.put(posErr, pos);
				}
				
				List<Integer> positionsInOrder = new ArrayList<>(positionsByError.values());
				Collections.reverse(positionsInOrder);
				
				Integer splitPos = null;
				RC splitRC = null;
				int numOldRCs = 0;
				List<RC> splitRCs = null;
				
				// split the rc with the worst error
				for (int pos1 : positionsInOrder) {
				
					splitPos = pos1;
					PositionConfSpace posConfSpace = splitWorld.confSpace.posFlex.get(pos1);
					numOldRCs = posConfSpace.RCs.size();
					int rc1 = minBoundConf[pos1];
					splitRC = posConfSpace.RCs.get(rc1);
					int rotNum = splitRC.rotNum;
					
					splitRCs = splitter.split(pos1, splitRC);
					if (splitRCs != null) {
						
						// TEMP
						System.out.println(String.format("split %d-%d", pos1, rotNum));
						splitRotsDesc += String.format("%d-%d ", pos1, rotNum);
					
						splitWorld.replaceRc(splitPos, splitRC, splitRCs);
						
						// TEMP
						System.out.println("split: " + splitRC.RCIndex + " -> " + getRcNums(splitRCs));
						
						break;
						
					} else {
						System.out.println(String.format("skipped rot %d-%d", pos1, rotNum));
					}
				}
				
				// no splits? nothing to do
				if (splitRCs == null) {
					break;
				}
				
				splitRotamers.add(splitRotsDesc);
				
				assert (splitPos != null);
				assert (splitRC != null);
				assert (numOldRCs > 0);
				
				System.out.println(String.format("prev conf:  %30s  bound=%12f  min=%12f  err=%12f",
					Arrays.toString(minBoundConf), minBoundEnergy, minimizedEnergy, (minimizedEnergy - minBoundEnergy)
				));
				
				// update confs list
				int numConfs = confs.size();
				for (int j=0; j<numConfs; j++) {
					int[] conf = confs.getConf(j);
					List<Integer> newRcs = splitWorld.getRcMap(splitPos).oldToNew(conf[splitPos]);
					if (newRcs.size() == 1) {
						conf[splitPos] = newRcs.get(0);
					} else {
						for (int k=0; k<newRcs.size(); k++) {
							int newRc = newRcs.get(k);
							if (k == 0) {
								
								// replace this conf
								conf[splitPos] = newRc;
								confs.removeEnergies(j);
								
							} else {
								
								// add new conf to the end
								int[] newConf = conf.clone();
								newConf[splitPos] = newRc;
								confs.add(newConf);
							}
						}
					}
				}
				
				// update energy and dof matrices
				System.out.println("updating energy and DoF matrices...");
				splitWorld.updateMatrices();
				
				// calculate energies for the new confs
				for (int j=0; j<confs.size(); j++) {
					int[] conf = confs.getConf(j);
					if (!confs.hasEnergies(j)) {
						confs.setEnergies(j, calcMinBoundEnergy(splitWorld.emat, conf), 0);
						/*
						System.out.println(String.format("new conf:   %30s  bound=%12f",
							Arrays.toString(conf), confs.getMinBound(j)
						));
						*/
					}
				}
				
				// find the min bound conf
				int bestConfi = confs.getMinBoundConfIndex();
				minBoundConf = confs.getConf(bestConfi);
				minBoundEnergy = confs.getMinBound(bestConfi);
				minimizedEnergy = splitWorld.minimize(minBoundConf);
				bestMinimizedEnergy = Math.min(bestMinimizedEnergy, minimizedEnergy);
				
				System.out.println(String.format("next conf:  %30s  bound=%12f  min=%12f  err=%12f",
					Arrays.toString(minBoundConf), minBoundEnergy, minimizedEnergy, (minimizedEnergy - minBoundEnergy)
				));
				
				// count the min bounds
				search.pruneMat = new PruningMatrix(search.confSpace, 1000);
				rcs = new RCs(search.pruneMat);
				order = new StaticScoreHMeanAStarOrder();
				hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 1, 0.0001);
				tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
				
				System.out.println("counting confs...");
				minBoundCount = 0;
				while (true) {
					ConfAStarNode leafNode = tree.nextLeafNode();
					if (leafNode == null) {
						break;
					}
					double bound = leafNode.getGScore();
					double betterBound = splitWorld.improveNode(leafNode).getGScore();
					if (betterBound < bestMinimizedEnergy) {
						minBoundCount++;
					}
					if (bound >= bestMinimizedEnergy) {
						break;
					}
				}
				System.out.println("confs to enumerate: " + minBoundCount);
				
				minBoundEnergies.add(minBoundEnergy);
				minimizedEnergies.add(minimizedEnergy);
				minBoundCounts.add(minBoundCount);
		
			/*
			} catch (Throwable t) {
				t.printStackTrace(System.out);
				break;
			}
			*/
		}
		
		// final report
		System.out.println("final improvements over " + numPartCrIters + " iterations:");
		System.out.println("min bound energies,minimized energies,num min bounds,split rots");
		for (int i=0; i<minBoundEnergies.size(); i++) {
			System.out.print(minBoundEnergies.get(i));
			System.out.print(",");
			System.out.print(minimizedEnergies.get(i));
			System.out.print(",");
			System.out.print(minBoundCounts.get(i));
			/*
			if (i < minBoundEnergies.size() - 1) {
				System.out.print(",");
				System.out.print(splitRotamers.get(i));
			}
			*/
			System.out.println();
		}
	}
	
	private static boolean areDofsInRange(double[] dofVals, RC splitRC) {
		assert (dofVals.length == splitRC.DOFs.size());
		for (int i=0; i<dofVals.length; i++) {
			double val = dofVals[i];
			double min = splitRC.DOFmin.get(i);
			double max = splitRC.DOFmax.get(i);
			if (val < min || val > max) {
				return false;
			}
		}
		return true;
	}

	private static int countBounds(ConfAStarTree tree, double minimizedEnergy) {
		int count = 0;
		while (true) {
			ConfAStarNode node = tree.nextLeafNode();
			if (node == null) {
				break;
			}
			count++;
			if (node.getGScore() >= minimizedEnergy) {
				break;
			}
		}
		return count;
	}

	private static ArrayList<Integer> getRcNums(Collection<RC> rcs) {
		ArrayList<Integer> nums = new ArrayList<>();
		for (RC rc : rcs) {
			nums.add(rc.RCIndex);
		}
		return nums;
	}
	
	private static double calcMinBoundEnergy(EnergyMatrix emat, int[] conf) {
		double energy = 0;
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			int rc1 = conf[pos1];
			energy += emat.getOneBody(pos1, rc1);
			for (int pos2=0; pos2<pos1; pos2++) {
				int rc2 = conf[pos2];
				energy += emat.getPairwise(pos1, rc1, pos2, rc2);
			}
		}
		return energy;
	}
	
	private static String dumpRc(RC rc, double[] dofBoundVals, int startDofIndex) {
		StringBuilder buf = new StringBuilder();
		for (int i=0; i<rc.DOFs.size(); i++) {
			double dofMinVal = rc.DOFs.get(i).getCurVal();
			double dofBoundVal = dofBoundVals[startDofIndex + i];
			double min = rc.DOFmin.get(i);
			double max = rc.DOFmax.get(i);
			String warnings = "";
			if (dofMinVal < min || dofMinVal > max) {
				warnings += " BADMIN";
			}
			if (dofBoundVal < min && dofBoundVal > max) {
				warnings += " BADBOUND";
			}
			buf.append(String.format("\tdof %d: min=%f, bound=%f, [%f,%f]%s\n", i, dofMinVal, dofBoundVal, min, max, warnings));
		}
		return buf.toString();
	}
	
	private static double getPosBoundEnergy(EnergyMatrix emat, int[] conf, int pos1) {
		
		int rc1 = conf[pos1];
		double energy = emat.getOneBody(pos1, rc1);
		
		for (int pos2=0; pos2<emat.getNumPos(); pos2++) {
			if (pos1 != pos2) {
				int rc2 = conf[pos2];
				energy += emat.getPairwise(pos1, rc1, pos2, rc2)/2;
			}
		}
		
		return energy;
	}
	
	private static double getPosMinEnergy(SimpleEnergyCalculator ecalc, int pos1) {
		
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
