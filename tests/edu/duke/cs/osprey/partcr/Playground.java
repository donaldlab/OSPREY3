package edu.duke.cs.osprey.partcr;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
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
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.EnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator;
import edu.duke.cs.osprey.ematrix.SimpleEnergyCalculator.Result;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.partcr.RCSplitter;
import edu.duke.cs.osprey.partcr.SmartRCSplitter;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.ObjectIO;

public class Playground extends TestBase {
	
	private static class DofMatrix extends TupleMatrixGeneric<double[]> {

		private static final long serialVersionUID = 7381812984847056950L;
		
		public DofMatrix(ConfSpace cSpace) {
			super (cSpace, 0, null);
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
		String flexRes = "38 39 40 41 42 43 44";
		ArrayList<String> flexResList = new ArrayList<>(Arrays.asList(flexRes.split(" ")));
		ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
		for (int i=0; i<flexResList.size(); i++) {
			allowedAAs.add(new ArrayList<>(Arrays.asList(aaNames.split(" "))));
		}
		boolean doMinimize = true;
		boolean addWt = true;
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
		
		// compute the energy matrix
		File ematFile = new File(String.format("/tmp/emat.partcr.dat"));
		if (ematFile.exists()) {
			System.out.println("\nReading energy matrix...");
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), true);
		}
		if (search.emat == null) {
			System.out.println("\nComputing energy matrix...");
			EnergyMatrixCalculator emCalc = new EnergyMatrixCalculator(search.confSpace, search.shellResidues, useERef, addResEntropy);
			emCalc.calcPEM();
			search.emat = emCalc.getEMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// compute the min dof values
		DofMatrix dofmat = null;
		File dofmatFile = new File(String.format("/tmp/dofmat.partcr.dat"));
		if (dofmatFile.exists()) {
			System.out.println("\nReading dof matrix...");
			dofmat = (DofMatrix)ObjectIO.readObject(dofmatFile.getAbsolutePath(), true);
		}
		if (dofmat == null) {
			System.out.println("\nComputing dof matrix...");
			dofmat = calcDofMat(search);
			ObjectIO.writeObject(dofmat, dofmatFile.getAbsolutePath());
		}
		
		// keep track of improvements over time
		List<Double> minBoundEnergies = new ArrayList<>();
		List<Double> minimizedEnergies = new ArrayList<>();
		List<String> splitRotamers = new ArrayList<>();
		List<Integer> minBoundCounts = new ArrayList<>();
		
		double bestEnergy = Double.POSITIVE_INFINITY;
		double pruningInterval = 10; // don't set this to 0, pruning is too good, only 6 conformations in the tree
		
		int numPartCrIters = 40;
		for (int i=0; i<numPartCrIters; i++) {
			
			// NEXTTIME: look at updating the pruning interval each iteration
			// read the paper! =)
			
			// do DEE
			boolean typeDep = false;
			double boundsThreshold = 100;
			int algoOption = 1;
			boolean useFlags = true;
			boolean useTriples = false;
			boolean useDacs = false;
			boolean useTupExp = false;
			double stericThreshold = 100;
			search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
			//
			new PruningControl(
				search, search.pruneMat.getPruningInterval(), typeDep, boundsThreshold,
				algoOption, useFlags, useTriples, useDacs, useEpic, useTupExp, stericThreshold
			).prune();
			//
			
			// get the min bound conformation
			System.out.println("Finding min bound GMEC among " + search.confSpace.getNumConformations().doubleValue() + " conformations...");
			RCs rcs = new RCs(search.pruneMat);
			AStarOrder order = new StaticScoreHMeanAStarOrder();
			AStarScorer hscorer = new MPLPPairwiseHScorer(new NodeUpdater(), search.emat, 4, 0.0001);
			ConfAStarTree tree = new ConfAStarTree(order, new PairwiseGScorer(search.emat), hscorer, rcs);
			tree.initProgress();
			int[] minBoundConf = tree.nextConf();
			
			// get the min bound and minimized energies
			System.out.println("Minimizing...");
			double minBoundEnergy = search.emat.getInternalEnergy(new RCTuple(minBoundConf));
			double minimizedEnergy = search.minimizedEnergy(minBoundConf);
			double boundError = minimizedEnergy - minBoundEnergy;
			System.out.println("min bound energy: " + minBoundEnergy);
			System.out.println("minimized energy: " + minimizedEnergy);
			System.out.println("Bound error:      " + boundError);
			
			// update best energy and pruning interval
			bestEnergy = Math.min(bestEnergy, minimizedEnergy);
			pruningInterval = bestEnergy - minBoundEnergy;
			System.out.println("interval: " + pruningInterval);
			
			/*
			// on the first and last iterations, enumerate lower bounds to check against PartCR
			if (i == 0 || i == numPartCrIters - 1) {
				System.out.println("lower bounds");
				System.out.println(minBoundEnergy);
				for (int j=1; j<numPartCrIters; j++) {
					System.out.println(tree.nextLeafNode().getGScore());
				}
			}
			*/
			
			// check how many bounds there are between the min bound and the min GMEC
			System.out.println("counting min bounds...");
			int numBounds = 0;
			while (true) {
				numBounds++;
				ConfAStarNode node = tree.nextLeafNode();
				if (node == null || node.getGScore() >= bestEnergy) {
					break;
				}
			}
			System.out.println("num min bounds: " + numBounds);
			minBoundCounts.add(numBounds);
			
			minBoundEnergies.add(minBoundEnergy);
			minimizedEnergies.add(minimizedEnergy);
			
			// on the last iteration, don't partition
			if (i == numPartCrIters - 1) {
				break;
			}
			System.out.println("PartCR iteration " + (i+1));

			// sort all the positions by errors
			TreeMap<Double,int[]> tuplesByError = new TreeMap<>();
			
			/* this heuristic for picking rotamers to split doesn't seem to work well
			   let's come up with something better
			for (int pos=0; pos<search.confSpace.numPos; pos++) {
				int rot = search.confSpace.posFlex.get(pos).RCs.get(minBoundConf[pos]).rotNum;
				double posBoundEnergy = getPosBoundEnergy(search.emat, minBoundConf, pos);
				double posMinEnergy = getPosMinEnergy(search, pos);
				double posErr = posMinEnergy - posBoundEnergy;
				System.out.println(String.format("pos=%2d, rot=%3d, bound=%10f, energy=%10f, err=%10f",
					pos, rot, posBoundEnergy, posMinEnergy, posErr
				));
				positionsByError.put(posErr, pos);
			}
			*/
			
			double minBoundEnergyCheck = 0;
			double minimizedEnergyCheck = 0;
			double boundErrorCheck = 0;
			EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
			for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
				Residue res1 = search.confSpace.posFlex.get(pos1).res;
				int rc1 = minBoundConf[pos1];
				
				// intra and shell energies
				double singleBoundEnergy = search.emat.getOneBody(pos1, rc1);
				
				double singleMinimizedEnergy = egen.singleResEnergy(res1).getEnergy();
				for (Residue shellResidue : search.shellResidues) {
					singleMinimizedEnergy += egen.resPairEnergy(res1, shellResidue).getEnergy();
				}
				
				double singleErr = singleMinimizedEnergy - singleBoundEnergy;
				//System.out.println(String.format("%5d: %f, %f, %f", pos1, singleBoundEnergy, singleMinimizedEnergy, singleErr));
				minBoundEnergyCheck += singleBoundEnergy;
				minimizedEnergyCheck += singleMinimizedEnergy;
				boundErrorCheck += singleErr;
				
				/* check DOFs
				RC rc1Obj = search.confSpace.posFlex.get(pos1).RCs.get(rc1);
				for (int j=0; j<rc1Obj.DOFs.size(); j++) {
					double dofMinVal = rc1Obj.DOFs.get(j).getCurVal();
					double dofBoundVal = dofmat.getOneBody(pos1, rc1)[j];
					System.out.println(String.format("\tdof %d: min=%f, bound=%f, [%f,%f]", j, dofMinVal, dofBoundVal, rc1Obj.DOFmin.get(j), rc1Obj.DOFmax.get(j)));
				}
				*/
				
				tuplesByError.put(singleErr, new int[] { pos1 });
				
				// pairwise energies
				for (int pos2=0; pos2<pos1; pos2++) {
					Residue res2 = search.confSpace.posFlex.get(pos2).res;
					int rc2 = minBoundConf[pos2];
					
					double pairwiseBoundEnergy = search.emat.getPairwise(pos1, rc1, pos2, rc2);
					double pairwiseMinimizedEnergy = egen.resPairEnergy(res1, res2).getEnergy();
					
					double pairwiseErr = pairwiseMinimizedEnergy - pairwiseBoundEnergy;
					//System.out.println(String.format("%2d,%2d: %f, %f, %f", pos1, pos2, pairwiseBoundEnergy, pairwiseMinimizedEnergy, pairwiseErr));
					minBoundEnergyCheck += pairwiseBoundEnergy;
					minimizedEnergyCheck += pairwiseMinimizedEnergy;
					boundErrorCheck += pairwiseErr;
					
					/* check DOFs
					RC rc2Obj = search.confSpace.posFlex.get(pos2).RCs.get(rc2);
					for (int j=0; j<rc1Obj.DOFs.size(); j++) {
						double dofMinVal = rc1Obj.DOFs.get(j).getCurVal();
						double dofBoundVal = dofmat.getOneBody(pos1, rc1)[j];
						System.out.println(String.format("\tdof %d: min=%f, bound=%f, [%f,%f]", j, dofMinVal, dofBoundVal, rc1Obj.DOFmin.get(j), rc1Obj.DOFmax.get(j)));
					}
					for (int j=0; j<rc2Obj.DOFs.size(); j++) {
						double dofMinVal = rc2Obj.DOFs.get(j).getCurVal();
						double dofBoundVal = dofmat.getOneBody(pos2, rc2)[j];
						System.out.println(String.format("\tdof %d: min=%f, bound=%f, [%f,%f]", j, dofMinVal, dofBoundVal, rc2Obj.DOFmin.get(j), rc2Obj.DOFmax.get(j)));
					}
					*/
					
					tuplesByError.put(pairwiseErr, new int[] { pos1, pos2 });
				}
			}
			//System.out.println("min bound energy err: " + (minBoundEnergy - minBoundEnergyCheck));
			//System.out.println("minimized energy err: " + (minimizedEnergy - minimizedEnergyCheck));
			//System.out.println("Bound error err:      " + (boundError - boundErrorCheck));
			
			// sort tuples by error
			List<int[]> tuplesInOrder = new ArrayList<>(tuplesByError.values());
			Collections.reverse(tuplesInOrder);
			
			// NOTE: the new A* implementation is super fast
			// which means we can bombard it with tons of new RCs and it won't slow down much
			// which means we should really aggressive with RC partitioning
			// since A* time is basically zero compared to energy calculation time
			int numTuplesToSplit = 1;
			
			// prep the mappings from new rcs to old rcs, to speed up energy matrix calculation
			List<List<Integer>> rcMaps = new ArrayList<>();
			for (int j=0; j<search.confSpace.numPos; j++) {
				rcMaps.add(null);
			}
			
			String splitRotsDesc = "";
			
			// NEXTTIME: make a smarter splitter
			// explicitly try to put a partition between the
			// conf picked by the energy timer minimizer
			// and the conf picked by the whole molecule minimization
			
			// split the tuples with the worst errors
			int numSplitTuples = 0;
			for (int[] tuple : tuplesInOrder) {
				
				// get a RC splitter
				RCSplitter splitter;
				if (tuple.length == 1) {
					splitter = new SmartRCSplitter(dofmat.getOneBody(tuple[0], minBoundConf[tuple[0]]));
				} else if (tuple.length == 2) {
					splitter = new SmartRCSplitter(dofmat.getPairwise(tuple[0], minBoundConf[tuple[0]], tuple[1], minBoundConf[tuple[1]]));
				} else {
					throw new Error("bad tuple length");
				}
				
				// split each pos in the tuple
				boolean splitSomething = false;
				for (int pos : tuple) {
					
					// get the rc at this pos
					PositionConfSpace posConfSpace = search.confSpace.posFlex.get(pos);
					int rc = minBoundConf[pos];
					RC rcObj = posConfSpace.RCs.get(rc);
					int rotNum = rcObj.rotNum;
					
					// will the splitter split this one?
					if (!splitter.willSplit(pos, rcObj)) {
						System.out.println(String.format("skipping rot %d-%d", pos, rotNum));
						continue;
					}
					
					// TEMP
					System.out.println(String.format("splitting %d-%d", pos, rotNum));
					splitRotsDesc += String.format("%d-%d-%d ", pos, rotNum, rc);
					
					// partition the rc and update the conf space
					List<RC> splitRCs = splitter.split(pos, rcObj);
					List<Integer> rcMap = posConfSpace.replaceRC(rcObj, splitRCs);
					rcMaps.set(pos, rcMap);
					
					splitSomething = true;
				}
				
				if (splitSomething) {
					numSplitTuples++;
				}
				
				if (numSplitTuples >= numTuplesToSplit) {
					break;
				}
			}
			if (numSplitTuples == 0) {
				
				// ran out of things to split! bail
				break;
			}
			
			splitRotamers.add(splitRotsDesc);
			
			// update energy and dof matrices
			System.out.println("updating energy and DoF matrices...");
			dofmat = updateMatrices(search, dofmat, rcMaps);
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
			if (i < minBoundEnergies.size() - 1) {
				System.out.print(",");
				System.out.print(splitRotamers.get(i));
			}
			System.out.println();
		}
	}
	
	private static DofMatrix calcDofMat(SearchProblem search) {
		
		DofMatrix dofmat = new DofMatrix(search.confSpace);
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(
			EnvironmentVars.curEFcnGenerator,
			search.confSpace,
			search.shellResidues
		);
		
    	for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
			System.out.println(String.format("Minimizing DOFs for residue %2d", pos1));
			for (int rc1=0; rc1<dofmat.getNumConfAtPos(pos1); rc1++) {

				// one-body
				dofmat.setOneBody(pos1, rc1, ecalc.calcSingle(pos1, rc1).getDofValues());
				
				// pairwise
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<dofmat.getNumConfAtPos(pos2); rc2++) {
						dofmat.setPairwise(pos1, rc1, pos2, rc2, ecalc.calcPair(pos1, rc1, pos2, rc2).getDofValues());
					}
				}
			}
    	}
		
		return dofmat;
	}

	private static DofMatrix updateMatrices(SearchProblem search, DofMatrix oldDofmat, List<List<Integer>> rcMaps) {
		
		EnergyMatrix oldEmat = search.emat;
		
		// copy as many energies as possible from the old matrix
		// to a new matrix that is sized for the new conf space (with the added rotamers)
		EnergyMatrix newEmat = new EnergyMatrix(search.confSpace, oldEmat.getPruningInterval());
		DofMatrix newDofmat = new DofMatrix(search.confSpace);
		
		SimpleEnergyCalculator ecalc = new SimpleEnergyCalculator(
			EnvironmentVars.curEFcnGenerator,
			search.confSpace,
			search.shellResidues
		);
		SimpleEnergyCalculator.ShellDistribution dist = SimpleEnergyCalculator.ShellDistribution.AllOnSingles;
		
    	for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
			for (int rc1=0; rc1<newEmat.getNumConfAtPos(pos1); rc1++) {

				// get the old rc1, if any
				Integer oldRc1 = rc1;
				if (rcMaps.get(pos1) != null) {
					oldRc1 = rcMaps.get(pos1).get(rc1);
				}
				
				// one-body
				if (oldRc1 == null) {
					Result result = ecalc.calcSingle(pos1, rc1);
					newDofmat.setOneBody(pos1, rc1, result.getDofValues());
					newEmat.setOneBody(pos1, rc1, result.getEnergy());
				} else {
					newDofmat.setOneBody(pos1, rc1, oldDofmat.getOneBody(pos1, oldRc1));
					newEmat.setOneBody(pos1, rc1, oldEmat.getOneBody(pos1, oldRc1));
				}
				
				// pairwise
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<newEmat.getNumConfAtPos(pos2); rc2++) {
						
						// get the old rc2, if any
						Integer oldRc2 = rc2;
						if (rcMaps.get(pos2) != null) {
							oldRc2 = rcMaps.get(pos2).get(rc2);
						}
						
						if (oldRc1 == null || oldRc2 == null) {
							Result result = ecalc.calcPair(pos1, rc1, pos2, rc2);
							newDofmat.setPairwise(pos1, rc1, pos2, rc2, result.getDofValues());
							newEmat.setPairwise(pos1, rc1, pos2, rc2, result.getEnergy());
						} else {
							newDofmat.setPairwise(pos1, rc1, pos2, rc2, oldDofmat.getPairwise(pos1, oldRc1, pos2, oldRc2));
							newEmat.setPairwise(pos1, rc1, pos2, rc2, oldEmat.getPairwise(pos1, oldRc1, pos2, oldRc2));
						}
					}
				}
			}	
    	}
    	
    	search.emat = newEmat;
    	return newDofmat;
	}
	
	private static String dumpRc(RC rc, String title) {
		StringBuilder buf = new StringBuilder();
		buf.append(title + ":");
		for (int i=0; i<rc.DOFs.size(); i++) {
			buf.append(String.format("\t%s: [%f,%f]",
				rc.DOFs.get(i).getClass().getSimpleName(),
				rc.DOFmin.get(i),
				rc.DOFmax.get(i)
			));
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
	
	private static double getPosMinEnergy(SearchProblem search, int pos1) {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		double energy = 0;
		
		// single energy
		Residue res1 = search.confSpace.posFlex.get(pos1).res;
		energy += egen.singleResEnergy(res1).getEnergy();
		
		// single to shell energy
		for (Residue shellResidue : search.shellResidues) {
			energy += egen.resPairEnergy(res1, shellResidue).getEnergy();
		}
		
		// pairwise energy
		for (PositionConfSpace flexPos : search.confSpace.posFlex) {
			Residue res2 = flexPos.res;
			if (res1 != res2) {
				energy += egen.resPairEnergy(res1, res2).getEnergy()/2;
			}
		}
		
		return energy;
	}
}
