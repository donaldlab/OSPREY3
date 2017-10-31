package edu.duke.cs.osprey.ematrix;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Progress;

public class SimplerEnergyMatrixCalculator {
	
	// NOTE: don't use GPUs on energy matrices, it's too slow
	// always use the CPU
	
	public static class Builder {
		
		/**
		 * How conformation energies should be calculated.
		 */
		private ConfEnergyCalculator confEcalc;
		
		/**
		 * Path to file where energy matrix should be saved between computations.
		 * 
		 * @note Energy matrix computation can take a long time, but often the results
		 * can be reused between computations. Use a cache file to skip energy matrix
		 * computation on the next Osprey run if the energy matrix has already been
		 * computed once before.
		 * 
		 * @warning If design settings are changed between runs, Osprey will make
		 * some effort to detect that the energy matrix cache is out-of-date and compute a
		 * new energy matrix instead of usng the cached, incorrect one. Osprey might not detect
		 * all design changes though, and incorrectly reuse a cached energy matrix, so it
		 * is best to manually delete the entry matrix cache file after changing design settings.
		 */
		private File cacheFile = null;
		
		public Builder(SimpleConfSpace confSpace, EnergyCalculator ecalc) {
			this(new ConfEnergyCalculator.Builder(confSpace, ecalc).build());
		}
		
		public Builder(ConfEnergyCalculator confEcalc) {
			this.confEcalc = confEcalc;
		}
		
		public Builder setCacheFile(File val) {
			cacheFile = val;
			return this;
		}
		
		public SimplerEnergyMatrixCalculator build() {
			return new SimplerEnergyMatrixCalculator(confEcalc, cacheFile);
		}
	}
	
	public final ConfEnergyCalculator confEcalc;
	public final File cacheFile;

	private SimplerEnergyMatrixCalculator(ConfEnergyCalculator confEcalc, File cacheFile) {
		this.confEcalc = confEcalc;
		this.cacheFile = cacheFile;
	}
	
	/**
	 * Computes a matrix of energies between pairs of residue conformations to be used by A* search.
	 */
	public EnergyMatrix calcEnergyMatrix() {
		
		if (cacheFile != null) {
			return ObjectIO.readOrMake(
				cacheFile,
				EnergyMatrix.class,
				"energy matrix",
				(emat) -> emat.matches(confEcalc.confSpace),
				(context) -> reallyCalcEnergyMatrix()
			);
		} else {
			return reallyCalcEnergyMatrix();
		}
	}
	
	private EnergyMatrix reallyCalcEnergyMatrix() {
		
		// allocate the new matrix
		EnergyMatrix emat = new EnergyMatrix(confEcalc.confSpace);
		
		// count how much work there is to do (roughly based on number of residue pairs)
		final int singleCost = confEcalc.makeSingleInters(0, 0).size();
		final int pairCost = confEcalc.makePairInters(0, 0, 0, 0).size();
		Progress progress = new Progress(confEcalc.confSpace.getNumResConfs()*singleCost + confEcalc.confSpace.getNumResConfPairs()*pairCost);
		
		// some fragments can be big and some can be small
		// try minimize thread sync overhead by not sending a bunch of small fragments in all separate tasks
		// ie, try to batch fragments together
		class Batch {
			
			List<RCTuple> fragments = new ArrayList<>();
			int cost = 0;
			
			void addSingle(int pos, int rc) {
				fragments.add(new RCTuple(pos, rc));
				cost += singleCost;
			}
			
			void addPair(int pos1, int rc1, int pos2, int rc2) {
				fragments.add(new RCTuple(pos1, rc1, pos2, rc2));
				cost += pairCost;
			}
			
			void submitTask() {
				confEcalc.tasks.submit(
					() -> {
						
						// calculate all the fragment energies
						List<Double> energies = new ArrayList<>();
						for (RCTuple frag : fragments) {
							if (frag.size() == 1) {
								energies.add(confEcalc.calcSingleEnergy(frag).energy);
							} else {
								energies.add(confEcalc.calcPairEnergy(frag).energy);
							}
						}
						
						return energies;
					},
					(List<Double> energies) -> {
						
						// update the energy matrix
						for (int i=0; i<fragments.size(); i++) {
							RCTuple frag = fragments.get(i);
							if (frag.size() == 1) {
								emat.setOneBody(frag.pos.get(0), frag.RCs.get(0), energies.get(i));
							} else {
								emat.setPairwise(frag.pos.get(0), frag.RCs.get(0), frag.pos.get(1), frag.RCs.get(1), energies.get(i));
							}
						}
						
						progress.incrementProgress(cost);
					}
				);
			}
		}
		
		final int CostThreshold = 100;
		
		class Batcher {
			
			Batch batch = null;
			
			Batch getBatch() {
				if (batch == null) {
					batch = new Batch();
				}
				return batch;
			}
			
			void submitIfFull() {
				if (batch != null && batch.cost >= CostThreshold) {
					submit();
				}
			}
			
			void submit() {
				if (batch != null) {
					batch.submitTask();
					batch = null;
				}
			}
		}
		Batcher batcher = new Batcher();
		
		// batch all the singles and pairs
		System.out.println("Calculating energy matrix with " + (confEcalc.confSpace.getNumResConfs() + confEcalc.confSpace.getNumResConfPairs()) + " entries...");
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
				
				// singles
				batcher.getBatch().addSingle(pos1, rc1);
				batcher.submitIfFull();
				
				// pairs
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {
						batcher.getBatch().addPair(pos1, rc1, pos2, rc2);
						batcher.submitIfFull();
					}
				}
			}
		}
		
		batcher.submit();
		confEcalc.tasks.waitForFinish();
		
		return emat;
	}
	
	
	/**
	 * Calculates a reference energy for each residue position and residue type
	 * based on the minimum energy of all residue conformations at that position
	 * and residue type.
	 */
	public SimpleReferenceEnergies calcReferenceEnergies() {
		
		SimpleReferenceEnergies eref = new SimpleReferenceEnergies();
		
		// send all the tasks
		Progress progress = new Progress(confEcalc.confSpace.getNumResConfs());
		System.out.println("Calculating reference energies for " + progress.getTotalWork() + " residue confs...");
		for (SimpleConfSpace.Position pos : confEcalc.confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
			
				String resType = rc.template.name;
				RCTuple frag = new RCTuple(pos.index, rc.index);
				confEcalc.calcEnergyAsync(
					frag,
					ResInterGen.of(confEcalc.confSpace).addIntra(pos.index).make(),
					(EnergyCalculator.EnergiedParametricMolecule epmol) -> {
						
						// keep the min energy for each pos,resType
						Double e = eref.get(frag.pos.get(0), resType);
						if (e == null || epmol.energy < e) {
							e = epmol.energy;
						}
						eref.set(frag.pos.get(0), resType, e);
						progress.incrementProgress();
					}
				);
			}
		}
		
		confEcalc.tasks.waitForFinish();
		
		return eref;
	}
}
