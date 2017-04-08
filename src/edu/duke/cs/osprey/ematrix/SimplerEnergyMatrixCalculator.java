package edu.duke.cs.osprey.ematrix;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.FragmentEnergyCalculator;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Progress;

public class SimplerEnergyMatrixCalculator {
	
	// NOTE: don't use GPUs on energy matrices, it's too slow
	// always use the CPU
	
	public static class Builder {
		
		/**
		 * The conformation space containing the strands to be designed.
		 * 
		 * If the strands are configured with continuous flexibility, the energy matrix will
		 * minimize residue conformation pairs before computing energies.
		 */
		private SimpleConfSpace confSpace;
		
		/**
		 * How conformation energies should be calculated.
		 */
		private FragmentEnergyCalculator.Async ecalc;
		
		/**
		 * How energies should be paritioned among single and pair fragments.
		 */
		private EnergyPartition epart = EnergyPartition.Traditional;
		
		private SimpleReferenceEnergies eref = null;
		
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
		
		public Builder(SimpleConfSpace confSpace, FragmentEnergyCalculator.Async ecalc) {
			this.confSpace = confSpace;
			this.ecalc = ecalc;
		}
		
		public Builder setCacheFile(File val) {
			cacheFile = val;
			return this;
		}
		
		public Builder setReferenceEnergies(SimpleReferenceEnergies val) {
			eref = val;
			return this;
		}
		
		public SimplerEnergyMatrixCalculator build() {
			return new SimplerEnergyMatrixCalculator(confSpace, ecalc, epart, eref, cacheFile);
		}
	}
	
	private final SimpleConfSpace confSpace;
	private final FragmentEnergyCalculator.Async ecalc;
	private final EnergyPartition epart;
	private final SimpleReferenceEnergies eref;
	private final File cacheFile;

	private SimplerEnergyMatrixCalculator(SimpleConfSpace confSpace, FragmentEnergyCalculator.Async ecalc, EnergyPartition epart, SimpleReferenceEnergies eref, File cacheFile) {
		this.confSpace = confSpace;
		this.ecalc = ecalc;
		this.epart = epart;
		this.eref = eref;
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
				(emat) -> emat.matches(confSpace),
				(context) -> reallyCalcEnergyMatrix()
			);
		} else {
			return reallyCalcEnergyMatrix();
		}
	}
	
	private EnergyMatrix reallyCalcEnergyMatrix() {
		
		// allocate the new matrix
		EnergyMatrix emat = new EnergyMatrix(confSpace);
		
		// count how much work there is to do (roughly based on number of residue pairs)
		// TODO: base cost on energy partition
		final int singleCost = confSpace.shellResNumbers.size() + 1;
		final int pairCost = 1;
		Progress progress = new Progress(confSpace.getNumResConfs()*singleCost + confSpace.getNumResConfPairs()*pairCost);
		
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
				ecalc.getTasks().submit(
					() -> {
						
						// calculate all the fragment energies
						List<Double> energies = new ArrayList<>();
						for (RCTuple frag : fragments) {
							if (frag.size() == 1) {
								energies.add(calcSingle(frag));
							} else {
								energies.add(calcPair(frag));
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
		System.out.println("Calculating energy matrix with " + progress.getTotalWork() + " entries...");
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
		ecalc.getTasks().waitForFinish();
		
		return emat;
	}
	
	public double calcSingle(int pos, int rc) {
		return calcSingle(new RCTuple(pos, rc));
	}
	
	public double calcSingle(RCTuple frag) {
		return ecalc.calcEnergy(frag, epart.makeSingleInters(confSpace, eref, frag.pos.get(0), frag.RCs.get(0)));
	}
	
	public double calcPair(int pos1, int rc1, int pos2, int rc2) {
		return calcPair(new RCTuple(pos1, rc1, pos2, rc2));
	}
	
	public double calcPair(RCTuple frag) {
		return ecalc.calcEnergy(frag, epart.makePairInters(confSpace, eref, frag.pos.get(0), frag.RCs.get(0), frag.pos.get(1), frag.RCs.get(1)));
	}
	
	/**
	 * Calculates a reference energy for each residue position and residue type
	 * based on the minimum energy of all residue conformations at that position
	 * and residue type.
	 */
	public SimpleReferenceEnergies calcReferenceEnergies() {
		
		SimpleReferenceEnergies eref = new SimpleReferenceEnergies();
		
		// send all the tasks
		Progress progress = new Progress(confSpace.getNumResConfs());
		System.out.println("Calculating reference energies for " + progress.getTotalWork() + " residue confs...");
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {
			
				String resType = rc.template.name;
				RCTuple frag = new RCTuple(pos.index, rc.index);
				ecalc.calcEnergyAsync(
					frag,
					ResInterGen.of(confSpace).addIntra(pos.index).make(),
					(Double energy) -> {
						
						// keep the min energy for each pos,resType
						Double e = eref.get(frag.pos.get(0), resType);
						if (e == null || energy < e) {
							e = energy;
						}
						eref.set(frag.pos.get(0), resType, e);
						progress.incrementProgress();
					}
				);
			}
		}
		
		ecalc.getTasks().waitForFinish();
		
		return eref;
	}
}
