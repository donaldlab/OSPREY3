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

package edu.duke.cs.osprey.ematrix;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tools.Progress;

import static edu.duke.cs.osprey.tools.Log.log;

public class SimplerEnergyMatrixCalculator {
	
	// NOTE: don't use GPUs on energy matrices, it's too slow
	// always use the CPU
	// (until we implement efficient batching in the GPU energy calculator that is...)
	// unless you're using AllOnPairs energy partition, then the GPU is pretty fast

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

		/**
		 * Compute energy corrections for all triples whose constituent single and pair energies
		 * are below the given threshold. ie. ignore triples with clashes.
		 */
		private Double tripleCorrectionThreshold = null;

		/**
		 * Compute energy corrections for all quads whose constituent single and pair energies
		 * are below the given threshold. ie. ignore quads with clashes.
		 */
		private Double quadCorrectionThreshold = null;
		
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

		public Builder setTripleCorrectionThreshold(Double val) {
			tripleCorrectionThreshold = val;
			return this;
		}

		public Builder setQuadCorrectionThreshold(Double val) {
			quadCorrectionThreshold = val;
			return this;
		}
		
		public SimplerEnergyMatrixCalculator build() {
			return new SimplerEnergyMatrixCalculator(confEcalc, cacheFile, tripleCorrectionThreshold, quadCorrectionThreshold);
		}
	}
	
	public final ConfEnergyCalculator confEcalc;
	public final File cacheFile;
	public final Double tripleCorrectionThreshold;
	public final Double quadCorrectionThreshold;

	private SimplerEnergyMatrixCalculator(ConfEnergyCalculator confEcalc, File cacheFile, Double tripleCorrectionThreshold, Double quadCorrectionThreshold) {
		this.confEcalc = confEcalc;
		this.cacheFile = cacheFile;
		this.tripleCorrectionThreshold = tripleCorrectionThreshold;
		this.quadCorrectionThreshold = quadCorrectionThreshold;
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
		Progress progress = new Progress(
			confEcalc.confSpace.getNumResConfs()*singleCost
				+ confEcalc.confSpace.getNumResConfPairs()*pairCost
		);
		
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

							double energy;

							// are there any RCs are from two different backbone states that can't connect?
							if (isParametricallyIncompatible(frag)) {

								// yup, give this frag an infinite energy so we never choose it
								energy = Double.POSITIVE_INFINITY;

							} else {

								// nope, calculate the usual fragment energy
								switch (frag.size()) {
									case 1: {
										energy = confEcalc.calcSingleEnergy(frag).energy;
									} break;
									case 2: {
										energy = confEcalc.calcPairEnergy(frag).energy;
									} break;
									default: {
										energy = confEcalc.calcEnergy(frag).energy;
									}
								}
							}

							energies.add(energy);
						}
						
						return energies;
					},
					(List<Double> energies) -> {
						
						// update the energy matrix
						for (int i=0; i<fragments.size(); i++) {
							RCTuple frag = fragments.get(i);
							if (frag.size() == 1) {
								emat.setOneBody(frag.pos.get(0), frag.RCs.get(0), energies.get(i));
							} else if (frag.size() == 2) {
								emat.setPairwise(frag.pos.get(0), frag.RCs.get(0), frag.pos.get(1), frag.RCs.get(1), energies.get(i));
							} else {
								emat.setTuple(frag, energies.get(i));
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
		
		// convert the workload into tasks for the task executor
		log("Calculating energy matrix with %d entries",
			confEcalc.confSpace.getNumResConfs()
			+ confEcalc.confSpace.getNumResConfPairs()
		);
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
				
				// single
				batcher.getBatch().addSingle(pos1, rc1);
				batcher.submitIfFull();
				
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {

						// pair
						batcher.getBatch().addPair(pos1, rc1, pos2, rc2);
						batcher.submitIfFull();
					}
				}
			}
		}
		
		batcher.submit();
		confEcalc.tasks.waitForFinish();

		// calc corrections if needed (but only use the highest-order corrections chosen)
		if (quadCorrectionThreshold != null) {
			calcQuadCorrections(emat);
		} else if (tripleCorrectionThreshold != null) {
			calcTripleCorrections(emat);
		}

		return emat;
	}

	// TODO: improve progress bar performance by pre-counting the tuples that pass the threshold

	private void calcTripleCorrections(EnergyMatrix emat) {

		Progress progress = new Progress(confEcalc.confSpace.getNumResConfTriples());
		log("calculating triple corrections for up to %d triples", progress.getTotalWork());
		int[] numCorrections = { 0 };

		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {
						for (int pos3=0; pos3<pos2; pos3++) {
							for (int rc3=0; rc3<emat.getNumConfAtPos(pos3); rc3++) {

								// if any of the components are too high, skip this triple
								if (emat.getOneBody(pos1, rc1) > tripleCorrectionThreshold
									|| emat.getOneBody(pos2, rc2) > tripleCorrectionThreshold
									|| emat.getOneBody(pos3, rc3) > tripleCorrectionThreshold
									|| emat.getPairwise(pos1, rc1, pos2, rc2) > tripleCorrectionThreshold
									|| emat.getPairwise(pos1, rc1, pos3, rc3) > tripleCorrectionThreshold
									|| emat.getPairwise(pos2, rc2, pos3, rc3) > tripleCorrectionThreshold) {

									synchronized (progress) {
										progress.incrementProgress();
									}
									continue;
								}

								final RCTuple triple = new RCTuple(pos3, rc3, pos2, rc2, pos1, rc1);

								// check the triple for parametric incompatibilities
								if (isParametricallyIncompatible(triple)) {
									synchronized (progress) {
										progress.incrementProgress();
									}
									continue;
								}

								ResidueInteractions inters = confEcalc.makeTripleCorrectionInters(pos1, rc1, pos2, rc2, pos3, rc3);
								double tripleEnergyOffset = confEcalc.epart.offsetTripleEnergy(pos1, rc1, pos2, rc2, pos3, rc3, emat);

								// calc the energy
								confEcalc.tasks.submit(
									() -> confEcalc.calcEnergy(triple, inters).energy,
									(tripleEnergy) -> {

										// convert the triple energy into a correction
										double correction = tripleEnergy - tripleEnergyOffset;

										// save the correction only if it's an improvement
										if (correction > 0) {
											emat.setTuple(triple, correction);
											numCorrections[0]++;
										}

										synchronized (progress) {
											progress.incrementProgress();
										}
									}
								);
							}
						}
					}
				}
			}
		}

		confEcalc.tasks.waitForFinish();

		log("calculated %d/%d useful triple corrections", numCorrections[0], progress.getTotalWork());
	}

	private void calcQuadCorrections(EnergyMatrix emat) {

		Progress progress = new Progress(confEcalc.confSpace.getNumResConfQuads());
		log("calculating quad corrections for up to %d quads", progress.getTotalWork());
		int[] numCorrections = { 0 };

		// TODO: this indentation is ridiculous...
		for (int pos1=0; pos1<emat.getNumPos(); pos1++) {
			for (int rc1=0; rc1<emat.getNumConfAtPos(pos1); rc1++) {
				for (int pos2=0; pos2<pos1; pos2++) {
					for (int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++) {
						for (int pos3=0; pos3<pos2; pos3++) {
							for (int rc3=0; rc3<emat.getNumConfAtPos(pos3); rc3++) {
								for (int pos4=0; pos4<pos3; pos4++) {
									for (int rc4=0; rc4<emat.getNumConfAtPos(pos4); rc4++) {

										// if any of the components are too high, skip this quad
										if (emat.getOneBody(pos1, rc1) > quadCorrectionThreshold
											|| emat.getOneBody(pos2, rc2) > quadCorrectionThreshold
											|| emat.getOneBody(pos3, rc3) > quadCorrectionThreshold
											|| emat.getOneBody(pos4, rc4) > quadCorrectionThreshold
											|| emat.getPairwise(pos1, rc1, pos2, rc2) > quadCorrectionThreshold
											|| emat.getPairwise(pos1, rc1, pos3, rc3) > quadCorrectionThreshold
											|| emat.getPairwise(pos1, rc1, pos4, rc4) > quadCorrectionThreshold
											|| emat.getPairwise(pos2, rc2, pos3, rc3) > quadCorrectionThreshold
											|| emat.getPairwise(pos2, rc2, pos4, rc4) > quadCorrectionThreshold
											|| emat.getPairwise(pos3, rc3, pos4, rc4) > quadCorrectionThreshold
										) {

											synchronized (progress) {
												progress.incrementProgress();
											}
											continue;
										}

										final RCTuple quad = new RCTuple(pos4, rc4, pos3, rc3, pos2, rc2, pos1, rc1);

										// check the quad for parametric incompatibilities
										if (isParametricallyIncompatible(quad)) {
											synchronized (progress) {
												progress.incrementProgress();
											}
											continue;
										}

										ResidueInteractions inters = confEcalc.makeQuadCorrectionInters(pos1, rc1, pos2, rc2, pos3, rc3, pos4, rc4);
										double quadEnergyOffset = confEcalc.epart.offsetQuadEnergy(pos1, rc1, pos2, rc2, pos3, rc3, pos4, rc4, emat);

										// calc the energy
										confEcalc.tasks.submit(
											() -> confEcalc.calcEnergy(quad, inters).energy,
											(quadEnergy) -> {

												// convert the quad energy into a correction
												double correction = quadEnergy - quadEnergyOffset;

												// save the correction only if it's an improvement
												if (correction > 0) {
													emat.setTuple(quad, correction);
													numCorrections[0]++;
												}

												synchronized (progress) {
													progress.incrementProgress();
												}
											}
										);
									}
								}
							}
						}
					}
				}
			}
		}

		confEcalc.tasks.waitForFinish();

		log("calculated %d/%d useful quad corrections", numCorrections[0], progress.getTotalWork());
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

	private boolean isParametricallyIncompatible(RCTuple tuple) {
		for (int i1=0; i1<tuple.size(); i1++) {
			SimpleConfSpace.ResidueConf rc1 = getRC(tuple, i1);
			for (int i2=0; i2<i1; i2++) {
				SimpleConfSpace.ResidueConf rc2 = getRC(tuple, i2);
				if (!isPairParametricallyCompatible(rc1, rc2)) {
					return true;
				}
			}
		}
		return false;
	}

	private SimpleConfSpace.ResidueConf getRC(RCTuple tuple, int index) {
		return confEcalc.confSpace.positions.get(tuple.pos.get(index)).resConfs.get(tuple.RCs.get(index));
	}

	private boolean isPairParametricallyCompatible(SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.ResidueConf rc2) {
		for(String dofName : rc1.dofBounds.keySet()){
			if(rc2.dofBounds.containsKey(dofName)){
				//shared DOF between the RCs; make sure the interval matches
				double[] interval1 = rc1.dofBounds.get(dofName);
				double[] interval2 = rc2.dofBounds.get(dofName);
				for(int a=0; a<2; a++){
					if( Math.abs(interval1[a] - interval2[a]) > 1e-8 ){
						return false;
					}
				}
			}
		}
		return true;//found no incompatibilities
	}
}
