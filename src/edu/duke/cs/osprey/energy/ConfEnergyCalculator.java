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

package edu.duke.cs.osprey.energy;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.Supplier;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.energy.approximation.ApproximatorMatrix;
import edu.duke.cs.osprey.energy.approximation.ResidueInteractionsApproximator;
import edu.duke.cs.osprey.minimization.MoleculeObjectiveFunction;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.parallelism.TaskExecutor.TaskListener;
import edu.duke.cs.osprey.tools.Progress;
import edu.duke.cs.osprey.tools.TimeTools;


/**
 * Calculate energy for molecules created from conformation spaces.
 * 
 * Provides support for applying conformation energy modifications,
 * such as reference energies, residue entropies, and energy partitions. 
 */
public class ConfEnergyCalculator {
	
	public static class Builder {
		
		private SimpleConfSpace confSpace;
		private EnergyCalculator ecalc;
		
		/**
		 * How energies should be partitioned among single and pair fragments.
		 */
		private EnergyPartition epart = EnergyPartition.Traditional;
		
		private SimpleReferenceEnergies eref = null;
		private boolean addResEntropy = false;

		/** The approximator matrix */
		private ApproximatorMatrix amat = null;

		/** How much error can be tolerated for the energy of a conformation? */
		private double approximationErrorBudget = 1e-2;
		
		public Builder(SimpleConfSpace confSpace, EnergyCalculator ecalc) {
			this.confSpace  = confSpace;
			this.ecalc = ecalc;
		}
		
		public Builder setEnergyPartition(EnergyPartition val) {
			this.epart = val;
			return this;
		}
		
		public Builder setReferenceEnergies(SimpleReferenceEnergies val) {
			this.eref = val;
			return this;
		}
		
		public Builder addResEntropy(boolean val) {
			this.addResEntropy = val;
			return this;
		}

		public Builder setApproximatorMatrix(ApproximatorMatrix val) {
			this.amat = val;
			return this;
		}

		public Builder setApproximationErrorBudget(double val) {
			this.approximationErrorBudget = val;
			return this;
		}
		
		public ConfEnergyCalculator build() {
			return new ConfEnergyCalculator(confSpace, ecalc, epart, eref, addResEntropy, amat, approximationErrorBudget);
		}
	}
	
	public final SimpleConfSpace confSpace;
	public final EnergyCalculator ecalc;
	public final EnergyPartition epart;
	public final SimpleReferenceEnergies eref;
	public final boolean addResEntropy;
	public final ApproximatorMatrix amat;
	public final double approximationErrorBudget;

	public final TaskExecutor tasks;

	protected final AtomicLong numCalculations = new AtomicLong(0L);
	protected final AtomicLong numConfDBReads = new AtomicLong(0L);

	protected ConfEnergyCalculator(SimpleConfSpace confSpace, TaskExecutor tasks) {
		this(confSpace, null, null, null, false, null, Double.NaN);
	}

	protected ConfEnergyCalculator(SimpleConfSpace confSpace, EnergyCalculator ecalc, EnergyPartition epart, SimpleReferenceEnergies eref, boolean addResEntropy, ApproximatorMatrix amat, double approximationErrorBudget) {

		this.confSpace = confSpace;
		this.ecalc = ecalc;
		this.epart = epart;
		this.eref = eref;
		this.addResEntropy = addResEntropy;
		this.amat = amat;
		this.approximationErrorBudget = approximationErrorBudget;

		this.tasks = ecalc.tasks;
	}

	protected ConfEnergyCalculator(ConfEnergyCalculator other) {
		this(other, other.ecalc);
	}

	public ConfEnergyCalculator(ConfEnergyCalculator other, EnergyCalculator ecalc) {
		this(other.confSpace, ecalc, other.epart, other.eref, other.addResEntropy, other.amat, other.approximationErrorBudget);
	}

	/**
	 * returns the number of requested energy calculations,
	 * including ones cached in a conf DB
	 */
	public long getNumRequests() {
		return numCalculations.get() + numConfDBReads.get();
	}

	/**
	 * returns the number of energy calculations performed,
	 * excluding values cached in a conf DB
	 */
	public long getNumCalculations() {
		return numCalculations.get();
	}

	/**
	 * returns the number of energies served from a conf DB
	 */
	public long getNumConfDBReads() {
		return numConfDBReads.get();
	}

	public void resetCounters() {
		numCalculations.set(0);
		numConfDBReads.set(0);
	}
	
	public ResidueInteractions makeFragInters(RCTuple frag) {
		return EnergyPartition.makeFragment(confSpace, eref, addResEntropy, frag);
	}
	
	public ResidueInteractions makeSingleInters(int pos, int rc) {
		return epart.makeSingle(confSpace, eref, addResEntropy, pos, rc);
	}

	public ResidueInteractions makePairInters(int pos1, int rc1, int pos2, int rc2) {
		return epart.makePair(confSpace, eref, addResEntropy, pos1, rc1, pos2, rc2);
	}

	public ResidueInteractions makeTupleInters(RCTuple tuple) {
		return epart.makeTuple(confSpace, eref, addResEntropy, tuple);
	}

	public ResidueInteractions makeTripleCorrectionInters(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {
		return epart.makeTripleCorrection(confSpace, eref, addResEntropy, pos1, rc1, pos2, rc2, pos3, rc3);
	}

	public ResidueInteractions makeQuadCorrectionInters(int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4, int rc4) {
		return epart.makeQuadCorrection(confSpace, eref, addResEntropy, pos1, rc1, pos2, rc2, pos3, rc3, pos4, rc4);
	}

	public EnergyCalculator.EnergiedParametricMolecule calcSingleEnergy(int pos, int rc) {
		return calcSingleEnergy(new RCTuple(pos, rc));
	}
	
	public EnergyCalculator.EnergiedParametricMolecule calcSingleEnergy(RCTuple frag) {
		return calcEnergy(frag, epart.makeSingle(confSpace, eref, addResEntropy, frag.pos.get(0), frag.RCs.get(0)));
	}
	
	public EnergyCalculator.EnergiedParametricMolecule calcPairEnergy(int pos1, int rc1, int pos2, int rc2) {
		return calcPairEnergy(new RCTuple(pos1, rc1, pos2, rc2));
	}
	
	public EnergyCalculator.EnergiedParametricMolecule calcPairEnergy(RCTuple frag) {
		return calcEnergy(frag, epart.makePair(confSpace, eref, addResEntropy, frag.pos.get(0), frag.RCs.get(0), frag.pos.get(1), frag.RCs.get(1)));
	}

	public EnergyCalculator.EnergiedParametricMolecule calcTupleEnergy(RCTuple frag) {
		return calcEnergy(frag, epart.makeTuple(confSpace, eref, addResEntropy, frag));
	}

	/**
	 * Calculate the energy of a molecule fragment generated from a conformation space using residue interactions generated by the energy partition.
	 * 
	 * @param frag The assignments of the conformation space 
	 * @return The energy of the resulting molecule fragment and its pose
	 */
	public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag) {
		return calcEnergy(frag, makeFragInters(frag));
	}
	
	/**
	 * Calculate the energy of a molecule fragment generated from a conformation space.
	 * 
	 * @param frag The assignments of the conformation space 
	 * @param inters The residue interactions
	 * @return The energy of the resulting molecule fragment and its pose
	 */
	public EnergyCalculator.EnergiedParametricMolecule calcEnergy(RCTuple frag, ResidueInteractions inters) {

		numCalculations.incrementAndGet();
		ParametricMolecule pmol = confSpace.makeMolecule(frag);

		ResidueInteractionsApproximator approximator = null;
		if (amat != null) {
			approximator = amat.get(frag, inters, approximationErrorBudget);
		}

		return ecalc.calcEnergy(pmol, inters, approximator);
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(RCTuple,ResidueInteractions)}.
	 * 
	 * @param frag The assignments of the conformation space 
	 * @param inters The residue interactions
	 * @param listener Callback function that will receive the energy and associated molecule pose.
	 *                 Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(RCTuple frag, ResidueInteractions inters, TaskListener<EnergyCalculator.EnergiedParametricMolecule> listener) {
		tasks.submit(() -> calcEnergy(frag, inters), listener);
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(RCTuple)}.
	 *
	 * @param frag The assignments of the conformation space
	 * @param listener Callback function that will receive the energy and associated molecule pose.
	 *                 Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(RCTuple frag, TaskListener<EnergyCalculator.EnergiedParametricMolecule> listener) {
		tasks.submit(() -> calcEnergy(frag), listener);
	}

	/**
	 * Version of {@link #calcEnergy(RCTuple)}
	 * using the specified ConfDB table as a cache.
	 *
	 * @param frag The assignments of the conformation space
	 * @param table the confDB table
	 * @return The energy of the resulting molecule fragment
	 */
	public double calcEnergy(RCTuple frag, ConfDB.ConfTable table) {
		return calcEnergy(frag, makeFragInters(frag), table);
	}

	/**
	 * Version of {@link #calcEnergy(RCTuple,ResidueInteractions)}
	 * using the specified ConfDB table as a cache.
	 *
	 * @param frag The assignments of the conformation space
	 * @param inters The residue interactions
	 * @param table the confDB table
	 * @return The energy of the resulting molecule fragment
	 */
	public double calcEnergy(RCTuple frag, ResidueInteractions inters, ConfDB.ConfTable table) {

		// no confDB? just compute the energy
		if (table == null) {
			return calcEnergy(frag, inters).energy;
		}

		// check the confDB for the energy
		int[] conf = Conf.make(confSpace, frag);
		ConfDB.Conf dbconf = table.get(conf);
		if (dbconf != null && dbconf.upper != null) {
			numConfDBReads.incrementAndGet();
			return dbconf.upper.energy;
		}

		// cache miss, compute the energy
		double energy = calcEnergy(frag, inters).energy;

		// update the ConfDB
		table.setUpperBound(conf, energy, TimeTools.getTimestampNs());
		table.flush();

		return energy;
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(RCTuple,ConfDB.ConfTable)}.
	 *
	 * @param frag The assignments of the conformation space
	 * @param table the confDB table
	 * @param listener Callback function that will receive the energy.
	 *                 Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(RCTuple frag, ConfDB.ConfTable table, TaskListener<Double> listener) {
		tasks.submit(() -> calcEnergy(frag, table), listener);
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(RCTuple,ResidueInteractions,ConfDB.ConfTable)}.
	 *
	 * @param frag The assignments of the conformation space
	 * @param inters The residue interactions
	 * @param table the confDB table
	 * @param listener Callback function that will receive the energy.
	 *                 Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(RCTuple frag, ResidueInteractions inters, ConfDB.ConfTable table, TaskListener<Double> listener) {
		tasks.submit(() -> calcEnergy(frag, inters, table), listener);
	}

	/**
	 * Calculate energy of a scored conformation. Residue interactions are generated from the energy partition.
	 * 
	 * @param conf The conformation to analyze
	 * @return The conformation with attached energy
	 */
	public EnergiedConf calcEnergy(ScoredConf conf) {
		return new EnergiedConf(conf, calcEnergy(new RCTuple(conf.getAssignments())).energy);
	}

	/**
	 * Calculate energy of a scored conformation, using the specified ConfDB table as a cache.
	 * Residue interactions are generated from the energy partition.
	 *
	 * @param conf The conformation to analyze
	 * @param table the confDB table
	 * @return The conformation with attached energy
	 */
	public EnergiedConf calcEnergy(ScoredConf conf, ConfDB.ConfTable table) {
		return calcEnergy(conf, table, () -> calcEnergy(conf));
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(ScoredConf)}.
	 * 
	 * @param conf The conformation to analyze
	 * @param listener Callback function that will receive the energy. Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(ScoredConf conf, TaskListener<EnergiedConf> listener) {
		tasks.submit(() -> calcEnergy(conf), listener);
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(ScoredConf)},
	 * using the specified ConfDB table as a cache.
	 *
	 * @param conf The conformation to analyze
	 * @param table the confDB table
	 * @param listener Callback function that will receive the energy. Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(ScoredConf conf, ConfDB.ConfTable table, TaskListener<EnergiedConf> listener) {
		calcEnergyAsync(conf, table, () -> calcEnergy(conf), listener);
	}

	/**
	 * Calculate energy of a scored conformation with specified residue interactions.
	 *
	 * @param conf The conformation to analyze
	 * @param inters The residue interactions
	 * @return The conformation with attached energy
	 */
	public EnergiedConf calcEnergy(ScoredConf conf, ResidueInteractions inters) {
		return new EnergiedConf(conf, calcEnergy(new RCTuple(conf.getAssignments()), inters).energy);
	}

	/**
	 * Calculate energy of a scored conformation with specified residue interactions,
	 * using the specified ConfDB table as a cache.
	 *
	 * @param conf The conformation to analyze
	 * @param inters The residue interactions
	 * @param table the confDB table
	 * @return The conformation with attached energy
	 */
	public EnergiedConf calcEnergy(ScoredConf conf, ResidueInteractions inters, ConfDB.ConfTable table) {
		return calcEnergy(conf, table, () -> calcEnergy(conf, inters));
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(ScoredConf)}.
	 *
	 * @param conf The conformation to analyze
	 * @param inters The residue interactions
	 * @param listener Callback function that will receive the energy. Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(ScoredConf conf, ResidueInteractions inters, TaskListener<EnergiedConf> listener) {
		tasks.submit(() -> calcEnergy(conf), listener);
	}

	/**
	 * Asynchronous version of {@link #calcEnergy(ScoredConf)},
	 * using the specified ConfDB table as a cache.
	 *
	 * @param conf The conformation to analyze
	 * @param inters The residue interactions
	 * @param table the confDB table
	 * @param listener Callback function that will receive the energy. Called on a listener thread which is separate from the calling thread.
	 */
	public void calcEnergyAsync(ScoredConf conf, ResidueInteractions inters, ConfDB.ConfTable table, TaskListener<EnergiedConf> listener) {
		calcEnergyAsync(conf, table, () -> calcEnergy(conf, inters), listener);
	}


	private EnergiedConf calcEnergy(ScoredConf conf, ConfDB.ConfTable table, Supplier<EnergiedConf> supplier) {

		// no confDB? just compute the energy
		if (table == null) {
			return supplier.get();
		}

		// check the confDB for the energy
		EnergiedConf econf = table.getEnergied(conf);
		if (econf != null) {
			numConfDBReads.incrementAndGet();
			return econf;
		}

		// cache miss, compute the energy
		econf = supplier.get();

		// update the ConfDB
		// NOTE: flushing the db every write might be noticeably slow at a high write rate
		// in testing so far, at about 20 writes/s, the performance hit is undetectable
		table.setBounds(econf, TimeTools.getTimestampNs());
		table.flush();

		return econf;
	}

	private void calcEnergyAsync(ScoredConf conf, ConfDB.ConfTable table, Supplier<EnergiedConf> supplier, TaskListener<EnergiedConf> listener) {
		tasks.submit(() -> calcEnergy(conf, table, supplier), listener);
	}


	public List<EnergiedConf> calcAllEnergies(List<ScoredConf> confs) {
		return calcAllEnergies(confs, false);
	}

	public List<EnergiedConf> calcAllEnergies(List<ScoredConf> confs, ConfDB.ConfTable table) {
		return calcAllEnergies(confs, false, table);
	}
	
	public List<EnergiedConf> calcAllEnergies(List<ScoredConf> confs, boolean reportProgress) {
		return calcAllEnergies(confs, reportProgress, null);
	}

	public List<EnergiedConf> calcAllEnergies(List<ScoredConf> confs, boolean reportProgress, ConfDB.ConfTable table) {
		
		// allocate space to hold the minimized values
		List<EnergiedConf> econfs = new ArrayList<>(confs.size());
		for (int i=0; i<confs.size(); i++) {
			econfs.add(null);
		}
		
		// track progress if desired
		final Progress progress;
		if (reportProgress) {
			progress = new Progress(confs.size());
		} else {
			progress = null;
		}
		
		// minimize them all
		for (int i=0; i<confs.size(); i++) {
			
			// capture i for the closure below
			final int fi = i;
			
			calcEnergyAsync(confs.get(i), table, (econf) -> {
				
				// save the minimized energy
				econfs.set(fi, econf);
				
				// update progress if needed
				if (progress != null) {
					progress.incrementProgress();
				}
			});
		}
		tasks.waitForFinish();
		
		return econfs;
	}

	//Making objective functions for EPIC fitting
	public MoleculeObjectiveFunction makeIntraShellObjFcn(int pos, int rc) {
		ParametricMolecule bpmol = confSpace.makeMolecule(new RCTuple(pos,rc));
		ResidueInteractions inters = EnergyPartition.Traditional.makeSingle(confSpace, eref, addResEntropy, pos, rc);
		return ecalc.makeEnergyObjFcn(bpmol, inters);
	}

	public MoleculeObjectiveFunction makePairwiseObjFcn(int pos1, int rc1, int pos2, int rc2) {
		ParametricMolecule bpmol = confSpace.makeMolecule(new RCTuple(pos1,rc1,pos2,rc2));
		ResidueInteractions inters = EnergyPartition.Traditional.makePair(confSpace, eref, addResEntropy, pos1, rc1, pos2, rc2);
		return ecalc.makeEnergyObjFcn(bpmol, inters);
	}
}
