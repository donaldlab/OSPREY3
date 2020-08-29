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

package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.KStar.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.*;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarBound;
import edu.duke.cs.osprey.sharkstar.SHARKSeqHScorer;
import edu.duke.cs.osprey.sharkstar.SingleSequenceSHARKStarBound;
import edu.duke.cs.osprey.tools.BigMath;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;
import java.util.function.BiFunction;
import java.util.function.Function;


/**
 * Implementation of the BBK* algorithm to predict protein sequence mutations that improve
 * binding affinity by computing provably accurate Boltzmann-weighted ensembles
 * {@cite Ojewole2008 Adegoke A. Ojewole, Jonathan D. Jou, Vance G. Fowler, Jr., and Bruce R. Donald, 2017.
 * BBK* (Branch and Bound over K*): A Provable and Efficient Ensemble-Based Algorithm to Optimize
 * Stability and Binding Affinity over Large Sequence Spaces
 * In Research in Computational Molecular Biology (accepted, in press).}.
 */
public class BBKStar {

	// TODO: update BBK* implementation to use new SeqAStarTree?

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		public static class Builder {

			/** The number of best (by K* score) sequences to evaluate before finishing */
			private int numBestSequences = 1;

			/**
			 * The number of conformations to evaluate per batch of partition function refinement
			 *
			 * For best results, make this a multiple of the available parallelism. e.g., if using 4 threads,
			 * try a batch size of 4, 8, 12, 16, etc.
			 */
			private int numConfsPerBatch = 8;

			/**
			 * The maximum number of conformations to evaluate per batch of partition function refinement
			 *
			 * For best results, make this a multiple of the available parallelism. e.g., if using 4 threads,
			 * try a batch size of 4, 8, 12, 16, etc.
			 */
			private int maxNumConfsPerBatch = 8;

			public Builder setNumBestSequences(int val) {
				numBestSequences = val;
				return this;
			}

			public Builder setNumConfsPerBatch(int val) {
				numConfsPerBatch = val;
				return this;
			}

			public Builder setMaxNumConfsPerBatch(int val) {
				maxNumConfsPerBatch = val;
				return this;
			}

			public Settings build() {
				return new Settings(numBestSequences, numConfsPerBatch, maxNumConfsPerBatch);
			}
		}

		public final int numBestSequences;
		public final int numConfsPerBatch;
		private final int maxNumConfsPerBatch;

		public Settings(int numBestSequences, int numConfsPerBatch) {
			this.numBestSequences = numBestSequences;
			this.numConfsPerBatch = numConfsPerBatch;
			this.maxNumConfsPerBatch = numConfsPerBatch;
		}
		public Settings(int numBestSequences, int numConfsPerBatch, int maxNumConfsPerBatch) {
			this.numBestSequences = numBestSequences;
			this.numConfsPerBatch = numConfsPerBatch;
			this.maxNumConfsPerBatch = maxNumConfsPerBatch;
		}
	}

	public class ConfSpaceInfo {

		public final SimpleConfSpace confSpace;
		public final KStar.ConfSpaceType type;
		public final String id;

		/** A ConfEnergyCalculator that computes minimized energies */
		public ConfEnergyCalculator confEcalcMinimized = null;

		/** A ConfSearch that minimizes over conformation scores from minimized tuples */
		public ConfSearchFactory confSearchFactoryMinimized = null;

		/** A ConfSearch that maximizes over conformation scores from rigid tuples */
		public ConfSearchFactory confSearchFactoryRigid = null;

		/** A class to manage and create the correct partition functions */
		public PartitionFunctionFactory pfuncFactory = null;

		public File confDBFile = null;
		public ConfEnergyCalculator confEcalcRigid;

		private BigDecimal stabilityThreshold = null;

		public ConfSpaceInfo(SimpleConfSpace confSpace, KStar.ConfSpaceType type) {
			this.confSpace = confSpace;
			this.type = type;
			this.id = type.name().toLowerCase();
		}

		private void check() {
			if (confEcalcMinimized == null) {
				throw new KStar.InitException(type, "confEcalcMinimized");
			}
			if (pfuncFactory == null) {
				throw new KStar.InitException(type, "pfuncFactory");
			}
		}

		public void setConfDBFile(String path) {
			confDBFile = new File(path);
		}
	}

	private class ConfDBs {
		public ConfDB protein = null;
		public ConfDB ligand = null;
		public ConfDB complex = null;
	}

	private static interface DBsUser {
		void use(ConfDBs confdbs);
	}

	public static enum PfuncsStatus {
		Estimating,
		Estimated,
		Blocked
	}

	private abstract class Node implements Comparable<Node> {

		public final Sequence sequence;
		public final ConfDB.DBs confDBs;

		/** for comparing in the tree, higher is first */
		public double score;

		/** signals whether or not partition function values are allowed the stability threshold */
		public boolean isUnboundUnstable;

		/** counts how many estimations the node has completed*/
		public long numComputeCycles = 0;

		protected Node(Sequence sequence, ConfDB.DBs confDBs) {
			this.sequence = sequence;
			this.confDBs = confDBs;
			this.score = 0;
		}

		@Override
		public int compareTo(Node other) {
			// negate for descending sort
			return -Double.compare(this.score, other.score);
		}

		public abstract void estimateScore();
	}

	public class MultiSequenceNode extends Node {

		public MultiSequenceNode(Sequence sequence, ConfDB.DBs confDBs) {
			super(sequence, confDBs);
		}

		public List<Node> makeChildren() {

			List<Node> children = new ArrayList<>();

			// pick the next design position
			// TODO: dynamic A*?
			List<SeqSpace.Position> positions = complex.confSpace.seqSpace.positions;
			SeqSpace.Position assignPos = positions.stream()
				.filter((pos) -> !sequence.isAssigned(pos.resNum))
				.findFirst()
				.orElseThrow(() -> new IllegalStateException("no design positions left to choose"));

			// get the possible assignments
			Set<SeqSpace.ResType> resTypes = new HashSet<>(assignPos.resTypes);

			// add wild-type option if mutations are limited
			if (kstarSettings.maxSimultaneousMutations < positions.size()) {
				resTypes.add(assignPos.wildType);
			}

			// for each assignment...
			for (SeqSpace.ResType resType : resTypes) {

				// update the sequence with this assignment
				Sequence s = sequence.copy().set(assignPos, resType);

				if (s.isFullyAssigned()) {

					// fully assigned, make single sequence node
					children.add(new SingleSequenceNode(s, confDBs));

				} else if (s.countMutations() == kstarSettings.maxSimultaneousMutations) {

					// mutation limit reached, fill unassigned positions with wild-type
					s.fillWildType();
					if (s.isFullyAssigned()) {
						children.add(new SingleSequenceNode(s, confDBs));
					}

					// NOTE: if we didn't fill the assignments, it means there aren't enough wild-types to do it
					// so don't explore that sequence

				} else {

					// still partial sequence, make multi-sequence node
					children.add(new MultiSequenceNode(s, confDBs));
				}
			}

			return children;
		}

		@Override
		public void estimateScore() {

			// TODO: expose setting?
			// NOTE: for the correctness of the bounds, the number of confs must be the same for every node
			// meaning, it might not be sound to do epsilon-based iterative approximations here
			final int numConfs = 1000;

			if (protein.stabilityThreshold != null) {

				// tank the sequence if the protein is unstable
				BigDecimal proteinUpperBound = calcUpperBound(protein, sequence, numConfs);
				if (MathTools.isLessThan(proteinUpperBound, protein.stabilityThreshold)) {
					score = Double.NEGATIVE_INFINITY;
					isUnboundUnstable = true;
					return;
				}
			}

			BigDecimal proteinLowerBound = calcLowerBoundByConf(protein, sequence, numConfs);

			// if the first few conf upper bound scores (for the pfunc lower bound) are too high,
			// then the K* upper bound is also too high
			if (MathTools.isZero(proteinLowerBound)) {
				score = Double.POSITIVE_INFINITY;
				isUnboundUnstable = false;
				return;
			}

			if (ligand.stabilityThreshold != null) {

				// tank the sequence if the ligand is unstable
				BigDecimal ligandUpperBound = calcUpperBound(ligand, sequence, numConfs);
				if (MathTools.isLessThan(ligandUpperBound, ligand.stabilityThreshold)) {
					score = Double.NEGATIVE_INFINITY;
					isUnboundUnstable = true;
					return;
				}
			}

			BigDecimal ligandLowerBound = calcLowerBoundByConf(ligand, sequence, numConfs);

			// if the first few conf upper bound scores (for the pfunc lower bound) are too high,
			// then the K* upper bound is also too high
			if (MathTools.isZero(ligandLowerBound)) {
				score = Double.POSITIVE_INFINITY;
				isUnboundUnstable = false;
				return;
			}

			BigDecimal complexUpperBound = calcUpperBound(complex, sequence, numConfs);

			// compute the node score
			score = MathTools.bigDivideDivide(
				complexUpperBound,
				proteinLowerBound,
				ligandLowerBound,
				PartitionFunction.decimalPrecision
			).doubleValue();
			isUnboundUnstable = false;
		}

		private BigDecimal calcUpperBoundBySumProduct(ConfSpaceInfo info, Sequence sequence) {
			RCs rcs = sequence.makeRCs(info.confSpace);
			EnergyMatrix rigidEmat = info.pfuncFactory.getOrMakeEmat(info.confEcalcRigid,
					info.id+"rigid");
			EnergyMatrix minimizedEmat = info.pfuncFactory.getOrMakeEmat(info.confEcalcMinimized,
					info.id+".minimized");
			SHARKSeqHScorer scorer = new SHARKSeqHScorer(info.confSpace.seqSpace, info.confSpace,
					rigidEmat, minimizedEmat);
			return scorer.calcBigDecimal(info.confSpace, sequence);
		}

		private BigDecimal calcLowerBoundByConf(ConfSpaceInfo info, Sequence sequence, int numConfs) {

			// to compute lower bounds on pfuncs, we'll do the usual lower bound calculation,
			// but use rigid energies instead of minimized energies

			RCs rcs = sequence.makeRCs(info.confSpace);
			ConfSearch astar = info.confSearchFactoryRigid.make(rcs);
			BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

			BigMath m = new BigMath(PartitionFunction.decimalPrecision)
				.set(0.0);
			for (int i=0; i<numConfs; i++) {
				ConfSearch.ScoredConf conf = astar.nextConf();
				if (conf == null) {
					break;
				}

				m.add(bcalc.calc(conf.getScore()));
			}
			return m.get();
		}

		private BigDecimal calcUpperBound(ConfSpaceInfo info, Sequence sequence, int numConfs) {
		    switch(info.pfuncFactory.getPfuncImpl()) {
				case SHARKStar:
					return calcUpperBoundBySumProduct(info, sequence);
				default:
					return calcUpperBoundByConf(info, sequence, numConfs);

			}
		}

		private BigDecimal calcUpperBoundByConf (ConfSpaceInfo info, Sequence sequence, int numConfs) {

			// to compute upper bounds on pfuncs,
			// we'll use the upper bound calculator in the usual way

			RCs rcs = sequence.makeRCs(info.confSpace);
			UpperBoundCalculator calc = new UpperBoundCalculator(
				info.confSearchFactoryMinimized.make(rcs),
				rcs.getNumConformations()
			);
			calc.run(numConfs);
			return calc.totalBound;
		}

		@Override
		public String toString() {
			return String.format("MultiSequenceNode[score=%12.6f, seq=%s]",
				score,
				sequence
			);
		}
	}

	public class SingleSequenceNode extends Node {

		public final PartitionFunction protein;
		public final PartitionFunction ligand;
		public final PartitionFunction complex;


		public SingleSequenceNode(Sequence sequence, ConfDB.DBs confDBs) {
			super(sequence, confDBs);

			// make the partition functions
			this.protein = makePfunc(proteinPfuncs, BBKStar.this.protein, confDBs.get(BBKStar.this.protein.confSpace));
			this.ligand = makePfunc(ligandPfuncs, BBKStar.this.ligand, confDBs.get(BBKStar.this.ligand.confSpace));
			this.complex = makePfunc(complexPfuncs, BBKStar.this.complex, confDBs.get(BBKStar.this.complex.confSpace));
		}

		private PartitionFunction makePfunc(Map<Sequence,PartitionFunction> pfuncCache, ConfSpaceInfo info, ConfDB confdb) {

			// filter the global sequence to this conf space
			Sequence sequence = this.sequence.filter(info.confSpace.seqSpace);

			// first check the cache
			PartitionFunction pfunc = pfuncCache.get(sequence);
			if (pfunc != null) {
				return pfunc;
			}

			// cache miss, need to compute the partition function

			// make the partition function
			RCs rcs = sequence.makeRCs(info.confSpace);

			pfunc = info.pfuncFactory.makePartitionFunctionFor(rcs, rcs.getNumConformations(), kstarSettings.epsilon, sequence);

			pfunc.setReportProgress(kstarSettings.showPfuncProgress);
			if (confdb != null) {
				PartitionFunction.WithConfTable.setOrThrow(pfunc, confdb.getSequence(sequence));
			}
			if (kstarSettings.useExternalMemory) {
				PartitionFunction.WithExternalMemory.setOrThrow(pfunc, true, rcs);
			}
			pfunc.setStabilityThreshold(info.stabilityThreshold);

			// update the cache
			pfuncCache.put(sequence, pfunc);
			return pfunc;
		}

		@Override
		public void estimateScore() {

			// tank the sequence if either unbound strand is unstable
			// yeah, we haven't refined any pfuncs yet this estimation,
			// but since pfuncs get cached, check before we do any more estimation
			if (protein.getStatus() == PartitionFunction.Status.Unstable
				|| ligand.getStatus() == PartitionFunction.Status.Unstable) {
				score = Double.NEGATIVE_INFINITY;
				isUnboundUnstable = true;
				return;
			}

			/* Determine how long we want to compute for.
			There is a trade-off -- computing for longer on bad sequences will take more time,
			but computing for less time wastes resources
			 */

			Function<Long, Integer> weightFunc = (a) ->{
				if (a >= 5)
					return bbkstarSettings.maxNumConfsPerBatch;
				else
					return bbkstarSettings.numConfsPerBatch;
			};

			int numConfsToCompute = weightFunc.apply(numComputeCycles);

			// refine the pfuncs if needed
			if (protein.getStatus().canContinue()) {
				//protein.compute(bbkstarSettings.numConfsPerBatch);
				protein.compute(numConfsToCompute);

				// tank the sequence if the unbound protein is unstable
				if (protein.getStatus() == PartitionFunction.Status.Unstable) {
					score = Double.NEGATIVE_INFINITY;
					isUnboundUnstable = true;
					return;
				}
			}

			if (ligand.getStatus().canContinue()) {
				//ligand.compute(bbkstarSettings.numConfsPerBatch);
				ligand.compute(numConfsToCompute);

				// tank the sequence if the unbound ligand is unstable
				if (ligand.getStatus() == PartitionFunction.Status.Unstable) {
					score = Double.NEGATIVE_INFINITY;
					isUnboundUnstable = true;
					return;
				}
			}

			if (complex.getStatus().canContinue()) {
				//complex.compute(bbkstarSettings.numConfsPerBatch);
				complex.compute(numConfsToCompute);
			}

			// update the score
			score = Math.log10(makeKStarScore().upperBound.doubleValue());
			isUnboundUnstable = false;

			// tank sequences that have no useful K* bounds, and are blocked
			if (getStatus() == PfuncsStatus.Blocked && score == Double.POSITIVE_INFINITY) {
				score = Double.NEGATIVE_INFINITY;
			}

			// Increment the numComputeCycles counter
			this.numComputeCycles+= numConfsToCompute / bbkstarSettings.numConfsPerBatch;
		}

		public KStarScore computeScore() {

			// refine the pfuncs until done
			while (protein.getStatus().canContinue()) {
				//protein.compute(bbkstarSettings.numConfsPerBatch);
				protein.compute();
			}
			while (ligand.getStatus().canContinue()) {
				//ligand.compute(bbkstarSettings.numConfsPerBatch);
				ligand.compute();
			}
			while (complex.getStatus().canContinue()) {
				complex.compute();
			}

			// update the score
			KStarScore kstarScore = makeKStarScore();
			score = Math.log10(kstarScore.upperBound.doubleValue());
			return kstarScore;
		}

		public KStarScore makeKStarScore() {
			return new KStarScore(protein.makeResult(), ligand.makeResult(), complex.makeResult());
		}

		public PfuncsStatus getStatus() {

			// aggregate pfunc statuses
			if (protein.getStatus() == PartitionFunction.Status.Estimated
				&& ligand.getStatus() == PartitionFunction.Status.Estimated
				&& complex.getStatus() == PartitionFunction.Status.Estimated) {
				return PfuncsStatus.Estimated;
			} else if (protein.getStatus() == PartitionFunction.Status.Estimating
				|| ligand.getStatus() == PartitionFunction.Status.Estimating
				|| complex.getStatus() == PartitionFunction.Status.Estimating) {
				return PfuncsStatus.Estimating;
			} else {
				return PfuncsStatus.Blocked;
			}
		}

		@Override
		public String toString() {
			return String.format("SingleSequenceNode[score=%12.6f, seq=%s, K*=%s]",
				score,
				sequence,
				makeKStarScore()
			);
		}
	}

	/** A configuration space containing just the protein strand */
	public final ConfSpaceInfo protein;

	/** A configuration space containing just the ligand strand */
	public final ConfSpaceInfo ligand;

	/** A configuration space containing both the protein and ligand strands */
	public final ConfSpaceInfo complex;

	/** Optional and overridable settings for BBK*, shared with K* */
	public final KStar.Settings kstarSettings;

	/** Optional and overridable settings for BBK* */
	public final Settings bbkstarSettings;

	/** Partition Function manager */
	PartitionFunctionFactory pfuncFactory;

	// TODO: caching these will keep lots of A* trees in memory. is that a problem? (oh yes, it definitely is)
	private final Map<Sequence,PartitionFunction> proteinPfuncs;
	private final Map<Sequence,PartitionFunction> ligandPfuncs;
	private final Map<Sequence,PartitionFunction> complexPfuncs;

	public MultiSequenceSHARKStarBound complexSHARK; // temporary variable for SHARK* test information
	public MultiSequenceSHARKStarBound proteinSHARK; // temporary variable for SHARK* test information
	public MultiSequenceSHARKStarBound ligandSHARK; // temporary variable for SHARK* test information

	public BBKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, KStar.Settings kstarSettings, Settings bbkstarSettings) {

		// BBK* doesn't work with external memory (never enough internal memory for all the priority queues)
		if (kstarSettings.useExternalMemory) {
			throw new IllegalArgumentException("BBK* is not compatible with external memory."
				+ " Please switch to regular K* with external memory, or keep using BBK* and disable external memory.");
		}

		this.protein = new ConfSpaceInfo(protein, KStar.ConfSpaceType.Protein);
		this.ligand = new ConfSpaceInfo(ligand, KStar.ConfSpaceType.Ligand);
		this.complex = new ConfSpaceInfo(complex, KStar.ConfSpaceType.Complex);
		this.kstarSettings = kstarSettings;
		this.bbkstarSettings = bbkstarSettings;

		proteinPfuncs = new HashMap<>();
		ligandPfuncs = new HashMap<>();
		complexPfuncs = new HashMap<>();
	}

	public Iterable<ConfSpaceInfo> confSpaceInfos() {
		return Arrays.asList(protein, ligand, complex);
	}

	public BBKStar.ConfSpaceInfo getConfSpaceInfo(SimpleConfSpace confSpace) {
		if (confSpace == protein.confSpace) {
			return protein;
		} else if (confSpace == ligand.confSpace) {
			return ligand;
		} else if (confSpace == complex.confSpace) {
			return complex;
		} else {
			throw new IllegalArgumentException("conf space does not match any known by this K* instance");
		}
	}

	public List<KStar.ScoredSequence> run() {

		protein.check();
		ligand.check();
		complex.check();

		// clear any previous state
		proteinPfuncs.clear();
		ligandPfuncs.clear();
		complexPfuncs.clear();

		List<KStar.ScoredSequence> scoredSequences = new ArrayList<>();
		PriorityQueue<Node> tree;
		List<SingleSequenceNode> finishedNodes = new ArrayList<>();
		List<SingleSequenceNode> unstableNodes = new ArrayList<>();

		// open the conf databases if needed
		try (ConfDB.DBs confDBs = new ConfDB.DBs()
			 .add(protein.confSpace, protein.confDBFile)
			.add(ligand.confSpace, ligand.confDBFile)
			.add(complex.confSpace, complex.confDBFile)
		) {

			// calculate wild-type first
			if (complex.confSpace.seqSpace.containsWildTypeSequence()) {
				System.out.println("computing K* score for the wild-type sequence...");
				SingleSequenceNode wildTypeNode = new SingleSequenceNode(complex.confSpace.makeWildTypeSequence(), confDBs);
				KStarScore wildTypeScore = wildTypeNode.computeScore();
				kstarSettings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
					-1,
					0,
					wildTypeNode.sequence,
					wildTypeScore
				));
				if (kstarSettings.stabilityThreshold != null) {
					BigDecimal stabilityThresholdFactor = new BoltzmannCalculator(PartitionFunction.decimalPrecision).calc(kstarSettings.stabilityThreshold);
					protein.stabilityThreshold = wildTypeScore.protein.values.calcLowerBound().multiply(stabilityThresholdFactor);
					ligand.stabilityThreshold = wildTypeScore.ligand.values.calcLowerBound().multiply(stabilityThresholdFactor);
				}
				// record the SHARK*bound for test information
				try{
					SingleSequenceSHARKStarBound complexBound = (SingleSequenceSHARKStarBound) wildTypeNode.complex;
					SingleSequenceSHARKStarBound proteinBound = (SingleSequenceSHARKStarBound) wildTypeNode.protein;
					SingleSequenceSHARKStarBound ligandBound = (SingleSequenceSHARKStarBound) wildTypeNode.ligand;
					this.complexSHARK = complexBound.multiSequenceSHARKStarBound;
					this.proteinSHARK = proteinBound.multiSequenceSHARKStarBound;
					this.ligandSHARK = ligandBound.multiSequenceSHARKStarBound;
				}catch(Exception ignored){}

			} else if (kstarSettings.stabilityThreshold != null) {
				System.out.println("Sequence space does not contain the wild type sequence, stability threshold is disabled");
			}


			// start the BBK* tree with the root node
			tree = new PriorityQueue<>();
			tree.add(new MultiSequenceNode(complex.confSpace.makeUnassignedSequence(), confDBs));

			// Ensemble stat hack
			PartitionFunction lastPfunc = null;
			// Hack for Time
			SingleSequenceNode lastNode = null;
			// To reduce output spam
			Sequence lastSequence = null;

			// start searching the tree
			System.out.println("computing K* scores for the " + bbkstarSettings.numBestSequences + " best sequences to epsilon = " + kstarSettings.epsilon + " ...");
			kstarSettings.scoreWriters.writeHeader();
			while (!tree.isEmpty()) {
				// termination criterion
				if(scoredSequences.size() >= bbkstarSettings.numBestSequences){
					// Since it is no longer guaranteed that partition function calculation is synchronous, we may return
					// sequences in the wrong order if we use a more asynchronous pfunc calculator. So, we check.
					scoredSequences.sort((a, b) -> Double.compare(b.score.upperBoundLog10(), a.score.upperBoundLog10()));
					scoredSequences.forEach(System.out::println);
					double worstFinishedUB = scoredSequences.get(bbkstarSettings.numBestSequences-1).score.upperBoundLog10();
					Node bestUnfinishedNode = tree.peek();
					double bestUnfinishedUB = bestUnfinishedNode.score; // this is actually the upper bound
					System.out.println(String.format("Worst finished UB: %s --> %.6f, Best unfinished UB: %s --> %.6f",
							scoredSequences.get(bbkstarSettings.numBestSequences-1).sequence,
							worstFinishedUB,
							bestUnfinishedNode.sequence,
							bestUnfinishedUB
							));
					if(worstFinishedUB >= bestUnfinishedUB)
						break;
				}

				// get the next node
				Node node = tree.poll();
				if(lastSequence != node.sequence)
					System.out.println("Refining sequence "+node.sequence);
				lastSequence = node.sequence;

				if (node instanceof SingleSequenceNode) {
					SingleSequenceNode ssnode = (SingleSequenceNode)node;

					// single-sequence node
					switch (ssnode.getStatus()) {
						case Estimated:

							// sequence is finished, return it!
							reportSequence(ssnode, scoredSequences);
							lastNode = ssnode;
							lastPfunc = ssnode.complex;
							finishedNodes.add(ssnode);

						break;
						case Estimating:

							// needs more estimation, catch-and-release
							ssnode.estimateScore();
							if (!ssnode.isUnboundUnstable) {
								tree.add(ssnode);
							}else{
								unstableNodes.add(ssnode);
							}

						break;
						case Blocked:

							// from here on out, it's all blocked sequences
							// so it's ok to put them in the sorted order now
							reportSequence(ssnode, scoredSequences);
					}

				} else if (node instanceof MultiSequenceNode) {
					MultiSequenceNode msnode = (MultiSequenceNode)node;

					// partial sequence, expand children
					// TODO: parallelize the multi-sequence node scoring here?
					for (Node child : msnode.makeChildren()) {
						child.estimateScore();
						if (!child.isUnboundUnstable) {
							tree.add(child);
						}else if(child instanceof SingleSequenceNode){
							unstableNodes.add((SingleSequenceNode) child);
						}else{
							System.out.println("Throwing away unstable multi-sequence node " + child.sequence);
						}
					}
				}
			}

			if (scoredSequences.size() < bbkstarSettings.numBestSequences) {
				if (tree.isEmpty()) {
					// all is well, we just don't have that many sequences in the design
					System.out.println("Tried to find " + bbkstarSettings.numBestSequences + " sequences,"
						+ " but design flexibility and sequence filters only allowed " + scoredSequences.size() + " sequences.");
				} else {
					throw new Error("BBK* ended, but the tree isn't empty and we didn't return enough sequences. This is a bug.");
				}
			}

			if(lastNode != null) {
				System.out.println("Trying to print ensemble stats for protein.");
				lastNode.protein.printStats();
				System.out.println("Trying to print ensemble stats for ligand.");
				lastNode.ligand.printStats();
				System.out.println("Trying to print ensemble stats for complex.");
				lastNode.complex.printStats();
				//lastPfunc.printStats();
			}
		}
		countCycles(tree, finishedNodes);
		System.out.println("Finished node information:");
		finishedNodes.forEach((SingleSequenceNode n) -> {
			System.out.println(String.format("%s protein:", n.sequence));
		    n.protein.printStats();
			System.out.println(String.format("%s ligand:", n.sequence));
			n.ligand.printStats();
			System.out.println(String.format("%s complex:", n.sequence));
			n.complex.printStats();
		});


		return scoredSequences;
	}

	private void reportSequence(SingleSequenceNode ssnode, List<KStar.ScoredSequence> scoredSequences) {

		KStarScore kstarScore = ssnode.makeKStarScore();
		scoredSequences.add(new KStar.ScoredSequence(ssnode.sequence, kstarScore));

		kstarSettings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
			scoredSequences.size() - 1,
			bbkstarSettings.numBestSequences,
			ssnode.sequence,
			kstarScore
		));
	}

	private void countCycles(PriorityQueue<Node> tree, List<SingleSequenceNode> finishedNodes){
		List<Node> nodes = new ArrayList<>(tree);
		nodes.sort(new Comparator<Node>() {
			@Override
			public int compare(Node o1, Node o2) {
				return Long.compare(o1.numComputeCycles, o2.numComputeCycles);
			}
		});
		nodes.forEach((n) -> System.out.println(String.format("%s: %d", n.sequence, n.numComputeCycles)));
		finishedNodes.forEach((n) -> System.out.println(String.format("%s: %d", n.sequence, n.numComputeCycles)));

	}
}
