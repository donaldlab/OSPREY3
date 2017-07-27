package edu.duke.cs.osprey.kstar;

import com.sun.xml.internal.bind.v2.runtime.unmarshaller.XsiNilLoader;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;

import java.io.File;
import java.math.BigDecimal;
import java.util.*;

public class BBKStar {

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		public static class Builder {

			private double epsilon = 0.683;
			private int maxSimultaneousMutations = 1;
			private boolean showPfuncProgress = false;
			private int numBestSequences = 1;
			private int numConfsPerBatch = 8;

			public Builder setEpsilon(double val) {
				epsilon = val;
				return this;
			}

			public Builder setMaxSimultaneousMutations(int val) {
				maxSimultaneousMutations = val;
				return this;
			}

			public Builder setShowPfuncProgress(boolean val) {
				showPfuncProgress = val;
				return this;
			}

			public Builder setNumBestSequences(int val) {
				numBestSequences = val;
				return this;
			}

			public Builder setNumConfsPerBatch(int val) {
				numConfsPerBatch = val;
				return this;
			}

			public Settings build() {
				return new Settings(epsilon, maxSimultaneousMutations, showPfuncProgress, numBestSequences, numConfsPerBatch);
			}
		}

		public final double epsilon;
		public final int maxSimultaneousMutations;
		public final boolean showPfuncProgress;
		public final int numBestSequences;
		public final int numConfsPerBatch;

		public Settings(double epsilon, int maxSimultaneousMutations, boolean dumpPfuncConfs, int numBestSequences, int numConfsPerBatch) {
			this.epsilon = epsilon;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.showPfuncProgress = dumpPfuncConfs;
			this.numBestSequences = numBestSequences;
			this.numConfsPerBatch = numConfsPerBatch;
		}
	}

	public static enum ConfSpaceType {
		Protein,
		Ligand,
		Complex
	}

	public class ConfSpaceInfo {

		public final ConfSpaceType type;
		public final SimpleConfSpace confSpace;
		public final ConfEnergyCalculator confEcalc;

		public EnergyMatrix emat = null;

		public ConfSpaceInfo(ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator confEcalc) {
			this.type = type;
			this.confSpace = confSpace;
			this.confEcalc = confEcalc;
		}

		public void calcEmat() {
			emat = new SimplerEnergyMatrixCalculator.Builder(confEcalc)
				.setCacheFile(new File(String.format("/tmp/emat.%s.dat", type))) // TEMP
				.build()
				.calcEnergyMatrix();
		}

		public void fillWildType(KStar.Sequence sequence) {
			for (SimpleConfSpace.Position pos : confSpace.positions) {
				if (sequence.get(pos.index) == null) {
					sequence.set(pos.index, pos.strand.mol.getResByPDBResNumber(pos.resNum).template.name);
				}
			}
		}

		public KStar.Sequence makeWildTypeSequence() {
			KStar.Sequence sequence = new KStar.Sequence(confSpace.positions.size());
			fillWildType(sequence);
			return sequence;
		}
	}

	public static enum PfuncsStatus {
		Estimating,
		Estimated,
		Blocked
	}

	private class Pfuncs {

		public PartitionFunction protein;
		public PartitionFunction ligand;
		public PartitionFunction complex;

		public void estimateMore(int numConfs) {
			if (protein.getStatus().canContinue()) {
				protein.compute(numConfs);
			}
			if (ligand.getStatus().canContinue()) {
				ligand.compute(numConfs);
			}
			if (complex.getStatus().canContinue()) {
				complex.compute(numConfs);
			}
		}

		public void compute() {
			if (protein.getStatus().canContinue()) {
				protein.compute();
			}
			if (ligand.getStatus().canContinue()) {
				ligand.compute();
			}
			if (complex.getStatus().canContinue()) {
				complex.compute();
			}
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

		public BigDecimal calcKStarUpperBound() {
			return PartitionFunction.calcKStarUpperBound(protein.makeResult(), ligand.makeResult(), complex.makeResult());
		}

		public BigDecimal calcKStarLowerBound() {
			return PartitionFunction.calcKStarLowerBound(protein.makeResult(), ligand.makeResult(), complex.makeResult());
		}

		public BigDecimal calcKStarScore() {
			return PartitionFunction.calcKStarScore(protein.makeResult(), ligand.makeResult(), complex.makeResult());
		}
	}

	private abstract class Node implements Comparable<Node> {

		/** for comparing in the tree, higher is first */
		public double score = 0;

		@Override
		public int compareTo(Node other) {
			// negate for descending sort
			return -Double.compare(this.score, other.score);
		}

		public abstract void estimateScore();
	}

	private abstract class InteriorNode extends Node {

		public final KStar.Sequence partialSequence;

		public InteriorNode(KStar.Sequence partialSequence) {
			this.partialSequence = partialSequence;
		}

		public List<Node> makeChildren() {

			List<Node> children = new ArrayList<>();

			List<SimpleConfSpace.Position> assignedPositions = getAssignedPositions();
			List<SimpleConfSpace.Position> unassignedPositions = getUnassignedPositions(assignedPositions);
			int numMutations = assignedPositions.size();

			// are we assigning the last mutation?
			if (numMutations == settings.maxSimultaneousMutations - 1) {

				// yup, make single-sequence nodes
				for (SimpleConfSpace.Position pos : unassignedPositions) {
					for (String resType : pos.resFlex.resTypes) {

						// we're choosing a mutation, so skip the wild type
						if (resType.equals(pos.resFlex.wildType)) {
							continue;
						}

						// assign this position
						KStar.Sequence sequence = new KStar.Sequence(partialSequence);
						sequence.set(pos.index, resType);

						// fill the rest of the positions with wild type
						complex.fillWildType(sequence);
						children.add(new SingleSequenceNode(sequence));
					}
				}

			} else {

				// nope, make multi-sequence nodes
				for (SimpleConfSpace.Position pos : unassignedPositions) {
					for (String resType : pos.resFlex.resTypes) {
						children.add(new MultiSequenceNode(this, pos, resType));
					}
				}
			}

			return children;
		}

		public abstract List<SimpleConfSpace.Position> getAssignedPositions();

		public List<SimpleConfSpace.Position> getUnassignedPositions(List<SimpleConfSpace.Position> assignedPositions) {
			List<SimpleConfSpace.Position> positions = new ArrayList<>();
			for (SimpleConfSpace.Position pos : complex.confSpace.positions) {
				if (!assignedPositions.contains(pos)) {
					positions.add(pos);
				}
			}
			return positions;
		}
	}

	private class RootNode extends InteriorNode {

		public RootNode() {
			super(new KStar.Sequence(complex.confSpace.positions.size()));
		}

		@Override
		public List<SimpleConfSpace.Position> getAssignedPositions() {
			return new ArrayList<>();
		}

		@Override
		public void estimateScore() {
			// nothing to do
		}

		@Override
		public String toString() {
			return "RootNode";
		}
	}

	private class MultiSequenceNode extends InteriorNode {

		public final InteriorNode parent;
		public final SimpleConfSpace.Position pos;
		public final String resType;

		public MultiSequenceNode(InteriorNode parent, SimpleConfSpace.Position pos, String resType) {
			super(new KStar.Sequence(parent.partialSequence));

			this.parent = parent;
			this.pos = pos;
			this.resType = resType;

			// assign this position
			this.partialSequence.set(pos.index, resType);
		}

		@Override
		public List<SimpleConfSpace.Position> getAssignedPositions() {
			List<SimpleConfSpace.Position> positions = parent.getAssignedPositions();
			positions.add(pos);
			return positions;
		}

		@Override
		public void estimateScore() {
			// TODO
		}

		@Override
		public String toString() {
			return String.format("MultiSequenceNode[score=%12.6f, seq=%s]", score, partialSequence.toString());
		}
	}

	private class SingleSequenceNode extends Node {

		public final KStar.Sequence sequence;
		public final Pfuncs kstarLower;

		public SingleSequenceNode(KStar.Sequence sequence) {

			this.sequence = sequence;

			// split complex sequence into protein/ligand sequences
			KStar.Sequence proteinSequence = protein.makeWildTypeSequence();
			KStar.Sequence ligandSequence = ligand.makeWildTypeSequence();
			for (SimpleConfSpace.Position pos : complex.confSpace.positions) {

				SimpleConfSpace.Position proteinPos = complexToProteinMap.get(pos);
				if (proteinPos != null) {
					proteinSequence.set(proteinPos.index, sequence.get(pos.index));
				}

				SimpleConfSpace.Position ligandPos = complexToLigandMap.get(pos);
				if (ligandPos != null) {
					ligandSequence.set(ligandPos.index, sequence.get(pos.index));
				}
			}

			this.kstarLower = new Pfuncs();
			this.kstarLower.protein = makePfunc(proteinPfuncs, protein, proteinSequence);
			this.kstarLower.ligand = makePfunc(ligandPfuncs, ligand, ligandSequence);
			this.kstarLower.complex = makePfunc(complexPfuncs, complex, sequence);
		}

		public PartitionFunction makePfunc(Map<KStar.Sequence,PartitionFunction> pfuncCache, ConfSpaceInfo info, KStar.Sequence sequence) {

			// first check the cache
			PartitionFunction pfunc = pfuncCache.get(sequence);
			if (pfunc != null) {
				return pfunc;
			}

			// cache miss, need to compute the partition function

			// prune down to just this sequence
			// TODO: inherit any pruning from e.g. DEE
			PruningMatrix pmat = new PruningMatrix(info.confSpace, 0.0);
			for (SimpleConfSpace.Position pos : info.confSpace.positions) {
				String resType = sequence.get(pos.index);

				for (SimpleConfSpace.ResidueConf rc : pos.resConfs) {

					// if this RC doesn't match the sequence, prune it
					if (!rc.template.name.equals(resType)) {
						pmat.setOneBody(pos.index, rc.index, true);
					}
				}
			}

			// make the partition function
			pfunc = new SimplePartitionFunction(info.emat, pmat, confSearchFactory, info.confEcalc);
			pfunc.setReportProgress(settings.showPfuncProgress);
			pfunc.init(settings.epsilon);
			return pfunc;
		}

		@Override
		public void estimateScore() {

			kstarLower.estimateMore(settings.numConfsPerBatch);

			/* TEMP
			System.out.println("single sequence score: " + sequence);
			System.out.println("\tprotein:   " + kstarLower.protein.getStatus() + " " + kstarLower.protein.getValues().qstar.doubleValue());
			System.out.println("\tligand:    " + kstarLower.ligand.getStatus() + " " + kstarLower.ligand.getValues().qstar.doubleValue());
			System.out.println("\tcomplex:   " + kstarLower.complex.getStatus() + " " + kstarLower.complex.getValues().qstar.doubleValue());
			System.out.println("\tK* upper:  " + PartitionFunction.scoreToLog10String(kstarLower.calcKStarUpperBound()));
			System.out.println("\tK* lower:  " + PartitionFunction.scoreToLog10String(kstarLower.calcKStarLowerBound()));
			*/

			// update the score
			BigDecimal kstarScore = kstarLower.calcKStarUpperBound();
			if (kstarScore == null) {
				// give the lowest score, so this node polls last
				score = Double.NEGATIVE_INFINITY;
			} else {
				score = Math.log10(kstarScore.doubleValue());
			}
		}

		@Override
		public String toString() {
			return String.format("SingleSequenceNode[score=%12.6f, seq=%s]", score, sequence.toString());
		}
	}

	public static class ScoredSequence {

		public final KStar.Sequence sequence;
		public final BigDecimal kstarScore;

		public ScoredSequence(KStar.Sequence sequence, BigDecimal kstarScore) {
			this.sequence = sequence;
			this.kstarScore = kstarScore;
		}
	}

	public final ConfSpaceInfo protein;
	public final ConfSpaceInfo ligand;
	public final ConfSpaceInfo complex;
	public final KStar.ConfEnergyCalculatorFactory confEcalcFactory;
	public final ConfSearchFactory confSearchFactory;
	public final Settings settings;

	private final Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToProteinMap;
	private final Map<SimpleConfSpace.Position,SimpleConfSpace.Position> complexToLigandMap;

	private final Map<KStar.Sequence,PartitionFunction> proteinPfuncs;
	private final Map<KStar.Sequence,PartitionFunction> ligandPfuncs;
	private final Map<KStar.Sequence,PartitionFunction> complexPfuncs;

	public BBKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator ecalc, KStar.ConfEnergyCalculatorFactory confEcalcFactory, ConfSearchFactory confSearchFactory, Settings settings) {

		this.protein = new ConfSpaceInfo(ConfSpaceType.Protein, protein, confEcalcFactory.make(protein, ecalc));
		this.ligand = new ConfSpaceInfo(ConfSpaceType.Ligand, ligand, confEcalcFactory.make(ligand, ecalc));
		this.complex = new ConfSpaceInfo(ConfSpaceType.Complex, complex, confEcalcFactory.make(complex, ecalc));
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;

		complexToProteinMap = this.complex.confSpace.mapPositionsTo(this.protein.confSpace);
		complexToLigandMap = this.complex.confSpace.mapPositionsTo(this.ligand.confSpace);

		proteinPfuncs = new HashMap<>();
		ligandPfuncs = new HashMap<>();
		complexPfuncs = new HashMap<>();
	}

	private void computeEmatsIfNeeded() {
		if (protein.emat == null) {
			protein.calcEmat();
		}
		if (ligand.emat == null) {
			ligand.calcEmat();
		}
		if (complex.emat == null) {
			complex.calcEmat();
		}
	}

	public ScoredSequence calcWildType() {

		computeEmatsIfNeeded();

		System.out.println("computing K* score for wild type sequence to epsilon = " + settings.epsilon + " ...");
		SingleSequenceNode wtnode = new SingleSequenceNode(complex.makeWildTypeSequence());
		wtnode.kstarLower.compute();
		return new ScoredSequence(wtnode.sequence, wtnode.kstarLower.calcKStarScore());
	}

	public List<ScoredSequence> calcMutants() {

		computeEmatsIfNeeded();

		List<ScoredSequence> scoredSequences = new ArrayList<>();

		// start the BBK* tree with the root node
		PriorityQueue<Node> tree = new PriorityQueue<>();
		tree.add(new RootNode());

		// start searching the tree
		System.out.println("computing K* scores for the " + settings.numBestSequences + " best sequences to epsilon = " + settings.epsilon + " ...");
		scoredSequences.clear();
		while (!tree.isEmpty() && scoredSequences.size() < settings.numBestSequences) {

			// get the next node
			Node node = tree.poll();

			if (node instanceof SingleSequenceNode) {
				SingleSequenceNode ssnode = (SingleSequenceNode)node;

				switch (ssnode.kstarLower.getStatus()) {
					case Estimated:

						// sequence is finished, return it!
						BigDecimal kstarScore = ssnode.kstarLower.calcKStarScore();
						scoredSequences.add(new ScoredSequence(ssnode.sequence, kstarScore));

						// TODO: logging
						System.out.println(String.format("sequence: %s  kstar: %s  bound:[%s,%s]",
							ssnode.sequence,
							Math.log10(kstarScore.doubleValue()),
							PartitionFunction.scoreToLog10String(ssnode.kstarLower.calcKStarLowerBound()),
							PartitionFunction.scoreToLog10String(ssnode.kstarLower.calcKStarUpperBound())
						));

					break;
					case Estimating:

						// needs more estimation, catch-and-release
						ssnode.estimateScore();
						tree.add(ssnode);

					break;
					case Blocked:

						// from here on out, it's all blocked sequences
						// so it's ok to put them in the sorted order now
						scoredSequences.add(new ScoredSequence(ssnode.sequence, null));

						// TODO: logging
						System.out.println(String.format("unscorable sequence: %s  bound:[%s,%s]",
							ssnode.sequence,
							PartitionFunction.scoreToLog10String(ssnode.kstarLower.calcKStarLowerBound()),
							PartitionFunction.scoreToLog10String(ssnode.kstarLower.calcKStarUpperBound())
						));
				}

			} else if (node instanceof InteriorNode) {
				InteriorNode inode = (InteriorNode)node;

				// partial sequence, expand children
				List<Node> children = inode.makeChildren();
				// TODO: parallelize?
				for (Node child : children) {
					child.estimateScore();
				}
				tree.addAll(children);

				// TODO: do we ever catch-and-release to refine a multi-sequence bound?
			}
		}

		return scoredSequences;
	}
}
