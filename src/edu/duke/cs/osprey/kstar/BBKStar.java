package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.kstar.KStar.ConfSearchFactory;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.kstar.pfunc.SimplePartitionFunction;

import java.io.File;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.*;
import java.util.function.Function;

public class BBKStar {

	// *sigh* Java makes this stuff so verbose to do...
	// Kotlin would make this so much easier
	public static class Settings {

		public static class Builder {

			private double epsilon = 0.683;
			private int maxSimultaneousMutations = Integer.MAX_VALUE;
			private boolean showPfuncProgress = false;
			private int numBestSequences = 1;
			private int numConfsPerBatch = 8;
			private KStarScoreWriter.Writers scoreWriters = new KStarScoreWriter.Writers();

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

			public Builder addScoreWriter(KStarScoreWriter val) {
				scoreWriters.add(val);
				return this;
			}

			public Builder addScoreConsoleWriter(KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToConsole(val));
			}

			public Builder addScoreConsoleWriter() {
				return addScoreConsoleWriter(new KStarScoreWriter.Formatter.SequenceKStarPfuncs());
			}

			public Builder addScoreFileWriter(File file, KStarScoreWriter.Formatter val) {
				return addScoreWriter(new KStarScoreWriter.ToFile(file, val));
			}

			public Builder addScoreFileWriter(File file) {
				return addScoreFileWriter(file, new KStarScoreWriter.Formatter.Log());
			}

			public Settings build() {
				return new Settings(epsilon, maxSimultaneousMutations, showPfuncProgress, numBestSequences, numConfsPerBatch, scoreWriters);
			}
		}

		public final double epsilon;
		public final int maxSimultaneousMutations;
		public final boolean showPfuncProgress;
		public final int numBestSequences;
		public final int numConfsPerBatch;
		public final KStarScoreWriter.Writers scoreWriters;

		public Settings(double epsilon, int maxSimultaneousMutations, boolean dumpPfuncConfs, int numBestSequences, int numConfsPerBatch, KStarScoreWriter.Writers scoreWriters) {
			this.epsilon = epsilon;
			this.maxSimultaneousMutations = maxSimultaneousMutations;
			this.showPfuncProgress = dumpPfuncConfs;
			this.numBestSequences = numBestSequences;
			this.numConfsPerBatch = numConfsPerBatch;
			this.scoreWriters = scoreWriters;
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
		public final ConfEnergyCalculator rigidConfEcalc;
		public final ConfEnergyCalculator minimizingConfEcalc;

		public EnergyMatrix rigidNegatedEmat = null;
		public EnergyMatrix minimizedEmat = null;

		public ConfSpaceInfo(ConfSpaceType type, SimpleConfSpace confSpace, ConfEnergyCalculator rigidConfEcalc, ConfEnergyCalculator minimizingConfEcalc) {
			this.type = type;
			this.confSpace = confSpace;
			this.rigidConfEcalc = rigidConfEcalc;
			this.minimizingConfEcalc = minimizingConfEcalc;
		}

		public void calcEmatsIfNeeded() {
			if (rigidNegatedEmat == null) {
				rigidNegatedEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc)
					.build()
					.calcEnergyMatrix();

				// negate the rigid energy matrix, so we can convince A* to return scores
				// in descending order instead of ascending order
				// (of course, we'll have to negate the scores returned by A* too)
				rigidNegatedEmat.negate();
			}
			if (minimizedEmat == null) {
				minimizedEmat = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc)
					.build()
					.calcEnergyMatrix();
			}
		}

		public RCs makeRCs(KStar.Sequence sequence) {
			return new RCs(confSpace, (pos, resConf) -> {

				// if there's an assignment here, only keep matching RCs
				String resType = sequence.get(pos.index);
				if (resType != null) {
					return resConf.template.name.equals(resType);
				}

				// otherwise, keep everything
				return true;
			});
		}
	}

	public static enum PfuncsStatus {
		Estimating,
		Estimated,
		Blocked
	}

	private abstract class Node implements Comparable<Node> {

		public final KStar.Sequence sequence;
		public final KStar.Sequence proteinSequence;
		public final KStar.Sequence ligandSequence;

		/** for comparing in the tree, higher is first */
		public double score;

		protected Node(KStar.Sequence sequence) {

			this.sequence = sequence;

			// split complex sequence into protein/ligand sequences
			proteinSequence = KStar.Sequence.makeWildType(BBKStar.this.protein.confSpace);
			ligandSequence = KStar.Sequence.makeWildType(BBKStar.this.ligand.confSpace);
			for (SimpleConfSpace.Position pos : BBKStar.this.complex.confSpace.positions) {

				SimpleConfSpace.Position proteinPos = complexToProteinMap.get(pos);
				if (proteinPos != null) {
					proteinSequence.set(proteinPos.index, sequence.get(pos.index));
				}

				SimpleConfSpace.Position ligandPos = complexToLigandMap.get(pos);
				if (ligandPos != null) {
					ligandSequence.set(ligandPos.index, sequence.get(pos.index));
				}
			}

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

		public MultiSequenceNode(KStar.Sequence sequence) {
			super(sequence);
		}

		public List<Node> makeChildren() {

			List<Node> children = new ArrayList<>();

			List<SimpleConfSpace.Position> positions = complex.confSpace.positions;

			// pick the next design position
			// TODO: dynamic A*?
			int posIndex;
			for (posIndex = 0; posIndex< sequence.size(); posIndex++) {
				if (sequence.get(posIndex) == null) {
					break;
				}
			}
			if (posIndex >= sequence.size()) {
				throw new IllegalStateException("no design positions left to choose");
			}
			SimpleConfSpace.Position assignPos = positions.get(posIndex);

			// get the possible assignments
			Set<String> resTypes = new HashSet<>(assignPos.resFlex.resTypes);

			// add wild-type option if mutations are limited
			if (settings.maxSimultaneousMutations < positions.size()) {
				resTypes.add(assignPos.resFlex.wildType);
			}

			// for each assignment...
			for (String resType : resTypes) {

				// update the sequence with this assignment
				KStar.Sequence s = sequence.makeWithAssignment(assignPos.index, resType);

				if (s.isFullyAssigned()) {

					// fully assigned, make single sequence node
					children.add(new SingleSequenceNode(s));

				} else if (s.countMutations(complex.confSpace) == settings.maxSimultaneousMutations) {

					// mutation limit reached, fill unassigned positions with wild-type
					s.fillWildType(complex.confSpace);
					children.add(new SingleSequenceNode(s));

				} else {

					// still partial sequence, make multi-sequence node
					children.add(new MultiSequenceNode(s));
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

			// to compute lower bounds on pfuncs, we'll use the upper bound calculator,
			// but feed it upper-bound scores in order of greatest upper bound
			// instead of the usual lower-bound scores in order of smallest lower bound
			// which means we need to feed A* with a rigid, negated energy matrix
			// and then negate all the scores returned by A*
			Function<ConfSearch,ConfSearch> astarNegater = (confSearch) -> new ConfSearch() {
				@Override
				public ScoredConf nextConf() {
					ScoredConf conf = confSearch.nextConf();
					conf.setScore(-conf.getScore());
					return conf;
				}

				@Override
				public BigInteger getNumConformations() {
					return confSearch.getNumConformations();
				}

				@Override
				public List<ScoredConf> nextConfs(double maxEnergy) {
					throw new UnsupportedOperationException("some lazy programmer didn't implement this =P");
				}
			};

			// calculate a lower bound on the protein partition function
			BigDecimal proteinLowerBound;
			{
				SimplePartitionFunction.UpperBoundCalculator calc = new SimplePartitionFunction.UpperBoundCalculator(
					astarNegater.apply(confSearchFactory.make(protein.rigidNegatedEmat, protein.makeRCs(proteinSequence)))
				);
				calc.run(numConfs);
				proteinLowerBound = calc.totalBound;
			}

			// if the first few conf upper bound scores (for the pfunc lower bound) are too high,
			// then the K* upper bound is also too high
			if (proteinLowerBound.compareTo(BigDecimal.ZERO) == 0) {
				score = Double.POSITIVE_INFINITY;
				return;
			}

			// calculate a lower bound on the ligand partition function
			BigDecimal ligandLowerBound;
			{
				SimplePartitionFunction.UpperBoundCalculator calc = new SimplePartitionFunction.UpperBoundCalculator(
					astarNegater.apply(confSearchFactory.make(ligand.rigidNegatedEmat, ligand.makeRCs(ligandSequence)))
				);
				calc.run(numConfs);
				ligandLowerBound = calc.totalBound;
			}

			// if the first few conf upper bound scores (for the pfunc lower bound) are too high,
			// then the K* upper bound is also too high
			if (ligandLowerBound.compareTo(BigDecimal.ZERO) == 0) {
				score = Double.POSITIVE_INFINITY;
				return;
			}

			// NOTE: to compute upper bounds on pfuncs,
			// we'll use the upper bound calculator in the usual way

			// calculate an upper bound on the complex partition function
			BigDecimal complexUpperBound;
			{
				SimplePartitionFunction.UpperBoundCalculator calc = new SimplePartitionFunction.UpperBoundCalculator(
					confSearchFactory.make(complex.minimizedEmat, complex.makeRCs(sequence))
				);
				calc.run(numConfs);
				complexUpperBound = calc.totalBound;
			}

			// if there are no low-energy confs, then make this node last in the queue
			if (ligandLowerBound.compareTo(BigDecimal.ZERO) == 0) {
				score = Double.NEGATIVE_INFINITY;
				return;
			}

			// compute the node score
			score = Math.log10(complexUpperBound
				.divide(proteinLowerBound, RoundingMode.HALF_UP)
				.divide(ligandLowerBound, RoundingMode.HALF_UP)
				.doubleValue());
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

		public SingleSequenceNode(KStar.Sequence sequence) {
			super(sequence);

			// make the partition functions
			this.protein = makePfunc(proteinPfuncs, BBKStar.this.protein, proteinSequence);
			this.ligand = makePfunc(ligandPfuncs, BBKStar.this.ligand, ligandSequence);
			this.complex = makePfunc(complexPfuncs, BBKStar.this.complex, sequence);
		}

		private PartitionFunction makePfunc(Map<KStar.Sequence,PartitionFunction> pfuncCache, ConfSpaceInfo info, KStar.Sequence sequence) {

			// first check the cache
			PartitionFunction pfunc = pfuncCache.get(sequence);
			if (pfunc != null) {
				return pfunc;
			}

			// cache miss, need to compute the partition function

			// get RCs for just this sequence
			RCs rcs = new RCs(info.confSpace, (pos, resConf) -> {
				return resConf.template.name.equals(sequence.get(pos.index));
			});

			// make the partition function
			pfunc = new SimplePartitionFunction(confSearchFactory.make(info.minimizedEmat, rcs), info.minimizingConfEcalc);
			pfunc.setReportProgress(settings.showPfuncProgress);
			pfunc.init(settings.epsilon);
			return pfunc;
		}

		@Override
		public void estimateScore() {

			// refine the pfuncs if needed
			if (protein.getStatus().canContinue()) {
				protein.compute(settings.numConfsPerBatch);
			}
			if (ligand.getStatus().canContinue()) {
				ligand.compute(settings.numConfsPerBatch);
			}
			if (complex.getStatus().canContinue()) {
				complex.compute(settings.numConfsPerBatch);
			}

			// update the score
			KStarScore kstarScore = makeKStarScore();
			if (kstarScore.upperBound == null) {
				// give the lowest score, so this node polls last
				score = Double.NEGATIVE_INFINITY;
			} else {
				score = Math.log10(kstarScore.upperBound.doubleValue());
			}
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

	public static class ScoredSequence {

		public final KStar.Sequence sequence;
		public final KStarScore score;

		public ScoredSequence(KStar.Sequence sequence, KStarScore score) {
			this.sequence = sequence;
			this.score = score;
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

	// TODO: caching these will keep lots of A* trees in memory. is that a problem?
	private final Map<KStar.Sequence,PartitionFunction> proteinPfuncs;
	private final Map<KStar.Sequence,PartitionFunction> ligandPfuncs;
	private final Map<KStar.Sequence,PartitionFunction> complexPfuncs;

	public BBKStar(SimpleConfSpace protein, SimpleConfSpace ligand, SimpleConfSpace complex, EnergyCalculator rigidEcalc, EnergyCalculator minimizingEcalc, KStar.ConfEnergyCalculatorFactory confEcalcFactory, ConfSearchFactory confSearchFactory, Settings settings) {

		this.protein = new ConfSpaceInfo(
			ConfSpaceType.Protein,
			protein,
			confEcalcFactory.make(protein, rigidEcalc),
			confEcalcFactory.make(protein, minimizingEcalc)
		);
		this.ligand = new ConfSpaceInfo(
			ConfSpaceType.Ligand,
			ligand,
			confEcalcFactory.make(ligand, rigidEcalc),
			confEcalcFactory.make(ligand, minimizingEcalc)
		);
		this.complex = new ConfSpaceInfo(
			ConfSpaceType.Complex,
			complex,
			confEcalcFactory.make(complex, rigidEcalc),
			confEcalcFactory.make(complex, minimizingEcalc)
		);
		this.confEcalcFactory = confEcalcFactory;
		this.confSearchFactory = confSearchFactory;
		this.settings = settings;

		complexToProteinMap = this.complex.confSpace.mapPositionsTo(this.protein.confSpace);
		complexToLigandMap = this.complex.confSpace.mapPositionsTo(this.ligand.confSpace);

		proteinPfuncs = new HashMap<>();
		ligandPfuncs = new HashMap<>();
		complexPfuncs = new HashMap<>();
	}

	public List<ScoredSequence> run() {

		protein.calcEmatsIfNeeded();
		ligand.calcEmatsIfNeeded();
		complex.calcEmatsIfNeeded();

		List<ScoredSequence> scoredSequences = new ArrayList<>();

		// start the BBK* tree with the root node
		PriorityQueue<Node> tree = new PriorityQueue<>();
		tree.add(new MultiSequenceNode(new KStar.Sequence(complex.confSpace.positions.size())));

		// start searching the tree
		System.out.println("computing K* scores for the " + settings.numBestSequences + " best sequences to epsilon = " + settings.epsilon + " ...");
		while (!tree.isEmpty() && scoredSequences.size() < settings.numBestSequences) {

			// get the next node
			Node node = tree.poll();

			if (node instanceof SingleSequenceNode) {
				SingleSequenceNode ssnode = (SingleSequenceNode)node;

				// single-sequence node
				switch (ssnode.getStatus()) {
					case Estimated:

						// sequence is finished, return it!
						reportSequence(ssnode, scoredSequences);

					break;
					case Estimating:

						// needs more estimation, catch-and-release
						ssnode.estimateScore();
						tree.add(ssnode);

					break;
					case Blocked:

						// from here on out, it's all blocked sequences
						// so it's ok to put them in the sorted order now
						reportSequence(ssnode, scoredSequences);
				}

			} else if (node instanceof MultiSequenceNode) {
				MultiSequenceNode msnode = (MultiSequenceNode)node;

				// partial sequence, expand children
				List<Node> children = msnode.makeChildren();
				// TODO: parallelize the multi-sequence node scoring here?
				for (Node child : children) {
					child.estimateScore();
				}
				tree.addAll(children);
			}
		}

		if (scoredSequences.size() < settings.numBestSequences) {
			if (tree.isEmpty()) {
				// all is well, we just don't have that many sequences in the design
				System.out.println("Tried to find " + settings.numBestSequences + " sequences,"
					+ " but design flexibility only allowed " + scoredSequences.size() + " sequences.");
			} else {
				throw new Error("BBK* ended, but the tree isn't empty and we didn't return enough sequences. This is a bug.");
			}
		}

		return scoredSequences;
	}

	private void reportSequence(SingleSequenceNode ssnode, List<ScoredSequence> scoredSequences) {

		KStarScore kstarScore = ssnode.makeKStarScore();
		scoredSequences.add(new ScoredSequence(ssnode.sequence, kstarScore));

		settings.scoreWriters.writeScore(new KStarScoreWriter.ScoreInfo(
			scoredSequences.size(),
			settings.numBestSequences,
			ssnode.sequence,
			complex.confSpace,
			kstarScore
		));
	}
}
