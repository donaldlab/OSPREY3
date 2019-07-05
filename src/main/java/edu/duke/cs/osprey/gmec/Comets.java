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

package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfSearchCache;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.astar.seq.SeqAStarTree;
import edu.duke.cs.osprey.astar.seq.order.SequentialSeqAStarOrder;
import edu.duke.cs.osprey.astar.seq.scoring.NOPSeqAStarScorer;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.formatBig;


/**
 * Implementation of the COMETS multi-state algorithm to predict protein sequence mutations that improve binding affinity.
 * {@cite Hallen2016 Mark A. Hallen, Bruce R. Donald, 2016.
 * Comets (Constrained Optimization of Multistate Energies by Tree Search): A provable and efficient protein design
 * algorithm to optimize binding affinity and specificity with respect to sequence.
 * Journal of Computational Biology, 23(5), 311-321.}.
 */
public class Comets {

	/**
	 * A state for a multi-state design
	 * (e.g., an unbound state, or a complex state)
	 *
	 * For COMETS, all states must share the same sequence space.
	 */
	public static class State {

		public static class InitException extends RuntimeException {

			public InitException(State state, String name) {
				super(String.format("set %s for state %s before running", name, state.name));
			}
		}


		/** A name for the state */
		public final String name;

		/** The conformation space for the state that describes the flexibility */
		public final SimpleConfSpace confSpace;

		public FragmentEnergies fragmentEnergies;
		public ConfEnergyCalculator confEcalc;
		public Function<RCs,ConfAStarTree> confTreeFactory;

		/**
		 * set this to a File if you want to use ConfDB with COMETS
		 */
		public File confDBFile = null;

		public State(String name, SimpleConfSpace confSpace) {
			this.name = name;
			this.confSpace = confSpace;
		}

		/**
		 * make sure the state is fully configured
		 */
		public void checkConfig() {
			if (fragmentEnergies == null) {
				throw new InitException(this, "fragmentEnergies");
			}
			if (confEcalc == null) {
				throw new InitException(this, "confEcalc");
			}
			if (confTreeFactory == null) {
				throw new InitException(this, "confTreeFactory");
			}
		}
	}

	public static class WeightedState {

		public final State state;
		public final double weight;

		public WeightedState(State state, double weight) {
			this.state = state;
			this.weight = weight;
		}

		public double getSingleEnergy(int pos, int rc) {
			return Math.abs(weight)*state.fragmentEnergies.getEnergy(pos, rc);
		}

		public double getPairEnergy(int pos1, int rc1, int pos2, int rc2) {
			return Math.abs(weight)*state.fragmentEnergies.getEnergy(pos1, rc1, pos2, rc2);
		}
	}

	/**
	 * A Linear Multi-state Energy
	 * (i.e., a weighted combination of state energies)
	 */
	public static class LME {

		public static class Builder {

			private double offset = 0.0;
			private final List<WeightedState> wstates = new ArrayList<>();

			public Builder setOffset(double val) {
				offset = val;
				return this;
			}

			/**
			 * Specify an upper bound constraint for this LME
			 * (i.e. set the offset to -val, since the constraint is satisfied iff the LME is < 0)
			 */
			public Builder constrainLessThan(double val) {
				return setOffset(-val);
			}

			public Builder addState(State state, double weight) {
				wstates.add(new WeightedState(state, weight));
				return this;
			}

			public LME build() {
				return new LME(offset, wstates);
			}
		}

		public final double offset;
		public final List<WeightedState> states;

		public LME(double offset, List<WeightedState> states) {
			this.offset = offset;
			this.states = states;
		}

		/**
		 * calculate a lower bound on the objective value for a fully-defined sequence node
		 * (i.e. A* heuristic from the COMETS paper SI section B.1)
		 *
		 * if the GMECs are known for all states, the lower bound will be tight (ie the true LME value)
		 */
		private double calc(SeqConfs confs) {
			double val = offset;
			for (WeightedState wstate : states) {
				StateConfs stateConfs = confs.statesConfs.get(wstate.state);
				if (wstate.weight > 0) {
					val += wstate.weight*stateConfs.getObjectiveLowerBound();
				} else {
					val += wstate.weight*stateConfs.getObjectiveUpperBound();
				}
			}
			return val;
		}

		/**
		 * calculate the exact value of the LME
		 */
		public double calc(Map<State,Double> stateEnergies) {
			double val = offset;
			for (WeightedState wstate : states) {
				val += wstate.weight*stateEnergies.get(wstate.state);
			}
			return val;
		}
	}

	/**
	 * implements A* heuristic for partially-defined sequences
	 * as described in COMETS paper, SI section B.2
	 */
	private class SeqHScorer implements SeqAStarScorer {

		MathTools.Optimizer opt = MathTools.Optimizer.Minimize;

		// collect conf space positions from all states
		List<SimpleConfSpace.Position> allPositions = new ArrayList<>();

		// the states associated with each position
		List<WeightedState> statesByAllPosition = new ArrayList<>();

		SeqHScorer() {
			for (WeightedState wstate : objective.states) {
				for (SimpleConfSpace.Position pos : wstate.state.confSpace.positions) {
					allPositions.add(pos);
					statesByAllPosition.add(wstate);
				}
			}
		}

		@Override
		public double calc(SeqAStarNode.Assignments assignments) {

			// TODO: if inner loops are independent of assignments, pre-compute them somehow?

			// sum over all positions
			double score = objective.offset;
			for (int i1=0; i1<allPositions.size(); i1++) {
				SimpleConfSpace.Position pos1 = allPositions.get(i1);

				// optimize over res types at pos1
				double bestPos1Energy = opt.initDouble();
				for (SeqSpace.ResType rt1 : getRTs(pos1, assignments)) {

					// sum over states (always just 1, due to how we build allPositions)
					WeightedState wstate = statesByAllPosition.get(i1);

					// min over RCs at (pos1,rt1,state)
					double bestRC1Energy = opt.initDouble();
					for (SimpleConfSpace.ResidueConf rc1 : getRCs(pos1, rt1, wstate.state)) {

						double rc1Energy = 0.0;

						// singles
						rc1Energy += wstate.getSingleEnergy(pos1.index, rc1.index);

						// pairs
						for (int i2=0; i2<pos1.index; i2++) {
							SimpleConfSpace.Position pos2 = wstate.state.confSpace.positions.get(i2);

							// min over RTs at pos2
							double bestRT2Energy = opt.initDouble();
							for (SeqSpace.ResType rt2 : getRTs(pos2, assignments)) {

								// min over RCs at (pos2,rt2,state)
								double bestRC2Energy = opt.initDouble();
								for (SimpleConfSpace.ResidueConf rc2 : getRCs(pos2, rt2, wstate.state)) {

									double rc2Energy = wstate.getPairEnergy(pos1.index, rc1.index, pos2.index, rc2.index);

									bestRC2Energy = opt.opt(bestRC2Energy, rc2Energy);
								}

								bestRT2Energy = opt.opt(bestRT2Energy, bestRC2Energy);
							}

							rc1Energy += bestRT2Energy;
						}

						bestRC1Energy = opt.opt(bestRC1Energy, rc1Energy);
					}

					bestPos1Energy = opt.opt(bestPos1Energy, bestRC1Energy);
				}

				score += bestPos1Energy;
			}

			return score;
		}

		List<SeqSpace.ResType> getRTs(SimpleConfSpace.Position confPos, SeqAStarNode.Assignments assignments) {

			// TODO: pre-compute this somehow?

			// map the conf pos to a sequence pos
			SeqSpace.Position seqPos = seqSpace.getPosition(confPos.resNum);
			if (seqPos != null) {

				Integer assignedRT = assignments.getAssignment(seqPos.index);
				if (assignedRT != null) {
					// use just the assigned res type
					return Collections.singletonList(seqPos.resTypes.get(assignedRT));
				} else {
					// use all the res types at the pos
					return seqPos.resTypes;
				}

			} else {

				// immutable position, use all the res types (should just be one)
				assert (confPos.resTypes.size() == 1);

				// use the null value to signal there's no res type here
				return Collections.singletonList(null);
			}
		}

		List<SimpleConfSpace.ResidueConf> getRCs(SimpleConfSpace.Position pos, SeqSpace.ResType rt, State state) {
			// TODO: pre-compute this somehow?
			if (rt != null) {
				// mutable pos, grab the RCs that match the RT
				return pos.resConfs.stream()
					.filter(rc -> rc.template.name.equals(rt.name))
					.collect(Collectors.toList());
			} else {
				// immutable pos, use all the RCs
				return pos.resConfs;
			}
		}
	}

	private class ConfDBs extends ConfDB.DBs {

		public Map<State,ConfDB.ConfTable> tables = new HashMap<>();

		public ConfDBs() {

			// open the DBs
			for (State state : states) {
				add(state.confSpace, state.confDBFile);
			}

			// make the tables
			for (State state : states) {
				ConfDB confdb = get(state.confSpace);
				if (confdb != null) {
					tables.put(state, confdb.new ConfTable("COMETS"));
				}
			}
		}
	}

	/**
	 * essentially, an iterative mini-GMEC-finder for a sequence and a state
	 */
	private static class StateConfs {

		private static class Key {

			final Sequence sequence;
			final State state;

			Key(Sequence sequence, State state) {
				this.sequence = sequence;
				this.state = state;
			}

			@Override
			public int hashCode() {
				return HashCalculator.combineHashes(
					sequence.hashCode(),
					state.hashCode()
				);
			}

			@Override
			public boolean equals(Object other) {
				return other instanceof Key && equals((Key)other);
			}

			public boolean equals(Key other) {
				return this.sequence.equals(other.sequence)
					&& this.state == other.state;
			}
		}

		final State state;
		final Sequence sequence; // filtered to the state

		ConfSearch confTree = null;
		ConfSearch.ScoredConf minScoreConf = null;
		ConfSearch.EnergiedConf minEnergyConf = null;
		ConfSearch.EnergiedConf gmec = null;

		List<ConfSearch.ScoredConf> confs = new ArrayList<>();

		StateConfs(Sequence sequence, State state, ConfSearchCache confTrees) {

			this.state = state;
			this.sequence = sequence;

			// make the conf tree
			RCs rcs = sequence.makeRCs(state.confSpace);
			confTree = confTrees.make(() -> state.confTreeFactory.apply(rcs));
		}

		void refineBounds(ConfDB.ConfTable confTable) {

			// already complete? no need to do more work
			if (gmec != null) {
				return;
			}

			// get the next batch of confs
			confs.clear();
			for (int i=0; i<state.confEcalc.tasks.getParallelism(); i++) {

				// get the next conf
				ConfSearch.ScoredConf conf = confTree.nextConf();
				if (conf == null) {
					break;
				}

				// "refine" the lower bound
				if (minScoreConf == null) {
					minScoreConf = conf;
				}

				confs.add(conf);
			}

			// no more confs? nothing to do
			if (confs.isEmpty()) {
				return;
			}

			for (ConfSearch.ScoredConf conf : confs) {

				// refine the upper bound
				state.confEcalc.calcEnergyAsync(conf, confTable, econf -> {

					// NOTE: don't need to lock here, since the main thread is waiting

					if (minEnergyConf == null || econf.getEnergy() < minEnergyConf.getEnergy()) {
						minEnergyConf = econf;
					}
				});
			}

			state.confEcalc.tasks.waitForFinish();

			// do we know the GMEC yet?
			ConfSearch.ScoredConf maxScoreConf = confs.get(confs.size() - 1);
			if (maxScoreConf.getScore() >= minEnergyConf.getEnergy()) {
				gmec = minEnergyConf;

				// release the resources used by the conf tree (could be a lot of memory)
				confTree = null;
			}

			confs.clear();
		}

		double getObjectiveLowerBound() {
			if (gmec != null) {
				return gmec.getEnergy();
			}
			return minScoreConf.getScore();
		}

		double getObjectiveUpperBound() {
			if (gmec != null) {
				return gmec.getEnergy();
			}
			return minEnergyConf.getEnergy();
		}
	}


	/**
	 * storage for conf trees at each sequence node
	 * also tracks GMECs for each state
	 */
	private class SeqConfs {

		final Map<State,StateConfs> statesConfs = new HashMap<>();

		SeqConfs(SeqAStarNode seqNode) {

			for (State state : states) {

				Sequence sequence = seqNode.makeSequence(seqSpace)
					.filter(state.confSpace.seqSpace);

				// get state confs from the global cache if possible, otherwise make a new one
				StateConfs.Key key = new StateConfs.Key(sequence, state);
				StateConfs stateConfs = stateConfsCache.get(key);
				if (stateConfs == null) {
					stateConfs = new StateConfs(sequence, state, confTrees);
					stateConfsCache.put(key, stateConfs);
				}

				statesConfs.put(state, stateConfs);
			}
		}

		boolean hasAllGMECs() {
			for (State state : states) {
				if (statesConfs.get(state).gmec == null) {
					return false;
				}
			}
			return true;
		}

		/**
		 * implements A* heuristic for fully-defined sequences
		 * as described in COMETS paper, SI section B.1
		 *
		 * returns the new score for the seqeunce node
		 *
		 * also flags that GMECs are found, when applicable
		 */
		public double refineBounds(ConfDBs confDBs) {

			// refine the GMEC bounds for each state
			for (State state : states) {
				statesConfs.get(state).refineBounds(confDBs.tables.get(state));
			}

			// if any constraints are violated, score the node +inf,
			// so it never gets enumerated again by A*
			for (LME constraint : constraints) {
				if (constraint.calc(this) > 0) {
					return Double.POSITIVE_INFINITY;
				}
			}

			// evaluate the objective function
			return objective.calc(this);
		}
	}

	public class SequenceInfo {

		public final Sequence sequence;
		public final Map<State,ConfSearch.EnergiedConf> GMECs = new HashMap<>();
		public final double objective;
		public final Map<LME,Double> constraints = new HashMap<>();

		public SequenceInfo(SeqAStarNode node, SeqConfs confs) {
			this.sequence = node.makeSequence(seqSpace);
			for (State state : states) {
				GMECs.put(state, confs.statesConfs.get(state).gmec);
			}
			objective = Comets.this.objective.calc(confs);
			for (LME constraint : Comets.this.constraints) {
				constraints.put(constraint, constraint.calc(confs));
			}
		}
	}


	public static class Builder {

		private final LME objective;
		private final List<LME> constraints = new ArrayList<>();

		/**
		 * The energy window is actually necessary for COMETS to finish in a reasonable
		 * amount of time in some edge cases. If the constraints restrict the sequence
		 * space to fewer than the desired number of best sequences, without an energy window,
		 * COMETS would enumerate every sequence while trying to find the desired number of
		 * sequences. This would effectively make COMETS run in time that is linear in the
		 * number of possible sequences, which could be very slow.
		 *
		 * The window size effectively adds an additional constraint that the difference between
		 * the objective function value and the objective function value of the best sequence
		 * must be below this value.
		 */
		private double objectiveWindowSize = 10.0;

		/**
		 * The window max effectively adds an additional constraint that the objective function value
		 * must be below this value.
		 */
		private double objectiveWindowMax = 0.0;

		/** The maximum number of simultaneous residue mutations to consider for each sequence mutant */
		private int maxSimultaneousMutations = 1;

		/**
		 * The minimum number of conformation trees to keep in memory at once.
		 *
		 * Defauls to null, which means keep all trees in memory at once.
		 *
		 * Positive values keep at least that number of trees in memory at once
		 * (where more frequently-used trees are preferred over less frequently-used trees),
		 * and any trees above that number may be deleted by the JVM's garbage collector.
		 *
		 * Any tree that is delete must be re-instantiated and reset to its
		 * previous state before it can be used again, which incurrs a performance penalty.
		 *
		 * Deleting less frequently-used trees and re-instantiating them when needed,
		 * along with using bounded-memory implementations of A* search, allows design
		 * algortihms to run within constant memory, while maintaining good performance.
		 */
		private Integer minNumConfTrees = null;

		private boolean printToConsole = true;

		/** File to which to log sequences as they are found */
		private File logFile = null;

		public Builder(LME objective) {
			this.objective = objective;
		}

		public Builder addConstraint(LME constraint) {
			constraints.add(constraint);
			return this;
		}

		public Builder setObjectiveWindowSize(double val) {
			objectiveWindowSize = val;
			return this;
		}

		public Builder setObjectiveWindowMax(double val) {
			objectiveWindowMax = val;
			return this;
		}

		public Builder setMaxSimultaneousMutations(int val) {
			maxSimultaneousMutations = val;
			return this;
		}

		public Builder setMinNumConfTrees(Integer val) {
			minNumConfTrees = val;
			return this;
		}

		public Builder setPrintToConsole(boolean val) {
			printToConsole = val;
			return this;
		}

		public Builder setLogFile(File val) {
			logFile = val;
			return this;
		}

		public Comets build() {
			return new Comets(objective, constraints, objectiveWindowSize, objectiveWindowMax, maxSimultaneousMutations, minNumConfTrees, printToConsole, logFile);
		}
	}


	public final LME objective;
	public final List<LME> constraints;
	public final double objectiveWindowSize;
	public final double objectiveWindowMax;
	public final int maxSimultaneousMutations;
	public final Integer minNumConfTrees;
	public final boolean printToConsole;
	public final File logFile;

	public final List<State> states;
	public final SeqSpace seqSpace;

	private final Map<StateConfs.Key,StateConfs> stateConfsCache = new HashMap<>();
	private final ConfSearchCache confTrees;

	private Comets(LME objective, List<LME> constraints, double objectiveWindowSize, double objectiveWindowMax, int maxSimultaneousMutations, Integer minNumConfTrees, boolean printToConsole, File logFile) {

		this.objective = objective;
		this.constraints = constraints;
		this.objectiveWindowSize = objectiveWindowSize;
		this.objectiveWindowMax = objectiveWindowMax;
		this.maxSimultaneousMutations = maxSimultaneousMutations;
		this.minNumConfTrees = minNumConfTrees;
		this.printToConsole = printToConsole;
		this.logFile = logFile;

		// collect all the states from the objective,constraints
		Set<State> statesSet = new LinkedHashSet<>();
		for (WeightedState wstate : objective.states) {
			statesSet.add(wstate.state);
		}
		for (LME constraint : constraints) {
			for (WeightedState wstate : constraint.states) {
				statesSet.add(wstate.state);
			}
		}
		states = new ArrayList<>(statesSet);
		if (states.isEmpty()) {
			throw new IllegalArgumentException("COMETS found no states");
		}

		// get the sequence space from the conf spaces
		seqSpace = SeqSpace.union(
			states.stream()
				.map(state -> state.confSpace.seqSpace)
				.collect(Collectors.toList())
		);

		confTrees = new ConfSearchCache(minNumConfTrees);

		log("sequence space has %s sequences\n%s", formatBig(new RTs(seqSpace).getNumSequences()), seqSpace);
	}

	/**
	 * find the best sequences as ranked by the objective function
	 *
	 * searches all sequences within the objective window
	 */
	public List<SequenceInfo> findBestSequences(int numSequences) {

		// reset any previous state
		stateConfsCache.clear();

		// make sure all the states are fully configured
		for (State state : states) {
			state.checkConfig();
		}

		// start the A* search over sequences
		SeqAStarTree seqTree = new SeqAStarTree.Builder(new RTs(seqSpace))
			.setHeuristics(
				new SequentialSeqAStarOrder(),
				new NOPSeqAStarScorer(),
				new SeqHScorer()
			)
			.setNumMutable(maxSimultaneousMutations)
			.build();

		List<SequenceInfo> infos = new ArrayList<>();

		log("\nCOMETS searching for the %d best sequences among %s with up to %d simultaneous mutations ...",
			numSequences,
			formatBig(new RTs(seqSpace).getNumSequences()),
			maxSimultaneousMutations
		);
		log("(up to objective function value %.6f kcal/mol, or +%.6f kcal/mol relative to the best sequence)",
			objectiveWindowMax,
			objectiveWindowSize
		);
		log("");

		// open the ConfDBs if needed
		try (ConfDBs confDBs = new ConfDBs()) {

			while (true) {

				// get the next sequence from the tree
				SeqAStarNode node = seqTree.nextLeafNode();
				if (node == null) {
					break;
				}

				// did we exhaust the sequences in the window?
				if (node.getScore() > objectiveWindowMax || (!infos.isEmpty() && node.getScore() > infos.get(0).objective + objectiveWindowSize)) {
					log("\nCOMETS exiting early: exhausted all conformations in energy window");
					break;
				}

				// how are the conf trees here looking?
				SeqConfs confs = (SeqConfs)node.getData();
				if (confs == null) {

					log("Discovered promising sequence: %s   objective lower bound: %12.6f",
						node.makeSequence(seqSpace),
						node.getScore()
					);

					// don't have them yet, make them
					confs = new SeqConfs(node);
					node.setData(confs);
				}

				// is this sequence finished already?
				if (confs.hasAllGMECs()) {

					// calc all the LMEs and output the sequence
					SequenceInfo info = new SequenceInfo(node, confs);
					infos.add(info);
					reportSequence(infos.size() == 1, info);

					// stop COMETS if we hit the desired number of sequences
					if (infos.size() >= numSequences) {
						break;
					}

				} else {

					// sequence needs more work, catch-and-release
					node.setHScore(confs.refineBounds(confDBs));

					if (node.getScore() == Double.POSITIVE_INFINITY) {
						// constraint violated, prune this conf
						continue;
					}

					// add the sequence back to the tree
					seqTree.add(node);
				}
			}
		}

		log("");
		if (infos.isEmpty()) {
			log("COMETS didn't find any sequences within the window that satisfy all the constraints.");
		} else {
			log("COMETS found the best %d within the window that satisfy all the constraints", infos.size());
		}

		return infos;
	}

	private void log(String msg, Object ... args) {
		if (printToConsole) {
			edu.duke.cs.osprey.tools.Log.log(msg, args);
		}
	}

	private void reportSequence(boolean isFirstSequence, SequenceInfo info) {

		if (printToConsole) {
			int cellSize = info.sequence.calcCellSize();
			log("\nSequence calculation complete: %s    objective: %12.6f   %s\n                               %s",
				info.sequence.toString(Sequence.Renderer.ResNum, cellSize),
				info.objective,
				info.sequence.isWildType() ? "Wild-type" : "",
				info.sequence.toString(Sequence.Renderer.ResTypeMutations, cellSize)
			);
			for (State state : states) {
				log("\tState: %-20s    GMEC Energy: %12.6f",
					state.name,
					info.GMECs.get(state).getEnergy()
				);
			}
		}

		if (logFile != null) {

			// build the log record in TSV format

			// make the header if needed
			StringBuilder header = null;
			if (isFirstSequence) {
				header = new StringBuilder();
				for (SeqSpace.Position pos : seqSpace.positions) {
					if (pos.index > 0) {
						header.append("\t");
					}
					header.append(pos.resNum);
				}
				header.append("\tObjective");
				for (int i=0; i<constraints.size(); i++) {
					header.append("\tConstraint ");
					header.append(i);
				}
				for (State state : states) {
					header.append("\t");
					header.append(state.name);
					header.append(" GMEC Energy\t");
					header.append(state.name);
					header.append(" GMEC Conf");
				}
			}

			// make the sequence entry
			StringBuilder buf = new StringBuilder();
			for (SeqSpace.Position pos : seqSpace.positions) {
				if (pos.index > 0) {
					buf.append("\t");
				}
				buf.append(info.sequence.get(pos).mutationName());
			}
			buf.append("\t");
			buf.append(String.format("%.6f", info.objective));
			for (int i=0; i<constraints.size(); i++) {
				buf.append("\t");
				buf.append(String.format("%.6f", info.constraints.get(constraints.get(i))));
			}
			for (State state : states) {
				ConfSearch.EnergiedConf gmec = info.GMECs.get(state);
				buf.append("\t");
				buf.append(String.format("%.6f", gmec.getEnergy()));
				buf.append("\t");
				buf.append(Conf.toString(gmec.getAssignments()));
			}

			// write to the log file, but don't keep the file open
			try (FileWriter out = new FileWriter(logFile, !isFirstSequence)) {

				if (header != null) {
					out.write(header.toString());
					out.write("\n");
				}
				out.write(buf.toString());
				out.write("\n");

			} catch (IOException ex) {

				// don't take down the whole job if we can't log the sequence for some reason
				// just report the error, and dump info to the console
				ex.printStackTrace(System.err);
				if (header != null) {
					System.err.println(header);
				}
				System.err.println(buf);
			}
		}
	}

	// for debugging
	@SuppressWarnings("unused")
	private void dump(SeqAStarNode node) {
		SeqConfs confs = (SeqConfs)node.getData();
		log("sequence %s", node.makeSequence(seqSpace));
		log("\tscore: %.6f   completed? %b", node.getScore(), confs.hasAllGMECs());
		log("\tobjective: %.6f", objective.calc(confs));
		for (LME constraint : constraints) {
			log("\tconstraint: %.3f", constraint.calc(confs));
		}
		for (StateConfs stateConfs : confs.statesConfs.values()) {
			log("\tstate %-20s GMEC bounds [%8.3f,%8.3f]    found gmec? %b",
				stateConfs.state.name, stateConfs.minScoreConf.getScore(), stateConfs.minEnergyConf.getEnergy(), stateConfs.gmec != null
			);
		}
	}
}
