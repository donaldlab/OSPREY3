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


import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.ConfSearchCache;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.SeqAStarTree;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.astar.seq.order.SequentialSeqAStarOrder;
import edu.duke.cs.osprey.astar.seq.scoring.NOPSeqAStarScorer;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.*;
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
 * basically COMETS, but with thermodynamic ensembles
 *
 * at least, that's the goal
 */
public class MSKStar {

	/**
	 * A state for a multi-state design
	 * (e.g., an unbound state, or a complex state)
	 */
	public static class State {

		public PartitionFunctionFactory pfuncFactory;

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
		 * set this to a File if you want to use ConfDB with MSK*
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
			if (pfuncFactory == null) {
				throw new InitException(this, "pfuncFactory");
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
	 * A Linear Multi-state Free Energy
	 * (i.e., a weighted combination of state free energies)
	 */
	public static class LMFE {

		public static class Builder {

			private double offset = 0.0;
			private final List<WeightedState> wstates = new ArrayList<>();

			public Builder setOffset(double val) {
				offset = val;
				return this;
			}

			/**
			 * Specify an upper bound constraint for this LMFE
			 * (i.e. set the offset to -val, since the constraint is satisfied iff the LMFE is < 0)
			 */
			public Builder constrainLessThan(double val) {
				return setOffset(-val);
			}

			public Builder addState(State state, double weight) {
				wstates.add(new WeightedState(state, weight));
				return this;
			}

			public LMFE build() {
				return new LMFE(offset, wstates);
			}
		}

		public final double offset;
		public final List<WeightedState> states;

		public LMFE(double offset, List<WeightedState> states) {
			this.offset = offset;
			this.states = states;
		}

		/**
		 * calculate bounds on the objective value for a fully-defined sequence node
		 */
		private MathTools.DoubleBounds calc(SeqConfs confs) {
			MathTools.DoubleBounds bound = new MathTools.DoubleBounds(offset, offset);
			for (WeightedState wstate : states) {
				StateConfs stateConfs = confs.statesConfs.get(wstate.state);
				if (wstate.weight > 0) {
					bound.lower += wstate.weight*stateConfs.freeEnergyBounds.lower;
					bound.upper += wstate.weight*stateConfs.freeEnergyBounds.upper;
				} else {
					bound.lower += wstate.weight*stateConfs.freeEnergyBounds.upper;
					bound.upper += wstate.weight*stateConfs.freeEnergyBounds.lower;
				}
			}
			return bound;
		}

		/**
		 * calculate the LMFE bounds based on bounds on free energies
		 */
		public MathTools.DoubleBounds calc(Map<State,MathTools.DoubleBounds> stateFreeEnergies) {
			MathTools.DoubleBounds bound = new MathTools.DoubleBounds(offset, offset);
			for (WeightedState wstate : states) {
				MathTools.DoubleBounds freeEnergyBounds = stateFreeEnergies.get(wstate.state);
				if (wstate.weight > 0) {
					bound.lower += wstate.weight*freeEnergyBounds.lower;
					bound.upper += wstate.weight*freeEnergyBounds.upper;
				} else {
					bound.lower += wstate.weight*freeEnergyBounds.upper;
					bound.upper += wstate.weight*freeEnergyBounds.lower;
				}
			}
			return bound;
		}
	}

	/**
	 * implements A* heuristic for partially-defined sequences
	 */
	private class SeqHScorer implements SeqAStarScorer {

		// TODO: make configurable
		private static final int upperBatchSize = 1000;
		private static final int numLowerBatches = 1;

		final ConfDBs confDBs;

		SeqHScorer(ConfDBs confDBs) {
			this.confDBs = confDBs;
		}

		BoltzmannCalculator bcalc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

		@Override
		public double calc(SeqAStarNode.Assignments assignments) {

			double lowerBound = 0.0;

			for (WeightedState wstate : objective.states) {

				// make a conf tree for all the sequences descibed by the partial assignments
				RCs rcs = assignments.makeRCs(seqSpace, wstate.state.confSpace);
				ConfAStarTree confTree = wstate.state.confTreeFactory.apply(rcs);

				if (wstate.weight > 0) {

					// compute an upper bound on the pfunc value of the multi-sequence conf space
					UpperBoundCalculator ubcalc = new UpperBoundCalculator(confTree, rcs.getNumConformations());
					ubcalc.run(upperBatchSize);

					// add the bound to our estimate of the objective LMFE
					lowerBound += wstate.weight*bcalc.freeEnergy(ubcalc.totalBound);

				} else {

					// compute a lower bound on the pfunc value of the multi-sequence conf space
					LowerBoundCalculator lbcalc = new LowerBoundCalculator(confTree, wstate.state.confEcalc);
					lbcalc.confTable = confDBs.tables.get(wstate.state);
					for (int i=0; i<numLowerBatches; i++) {
						lbcalc.run(wstate.state.confEcalc.tasks.getParallelism());
					}

					// add the weighted bound to our estimate of the objective LMFE
					lowerBound += wstate.weight*bcalc.freeEnergy(lbcalc.weightedEnergySum);
				}
			}

			return lowerBound;
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
					tables.put(state, confdb.new ConfTable("MSK*"));
				}
			}
		}
	}

	/**
	 * essentially, an iterative partition function calculator for a sequence and a state
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

		final MathTools.DoubleBounds freeEnergyBounds = new MathTools.DoubleBounds();

		PartitionFunction pfunc = null;
		PartitionFunction.Result pfuncResult = null;

		StateConfs(Sequence sequence, State state, double epsilon, ConfDB.ConfTable confTable, ConfSearchCache confTrees) {

			this.state = state;
			this.sequence = sequence;

			// init pfunc calculation
			RCs rcs = sequence.makeRCs(state.confSpace);
			pfunc = state.pfuncFactory.makePartitionFunctionFor(rcs, rcs.getNumConformations(), epsilon);

			if (confTable != null) {
				PartitionFunction.WithConfTable.setOrThrow(pfunc, confTable);
			}
		}

		void refineBounds() {

			// already done? no need to be an over-achiever
			if (pfuncResult != null) {
				return;
			}

			// update the bounds
			pfunc.compute(state.confEcalc.tasks.getParallelism());
			pfunc.getValues().calcFreeEnergyBounds(freeEnergyBounds);

			// are we there yet?
			if (!pfunc.getStatus().canContinue()) {
				pfuncResult = pfunc.makeResult();

				// release the resources used by the pfunc (e.g., the memory for the A* tree)
				pfunc = null;
			}
		}

		boolean isRefinementComplete() {
			return pfuncResult != null;
		}
	}


	/**
	 * storage for pfunc calculators at each sequence node
	 */
	private class SeqConfs {

		final Map<State,StateConfs> statesConfs = new HashMap<>();

		SeqConfs(SeqAStarNode seqNode, ConfDBs confDBs) {

			for (State state : states) {

				Sequence sequence = seqNode.makeSequence(seqSpace)
					.filter(state.confSpace.seqSpace);

				// get state confs from the global cache if possible, otherwise make a new one
				StateConfs.Key key = new StateConfs.Key(sequence, state);
				StateConfs stateConfs = stateConfsCache.get(key);
				if (stateConfs == null) {
					stateConfs = new StateConfs(sequence, state, epsilon, confDBs.tables.get(state), confTrees);
					stateConfsCache.put(key, stateConfs);
				}

				statesConfs.put(state, stateConfs);
			}
		}

		boolean isRefinementComplete() {
			for (State state : states) {
				if (!statesConfs.get(state).isRefinementComplete()) {
					return false;
				}
			}
			return true;
		}

		/**
		 * implements A* heuristic for fully-defined sequences
		 *
		 * returns bounds on the objective function for the seqeunce node
		 */
		public MathTools.DoubleBounds refineBounds() {

			// refine the GMEC bounds for each state
			for (State state : states) {
				statesConfs.get(state).refineBounds();
			}

			// if any constraints are violated, score the node +inf,
			// so it never gets enumerated again by A*
			for (LMFE constraint : constraints) {
				if (constraint.calc(this).lower > 0) {
					return new MathTools.DoubleBounds(Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY);
				}
			}

			// evaluate the objective function
			return objective.calc(this);
		}
	}

	public class SequenceInfo {

		public final Sequence sequence;
		public final Map<State,PartitionFunction.Result> pfuncResults = new HashMap<>();
		public final MathTools.DoubleBounds objective;
		public final Map<LMFE,MathTools.DoubleBounds> constraints = new HashMap<>();

		public SequenceInfo(SeqAStarNode node, SeqConfs confs) {
			this.sequence = node.makeSequence(seqSpace);
			for (State state : states) {
				pfuncResults.put(state, confs.statesConfs.get(state).pfuncResult);
			}
			objective = MSKStar.this.objective.calc(confs);
			for (LMFE constraint : MSKStar.this.constraints) {
				constraints.put(constraint, constraint.calc(confs));
			}
		}
	}

	public static class Builder {

		private final LMFE objective;
		private final List<LMFE> constraints = new ArrayList<>();

		/**
		 * Epsilon value for partition function approximation.
		 */
		private double epsilon = 0.68;

		/**
		 * The energy window is actually necessary for MSK* to finish in a reasonable
		 * amount of time in some edge cases. If the constraints restrict the sequence
		 * space to fewer than the desired number of best sequences, without an energy window,
		 * MSK* would enumerate every sequence while trying to find the desired number of
		 * sequences. This would effectively make MSK* run in time that is linear in the
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

		/** Temporary MARKStar fields because this interface doesn't cleanly support
		 * non-ConfSearch PartitionFunction implementations...
		 */
		public EnergyMatrix rigidEmat;
		public EnergyMatrix minimizingEmat;


		public Builder(LMFE objective) {
			this.objective = objective;
		}

		public Builder addConstraint(LMFE constraint) {
			constraints.add(constraint);
			return this;
		}

		public Builder setEpsilon(double val) {
			epsilon = val;
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

		public MSKStar build() {
		    MSKStar mskStar = new MSKStar(objective, constraints, epsilon, objectiveWindowSize, objectiveWindowMax, maxSimultaneousMutations, minNumConfTrees, printToConsole, logFile);

			return mskStar;
		}
	}


	public final LMFE objective;
	public final List<LMFE> constraints;
	public final double epsilon;
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


	private MSKStar(LMFE objective, List<LMFE> constraints, double epsilon, double objectiveWindowSize, double objectiveWindowMax, int maxSimultaneousMutations, Integer minNumConfTrees, boolean printToConsole, File logFile) {

		this.objective = objective;
		this.constraints = constraints;
		this.epsilon = epsilon;
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
		for (LMFE constraint : constraints) {
			for (WeightedState wstate : constraint.states) {
				statesSet.add(wstate.state);
			}
		}
		states = new ArrayList<>(statesSet);
		if (states.isEmpty()) {
			throw new IllegalArgumentException("MSK* found no states");
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

		List<SequenceInfo> infos = new ArrayList<>();

		// open the ConfDBs if needed
		try (ConfDBs confDBs = new ConfDBs()) {

			// start the A* search over sequences
			SeqAStarTree seqTree = new SeqAStarTree.Builder(new RTs(seqSpace))
				.setHeuristics(
					new SequentialSeqAStarOrder(),
					new NOPSeqAStarScorer(),
					new SeqHScorer(confDBs)
				)
				.setNumMutable(maxSimultaneousMutations)
				.build();

			log("\nMSK* searching for the %d best sequences among %s with up to %d simultaneous mutations ...",
				numSequences,
				formatBig(new RTs(seqSpace).getNumSequences()),
				maxSimultaneousMutations
			);
			log("(up to objective function value %.6f kcal/mol, or +%.6f kcal/mol relative to the best sequence)",
				objectiveWindowMax,
				objectiveWindowSize
			);
			log("");

			while (true) {

				// get the next sequence from the tree
				SeqAStarNode node = seqTree.nextLeafNode();
				if (node == null) {
					break;
				}

				// did we exhaust the sequences in the window?
				if (node.getScore() > objectiveWindowMax || (!infos.isEmpty() && node.getScore() > infos.get(0).objective.upper + objectiveWindowSize)) {
					log("\nMSK* exiting early: exhausted all conformations in energy window");
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
					confs = new SeqConfs(node, confDBs);
					node.setData(confs);
				}

				// is this sequence finished already?
				if (confs.isRefinementComplete()) {

					// calc all the LMEs and output the sequence
					SequenceInfo info = new SequenceInfo(node, confs);
					infos.add(info);
					reportSequence(infos.size() == 1, info);

					// stop MSK* if we hit the desired number of sequences
					if (infos.size() >= numSequences) {
						break;
					}

				} else {

					// sequence needs more work, catch-and-release
					node.setHScore(confs.refineBounds().lower);

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
			log("MSK* didn't find any sequences within the window that satisfy all the constraints.");
		} else {
			log("MSK* found the best %d within the window that satisfy all the constraints", infos.size());
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
			log("\nSequence calculation complete: %s    objective: %s   %s\n                               %s",
				info.sequence.toString(Sequence.Renderer.ResNum, cellSize),
				info.objective.toString(4, 9),
				info.sequence.isWildType() ? "Wild-type" : "",
				info.sequence.toString(Sequence.Renderer.ResTypeMutations, cellSize)
			);
			for (State state : states) {
				log("\tState: %-20s    Free Energy: %s",
					state.name,
					info.pfuncResults.get(state).values.calcFreeEnergyBounds().toString(4, 9)
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
				header.append("\tObjective Min");
				header.append("\tObjective Max");
				for (int i=0; i<constraints.size(); i++) {
					header.append("\tConstraint ");
					header.append(i);
					header.append(" Min\tConstraint ");
					header.append(i);
					header.append(" Max");
				}
				for (State state : states) {
					header.append("\t");
					header.append(state.name);
					header.append(" Free Energy Min\t");
					header.append(state.name);
					header.append(" Free Energy Max");
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
			buf.append(String.format("\t%.4f\t%.4f", info.objective.lower, info.objective.upper));
			for (int i=0; i<constraints.size(); i++) {
				MathTools.DoubleBounds constraint = info.constraints.get(constraints.get(i));
				buf.append(String.format("\t%.4f\t%.4f", constraint.lower, constraint.upper));
			}
			for (State state : states) {
				PartitionFunction.Result pfuncResult = info.pfuncResults.get(state);
				buf.append(String.format("\t%.4f\t%.4f", pfuncResult.values.calcFreeEnergyLowerBound(), pfuncResult.values.calcFreeEnergyUpperBound()));
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
}
