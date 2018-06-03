package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.SeqAStarNode;
import edu.duke.cs.osprey.astar.seq.SeqAStarTree;
import edu.duke.cs.osprey.astar.seq.order.SequentialSeqAStarOrder;
import edu.duke.cs.osprey.astar.seq.scoring.NOPSeqAStarScorer;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.FragmentEnergies;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.tools.MathTools;
import edu.duke.cs.osprey.tools.UnpossibleError;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.formatBig;
import static edu.duke.cs.osprey.tools.Log.log;


// TODO: confDB
// TODO: sequence logging
// TODO: enforce constraints

public class Comets {

	/**
	 * e.g., an unbound state, or a complex state
	 */
	public static class State {

		public static class InitException extends RuntimeException {

			public InitException(State state, String name) {
				super(String.format("set %s for state %s before running", name, state.name));
			}
		}


		public final String name;
		public final SimpleConfSpace confSpace;

		public FragmentEnergies fragmentEnergies;
		public ConfEnergyCalculator confEcalc;
		public Function<RCs,ConfAStarTree> confTreeFactory;

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

	/** linear multi-state energy */
	public static class LME {

		public static class Builder {

			private double offset = 0.0;
			private final List<WeightedState> wstates = new ArrayList<>();

			public Builder setOffset(double val) {
				offset = val;
				return this;
			}

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
				SeqConfs.StateConfs stateConfs = confs.statesConfs.get(wstate.state);
				if (wstate.weight > 0) {
					val += wstate.weight*stateConfs.getObjectiveLowerBound();
				} else {
					val += wstate.weight*stateConfs.getObjectiveUpperBound();
				}
			}
			return val;
		}
	}

	/**
	 * implements A* heuristic for partially-defined sequences
	 * as described in COMETS paper, SI section B.2
	 */
	private class SeqHScorer implements SeqAStarScorer {

		SeqAStarNode.Assignments assignments = new SeqAStarNode.Assignments(rts.numMutablePos);
		MathTools.Optimizer opt = MathTools.Optimizer.Minimize;

		// collect positions from all states
		// ie, the union of positions where corresponding mutable positions across states are identified
		List<SimpleConfSpace.Position> allPositions = new ArrayList<>();

		// the states associated with each position
		// ie, mutable positions associate with all states
		// immutable positions associate only with their originating state
		List<List<WeightedState>> statesByAllPosition = new ArrayList<>();

		SeqHScorer() {
			allPositions.addAll(minimalState.confSpace.positions);
			for (int apos=0; apos<allPositions.size(); apos++) {
				statesByAllPosition.add(objective.states);
			}
			for (WeightedState wstate : objective.states) {
				allPositions.addAll(wstate.state.confSpace.immutablePositions);
				statesByAllPosition.add(Arrays.asList(wstate));
			}
		}

		@Override
		public double calc(SeqAStarNode node) {

			node.getAssignments(assignments);

			// TODO: in inner loops are independent of assignments, pre-compute them somehow?

			// sum over unassigned positions
			double score = objective.offset;
			for (int i1=0; i1<allPositions.size(); i1++) {
				SimpleConfSpace.Position pos1 = allPositions.get(i1);

				// optimize over res types at pos1
				double bestPos1Energy = opt.initDouble();
				for (String rt1 : getRTs(pos1, assignments)) {

					// sum over states
					double pos1Energy = 0.0;
					for (WeightedState wstate : statesByAllPosition.get(i1)) {

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
								for (String rt2 : getRTs(pos2, assignments)) {

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

						pos1Energy += bestRC1Energy;
					}

					bestPos1Energy = opt.opt(bestPos1Energy, pos1Energy);
				}

				score += bestPos1Energy;
			}

			// TEMP
			log("objective lower bound for seq %s = %.3f", makeSequence(node), score);

			return score;
		}

		List<String> getRTs(SimpleConfSpace.Position pos, SeqAStarNode.Assignments assignments) {
			// TODO: pre-compute this somehow?
			if (pos.hasMutations()) {
				Integer assignedRT = assignments.getAssignment(pos.mindex);
				if (assignedRT != null) {
					// use just the assigned res type
					return Arrays.asList(rts.getName(pos.mindex, assignedRT));
				} else {
					// use all the res types from the RTs collection
					return Arrays.stream(rts.indicesAt(pos.mindex))
						.mapToObj(rt -> rts.getName(pos.mindex, rt))
						.collect(Collectors.toList());
				}
			} else {
				// immutable position, use all the res types (should just be one)
				assert (pos.resTypes.size() == 1);
				return pos.resTypes;
			}
		}

		List<SimpleConfSpace.ResidueConf> getRCs(SimpleConfSpace.Position pos, String rt, State state) {
			// TODO: pre-compute this somehow?
			return pos.resConfs.stream()
				.filter(rc -> rc.template.name.equals(rt))
				.collect(Collectors.toList());
		}
	}

	/**
	 * storage for conf trees at each sequence node
	 * also tracks GMECs for each state
	 */
	private class SeqConfs {

		private class StateConfs {

			final State state;
			final ConfAStarTree confTree;

			ConfSearch.ScoredConf minScoreConf = null;
			ConfSearch.EnergiedConf minEnergyConf = null;

			ConfSearch.EnergiedConf gmec = null;

			StateConfs(SeqAStarNode seqNode, State state) {

				this.state = state;

				// make the conf tree
				RCs rcs = seqNode.makeSequence(state.confSpace).makeRCs();
				confTree = state.confTreeFactory.apply(rcs);
			}

			void refineBounds() {

				// already complete? no need to do more work
				if (gmec != null) {
					return;
				}

				// get the next conf
				ConfSearch.ScoredConf conf = confTree.nextConf();
				if (conf == null) {
					return;
				}

				// "refine" the lower bound
				if (minScoreConf == null) {
					minScoreConf = conf;
				}

				// refine the upper bound
				// TODO: parallelize: maybe do one batch in parallel each iteration?
				ConfSearch.EnergiedConf econf = state.confEcalc.calcEnergy(conf);
				if (minEnergyConf == null || econf.getEnergy() < minEnergyConf.getEnergy()) {
					minEnergyConf = econf;
				}

				// do we know the GMEC yet?
				if (conf.getScore() >= minEnergyConf.getEnergy()) {
					gmec = minEnergyConf;
				}
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

		final Map<State,StateConfs> statesConfs = new HashMap<>();

		SeqConfs(SeqAStarNode seqNode) {
			for (State state : states) {
				statesConfs.put(state, new StateConfs(seqNode, state));
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
		public double refineBounds() {

			// refine the GMEC bounds for each state
			for (State state : states) {
				statesConfs.get(state).refineBounds();
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

		public SequenceInfo(Sequence sequence, SeqConfs confs) {
			this.sequence = sequence;
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

		private double objectiveWindowSize = 10.0;
		private double objectiveWindowMax = 0.0;

		public Builder(LME objective) {
			this.objective = objective;
		}

		public Builder addConstraint(LME constraint) {
			constraints.add(constraint);
			return this;
		}

		/**
		 * The energy window is actually necessary for COMETS to finish in a reasonable
		 * amount of time in some edge cases. If the constraints restrict the sequence
		 * space to fewer than the desired number of best sequences, without an energy window,
		 * COMETS would enumerate every sequence while trying to find the desired number of
		 * sequences. This would effectively make COMETS linear in the number of possible
		 * sequences, which could be very slow.
		 *
		 * @param size limits the window relative to the objective value of the best sequence
		 * @param max absolute limit on the value of the objective function
		 */
		public Builder setObjectiveWindow(double size, double max) {
			objectiveWindowSize = size;
			objectiveWindowMax = max;
			return this;
		}

		public Comets build() {
			return new Comets(objective, constraints, objectiveWindowSize, objectiveWindowMax);
		}
	}


	public final LME objective;
	public final List<LME> constraints;
	public final double objectiveWindowSize;
	public final double objectiveWindowMax;

	public final List<State> states;

	private final State minimalState;
	private final RTs rts;

	private Comets(LME objective, List<LME> constraints, double objectiveWindowSize, double objectiveWindowMax) {

		this.objective = objective;
		this.constraints = constraints;
		this.objectiveWindowSize = objectiveWindowSize;
		this.objectiveWindowMax = objectiveWindowMax;

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

		// pick the state whos conf space has the fewest non-mutable positions as our canonical state
		minimalState = states.stream()
			.min(Comparator.comparing(state -> state.confSpace.positions.size()))
			.orElseThrow(() -> new UnpossibleError());

		// get the sequence space from the conf space
		rts = new RTs(minimalState.confSpace);
		log("sequence space has %s sequences\n%s", formatBig(rts.getNumSequences()), rts);

		// make sure all the other states match
		for (State state : states) {
			RTs stateRTs = new RTs(state.confSpace);
			if (!stateRTs.equals(rts)) {
				throw new IllegalArgumentException(String.format(
					"states have different sequence spaces:\n%s %s\n%s %s",
					minimalState, rts,
					state.name, stateRTs
				));
			}
		}
	}

	private Sequence makeSequence(SeqAStarNode node) {
		return node.makeSequence(minimalState.confSpace);
	}

	/**
	 * find the best sequences as ranked by the objective function
	 *
	 * searches all sequences within the objective window
	 */
	public List<SequenceInfo> findBestSequences(int numSequences) {

		// make sure all the states are fully configured
		for (State state : states) {
			state.checkConfig();
		}

		// start the A* search over sequences
		SeqAStarTree seqTree = new SeqAStarTree.Builder(rts)
			.setHeuristics(
				new SequentialSeqAStarOrder(),
				new NOPSeqAStarScorer(),
				new SeqHScorer()
			)
			.build();

		List<SequenceInfo> infos = new ArrayList<>();

		while (true) {

			// get the next sequence from the tree
			SeqAStarNode node = seqTree.nextLeafNode();
			if (node == null) {
				break;
			}

			// did we exhaust the sequences in the window?
			if (node.getScore() > objectiveWindowMax || (!infos.isEmpty() && node.getScore() > infos.get(0).objective + objectiveWindowSize)) {
				log("COMETS exiting early: exhausted all conformations in energy window");
				break;
			}

			// how are the conf trees here looking?
			SeqConfs confs = (SeqConfs)node.getData();
			if (confs == null) {

				// don't have them yet, make them
				confs = new SeqConfs(node);
				node.setData(confs);
			}

			// TEMP
			log("next sequence to check: %s   score: %.3f", makeSequence(node), node.getScore());

			// is this sequence finished already?
			if (confs.hasAllGMECs()) {

				// calc all the LMEs and output the sequence
				infos.add(new SequenceInfo(makeSequence(node), confs));

				// TEMP
				log("completed!");
				dump(node);

				// stop COMETS if we hit the desired number of sequences
				if (infos.size() >= numSequences) {
					break;
				}

			} else {

				// sequence needs more work, catch-and-release
				node.setHScore(confs.refineBounds());

				// TEMP
				log("refined!");
				dump(node);

				if (node.getScore() == Double.POSITIVE_INFINITY) {
					// constraint violated, prune this conf
					continue;
				}
				seqTree.add(node);
			}
		}

		if (infos.isEmpty()) {
			log("COMETS didn't find any sequences within the window that satisfy all the constraints.");
		} else {
			log("COMETS found the best %d within the window that satisfy all the constraints", infos.size());
		}

		return infos;
	}

	private void dump(SeqAStarNode node) {
		SeqConfs confs = (SeqConfs)node.getData();
		log("sequence %s", makeSequence(node));
		log("\tscore: %.6f   completed? %b", node.getScore(), confs.hasAllGMECs());
		log("\tobjective: %.6f", objective.calc(confs));
		for (LME constraint : constraints) {
			log("\tconstraint: %.3f", constraint.calc(confs));
		}
		for (SeqConfs.StateConfs stateConfs : confs.statesConfs.values()) {
			log("\tstate %-20s GMEC bounds [%8.3f,%8.3f]    found gmec? %b",
				stateConfs.state.name, stateConfs.minScoreConf.getScore(), stateConfs.minEnergyConf.getEnergy(), stateConfs.gmec != null
			);
		}
	}
}
