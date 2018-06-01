package edu.duke.cs.osprey.gmec;

import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.SeqAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.*;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.log;


public class Comets {

	/**
	 * e.g., an unbound state, or a complex state
	 */
	public static class State {

		public static enum Type {

			/** desired binding */
			Positive(1.0),

			/** undesired binding */
			Negative(-1.0);

			public final double factor;

			private Type(double factor) {
				this.factor = factor;
			}
		}

		public final String name;
		public final SimpleConfSpace confSpace;
		public final Type type;

		public State(String name, SimpleConfSpace confSpace, Type type) {
			this.name = name;
			this.confSpace = confSpace;
			this.type = type;
		}
	}

	public static class WeightedState {

		public final State state;
		public final double weight;

		public WeightedState(State state, double weight) {
			this.state = state;
			this.weight = weight;
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

			public Builder addState(String name, SimpleConfSpace confSpace, State.Type type) {
				return addState(name, confSpace, type, 1.0);
			}

			public Builder addState(String name, SimpleConfSpace confSpace, State.Type type, double weight) {
				wstates.add(new WeightedState(new State(name, confSpace, type), weight));
				return this;
			}

			public Builder addPositiveState(String name, SimpleConfSpace confSpace) {
				return addState(name, confSpace, State.Type.Positive);
			}

			public Builder addNegativeState(String name, SimpleConfSpace confSpace) {
				return addState(name, confSpace, State.Type.Negative);
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
	}

	public static class SequenceInfo {

		public final Sequence sequence;
		public final ConfSearch.EnergiedConf minGMEC;

		public SequenceInfo(Sequence sequence, ConfSearch.EnergiedConf minGMEC) {
			this.sequence = sequence;
			this.minGMEC = minGMEC;
		}
	}


	public static class Builder {

		private final LME objective;
		private final List<LME> constraints = new ArrayList<>();

		public Builder(LME objective) {
			this.objective = objective;
		}

		public Builder addConstraint(LME constraint) {
			constraints.add(constraint);
			return this;
		}

		public Comets build() {
			return new Comets(objective, constraints);
		}
	}


	public final LME objective;
	public final List<LME> constraints;

	private final RTs rts;

	private Comets(LME objective, List<LME> constraints) {

		this.objective = objective;
		this.constraints = constraints;

		// collect all the states
		Set<State> statesSet = new LinkedHashSet<>();
		for (WeightedState wstate : objective.states) {
			statesSet.add(wstate.state);
		}
		for (LME constraint : constraints) {
			for (WeightedState wstate : constraint.states) {
				statesSet.add(wstate.state);
			}
		}
		List<State> states = new ArrayList<>(statesSet);

		// get the residue types for each state
		List<RTs> statesRTs = states.stream()
			.map(state -> new RTs(state.confSpace))
			.collect(Collectors.toList());

		// TEMP: dump all the RTs
		for (int i=0; i<states.size(); i++) {
			State state = states.get(i);
			log("%s %s", state.name, statesRTs.get(i).toString(state.confSpace));
		}

		// TODO: make sure they're all the same
		// except the bound state will have more single-restype positions
		// NEXTTIME: how to deal with that?

		rts = statesRTs.get(0);

		/* NOTES
			a state is a conf space with a weight
			all states must have the same sequence space
			states can have different flexibilities though
			LME: linear combination of state min GMEC energies
			objective function and constraints are all LMEs

			seq node g+h score is lower bound on objective function
		*/
	}

	public List<SequenceInfo> findBestSequences(int numSequences) {

		SeqAStarTree seqAstar = new SeqAStarTree.Builder(rts)
			.build();

		// TEMP
		throw new Error("not done yet");
	}
}
