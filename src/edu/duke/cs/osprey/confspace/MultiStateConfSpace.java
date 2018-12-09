package edu.duke.cs.osprey.confspace;


import java.util.*;


public class MultiStateConfSpace {

	public static class State {

		public final int index;
		public final int sequencedIndex;
		public final int unsequencedIndex;
		public final String name;
		public final SimpleConfSpace confSpace;
		public final boolean isSequenced;

		private State(int index, int sequencedIndex, int unsequencedIndex, String name, SimpleConfSpace confSpace) {
			this.index = index;
			this.sequencedIndex = sequencedIndex;
			this.unsequencedIndex = unsequencedIndex;
			this.name = name;
			this.confSpace = confSpace;
			this.isSequenced = sequencedIndex >= 0;
		}

		@Override
		public int hashCode() {
			return index;
		}

		@Override
		public boolean equals(Object other) {
			return other instanceof State && equals((State)other);
		}

		public boolean equals(State other) {
			return this.index == other.index;
		}
	}

	public static class Builder {

		private final State sequencedState;
		private final List<State> states = new ArrayList<>();
		private final List<State> sequencedStates = new ArrayList<>();
		private final List<State> unsequencedStates = new ArrayList<>();
		private Double fixedResiduesCutoff = null;

		/** adds the inital mutable state that defines the sequence space */
		public Builder(String name, SimpleConfSpace confSpace) {
			sequencedState = new State(states.size(), sequencedStates.size(), -1, name, confSpace);
			states.add(sequencedState);
			sequencedStates.add(sequencedState);
		}

		public Builder addMutableState(String name, SimpleConfSpace confSpace) {

			// make sure this conf space matches the sequence space
			if (!confSpace.seqSpace.equals(sequencedState.confSpace.seqSpace)) {
				throw new IllegalArgumentException(String.format(
					"sequence space for state \"%s\" doesn't match state \"%s\"\nexpected:%s\nobserved:%s",
					name,
					sequencedState.name,
					sequencedState.confSpace.seqSpace,
					confSpace.seqSpace
				));
			}

			State state = new State(states.size(), sequencedStates.size(), -1, name, confSpace);
			states.add(state);
			sequencedStates.add(state);

			return this;
		}

		public Builder addUnmutableState(String name, SimpleConfSpace confSpace) {

			// this conf space should have no mutants
			if (confSpace.seqSpace.hasMutants()) {
				throw new IllegalArgumentException("unmutable conf space can't have any mutants");
			}

			State state = new State(states.size(), -1, unsequencedStates.size(), name, confSpace);
			states.add(state);
			unsequencedStates.add(state);

			return this;
		}

		public MultiStateConfSpace build() {
			return new MultiStateConfSpace(sequencedState.confSpace.seqSpace, states, sequencedStates, unsequencedStates);
		}
	}


	public final SeqSpace seqSpace;
	public final List<State> states;
	public final List<State> sequencedStates;
	public final List<State> unsequencedStates;

	private MultiStateConfSpace(SeqSpace seqSpace, List<State> states, List<State> sequencedStates, List<State> unsequencedStates) {
		this.seqSpace = seqSpace;
		this.states = Collections.unmodifiableList(states);
		this.sequencedStates = sequencedStates;
		this.unsequencedStates = unsequencedStates;
	}

	public State getState(String name) {
		return states.stream()
			.filter(state -> state.name.equals(name))
			.findFirst()
			.orElseThrow(() -> new NoSuchElementException("no state with name " + name));
	}

	/**
	 * A linear multi-state free energy.
	 * Used as an optimization objective function for many design algorithms.
	 */
	public class LMFE {

		public final MultiStateConfSpace confSpace = MultiStateConfSpace.this;

		private final Map<State,Double> weightsByState;

		private LMFE(Map<State,Double> weightsByState) {
			this.weightsByState = Collections.unmodifiableMap(weightsByState);
		}

		public Set<State> states() {
			return weightsByState.keySet();
		}

		public double getWeight(State state) {
			return weightsByState.get(state);
		}
	}

	public class LMFEBuilder {

		private final Map<State,Double> weightsByState = new HashMap<>();

		public LMFEBuilder add(State state, double weight) {
			weightsByState.put(state, weight);
			return this;
		}

		public LMFEBuilder add(String stateName, double weight) {
			weightsByState.put(getState(stateName), weight);
			return this;
		}

		public LMFEBuilder addPositive(State state) {
			return add(state, 1.0);
		}

		public LMFEBuilder addPositive(String stateName) {
			return addPositive(getState(stateName));
		}

		public LMFEBuilder addNegative(State state) {
			return add(state, -1.0);
		}

		public LMFEBuilder addNegative(String stateName) {
			return addNegative(getState(stateName));
		}

		public LMFE build() {
			return new LMFE(weightsByState);
		}
	}


	/** build a linear multi-state free energy */
	public LMFEBuilder lmfe() {
		return new LMFEBuilder();
	}
}
