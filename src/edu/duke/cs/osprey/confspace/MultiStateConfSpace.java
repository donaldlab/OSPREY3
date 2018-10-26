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
	}

	public static class Builder {

		private final State sequencedState;
		private final List<State> states = new ArrayList<>();
		private final List<State> sequencedStates = new ArrayList<>();
		private final List<State> unsequencedStates = new ArrayList<>();

		/** adds the inital state that defines the sequence space */
		public Builder(String name, SimpleConfSpace confSpace) {
			sequencedState = new State(states.size(), sequencedStates.size(), -1, name, confSpace);
			states.add(sequencedState);
			sequencedStates.add(sequencedState);
		}

		public Builder addSequencedState(String name, SimpleConfSpace confSpace) {

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

		public Builder addState(String name, SimpleConfSpace confSpace) {

			// this conf space should have no mutants
			if (confSpace.seqSpace.hasMutants()) {
				throw new IllegalArgumentException("flexible conf space can't have any mutants");
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
}
