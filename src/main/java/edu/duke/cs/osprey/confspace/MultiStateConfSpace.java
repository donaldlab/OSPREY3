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

package edu.duke.cs.osprey.confspace;


import edu.duke.cs.osprey.tools.MathTools.DoubleBounds;

import java.util.*;
import java.util.function.Function;


/**
 * A composite conformation space composed of other conformation spaces.
 * Each conformation space is a "state" and is associated with a descriptive name, like "protein" or "ligand".
 *
 * Each state can be mutable and have an associated sequence space, or be unmutable.
 * The first conformation space must be mutable, and it defines the sequence space
 * for the entire multi-state conformation space.
 */
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

		public Double getWeight(State state) {
			return weightsByState.get(state);
		}

		public DoubleBounds[] collectFreeEnergies(Function<State,DoubleBounds> f) {
			DoubleBounds[] boundsByState = new DoubleBounds[confSpace.states.size()];
			for (State state : states()) {
				boundsByState[state.index] = f.apply(state);
			}
			return boundsByState;
		}

		public LMFECalculator calc() {
			return new LMFECalculator(this);
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

	public class LMFECalculator {

		public final LMFE lmfe;

		public final DoubleBounds bounds = new DoubleBounds(0, 0);

		public LMFECalculator(LMFE lmfe) {
			this.lmfe = lmfe;
		}

		public LMFECalculator add(State state, DoubleBounds freeEnergy) {

			// get the weight, if any
			Double weight = lmfe.getWeight(state);
			if (weight == null) {
				return this;
			}

			// if either bound is [+inf,+inf], assume the result is [+inf,+inf]
			if ((bounds.lower == Double.POSITIVE_INFINITY && bounds.upper == Double.POSITIVE_INFINITY)
				|| (freeEnergy.lower == Double.POSITIVE_INFINITY && freeEnergy.upper == Double.POSITIVE_INFINITY)) {

				bounds.lower = Double.POSITIVE_INFINITY;
				bounds.upper = Double.POSITIVE_INFINITY;

			// otherwise, do the usual arithmetic
			} else {

				if (weight < 0) {
					bounds.lower += freeEnergy.upper*weight;
					bounds.upper += freeEnergy.lower*weight;
				} else if (weight > 0) {
					bounds.lower += freeEnergy.lower*weight;
					bounds.upper += freeEnergy.upper*weight;
				}
			}

			return this;
		}

		public LMFECalculator addAll(DoubleBounds[] freeEnergiesByState) {

			for (State state : lmfe.confSpace.states) {
				add(state, freeEnergiesByState[state.index]);
			}

			return this;
		}
	}
}
