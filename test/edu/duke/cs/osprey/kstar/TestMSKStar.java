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

import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.Test;

import java.util.*;
import java.util.stream.Collectors;


public class TestMSKStar {

	private static class Problem {

		final MSKStar mskstar;
		final ForcefieldParams ffparams;

		public Problem(MSKStar mskstar, ForcefieldParams ffparams) {
			this.mskstar = mskstar;
			this.ffparams = ffparams;
		}
	}

	private static Problem make2RL0PPI(boolean boundedMemory) {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0();
		final double epsilon = 0.95;

		MSKStar.State protein = new MSKStar.State("Protein", confSpaces.protein);
		MSKStar.State ligand = new MSKStar.State("Ligand", confSpaces.ligand);
		MSKStar.State complex = new MSKStar.State("Complex", confSpaces.complex);

		MSKStar.LMFE objective = new MSKStar.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();
		MSKStar mskstar = new MSKStar.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowSize(100.0) // need a big window to get all the sequences
			.setObjectiveWindowMax(100.0)
			.setMinNumConfTrees(boundedMemory ? 5 : null)
			.build();

		initStates(mskstar.states, confSpaces.ffparams, boundedMemory);

		return new Problem(mskstar, confSpaces.ffparams);
	}

	private static void check2RL0PPI(MSKStar mskstar) {

		List<MSKStar.SequenceInfo> sequences = mskstar.findBestSequences(25);

		// check the sequences (values collected with e = 0.01 and 64 digits precision)
		assertSequence(sequences, "phe asp glu thr phe lys ile thr",   -20.947212,   -20.935187);
		assertSequence(sequences, "TYR asp glu thr phe lys ile thr",   -15.810058,   -15.797771);
		assertSequence(sequences, "ALA asp glu thr phe lys ile thr",   -16.105190,   -16.094644);
		assertSequence(sequences, "VAL asp glu thr phe lys ile thr",   -17.246337,   -17.235435);
		assertSequence(sequences, "ILE asp glu thr phe lys ile thr",   -17.589388,   -17.576462);
		assertSequence(sequences, "LEU asp glu thr phe lys ile thr",   -16.887802,   -16.875529);
		assertSequence(sequences, "phe GLU glu thr phe lys ile thr",   -19.241928,   -19.228557);
		assertSequence(sequences, "phe asp ASP thr phe lys ile thr",   -18.347308,   -18.335335);
		assertSequence(sequences, "phe asp glu SER phe lys ile thr",   -21.863663,   -21.847602);
		assertSequence(sequences, "phe asp glu ASN phe lys ile thr",   -22.056226,   -22.040606);
		assertSequence(sequences, "phe asp glu GLN phe lys ile thr",   -22.480570,   -22.465154);
		assertSequence(sequences, "phe asp glu thr TYR lys ile thr",   -20.557311,   -20.545232);
		assertSequence(sequences, "phe asp glu thr ALA lys ile thr",   -19.361427,   -19.350788);
		assertSequence(sequences, "phe asp glu thr VAL lys ile thr",   -20.040577,   -20.029484);
		assertSequence(sequences, "phe asp glu thr ILE lys ile thr",   -20.522378,   -20.511615);
		assertSequence(sequences, "phe asp glu thr LEU lys ile thr",   -19.998842,   -19.987029);
		assertSequence(sequences, "phe asp glu thr phe ASP ile thr",   -14.787522,   -14.770574);
		assertSequence(sequences, "phe asp glu thr phe GLU ile thr",   -13.699482,   -13.689671);
		assertSequence(sequences, "phe asp glu thr phe lys ALA thr",   -19.129036,   -19.117935);
		assertSequence(sequences, "phe asp glu thr phe lys VAL thr",   -19.915711,   -19.903392);
		assertSequence(sequences, "phe asp glu thr phe lys LEU thr",    -5.361432,    -5.353948);
		assertSequence(sequences, "phe asp glu thr phe lys PHE thr",    39.735126,    39.741235);
		assertSequence(sequences, "phe asp glu thr phe lys TYR thr",    42.064914,    42.071749);
		assertSequence(sequences, "phe asp glu thr phe lys ile SER",   -21.729311,   -21.716258);
		assertSequence(sequences, "phe asp glu thr phe lys ile ASN",   -20.413312,   -20.403102);

		assertThat(sequences.size(), is(25));
		assertSequenceOrder(sequences);
	}

	private static Problem make2RL0OnlyOneMutant() {

		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0OnlyOneMutant();
		final double epsilon = 0.95;

		MSKStar.State protein = new MSKStar.State("Protein", confSpaces.protein);
		MSKStar.State ligand = new MSKStar.State("Ligand", confSpaces.ligand);
		MSKStar.State complex = new MSKStar.State("Complex", confSpaces.complex);

		MSKStar.LMFE objective = new MSKStar.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();
		MSKStar mskstar = new MSKStar.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(1)
			.setObjectiveWindowSize(100.0) // need a big window to get all the sequences
			.setObjectiveWindowMax(100.0)
			.build();

		initStates(mskstar.states, confSpaces.ffparams, false);

		return new Problem(mskstar, confSpaces.ffparams);
	}

	private static void check2RL0OnlyOneMutant(MSKStar mskstar) {

		List<MSKStar.SequenceInfo> sequences = mskstar.findBestSequences(1);

		// check the sequences (values collected with e = 0.01 and 64 digits precision)
		assertSequence(sequences, "VAL",    -4.541643,    -4.541643);
		
		assertThat(sequences.size(), is(1));
		assertSequenceOrder(sequences);
	}

	private static Problem make2RL0SpaceWithoutWildType() {
	
		TestKStar.ConfSpaces confSpaces = TestKStar.make2RL0SpaceWithoutWildType();
		final double epsilon = 0.95;

		MSKStar.State protein = new MSKStar.State("Protein", confSpaces.protein);
		MSKStar.State ligand = new MSKStar.State("Ligand", confSpaces.ligand);
		MSKStar.State complex = new MSKStar.State("Complex", confSpaces.complex);

		MSKStar.LMFE objective = new MSKStar.LMFE.Builder()
			.addState(complex, 1.0)
			.addState(protein, -1.0)
			.addState(ligand, -1.0)
			.build();
		MSKStar mskstar = new MSKStar.Builder(objective)
			.setEpsilon(epsilon)
			.setMaxSimultaneousMutations(2)
			.setObjectiveWindowSize(100.0) // need a big window to get all the sequences
			.setObjectiveWindowMax(100.0)
			.build();

		initStates(mskstar.states, confSpaces.ffparams, false);

		return new Problem(mskstar, confSpaces.ffparams);
	}

	private static void check2RL0WithoutWildType(MSKStar mskstar) {

		List<MSKStar.SequenceInfo> sequences = mskstar.findBestSequences(2);

		// check the sequences (values collected with e = 0.01 and 64 digits precision)
		assertSequence(sequences, "thr VAL",    -4.541643,    -4.541643);
		assertSequence(sequences, "VAL VAL",    -4.674453,    -4.674453);
		
		assertThat(sequences.size(), is(2));
		assertSequenceOrder(sequences);
	}

	private static void initStates(List<MSKStar.State> states, ForcefieldParams ffparams, boolean boundedMemory) {

		List<SimpleConfSpace> confSpaceList = states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaceList, ffparams)
			.setParallelism(Parallelism.makeCpu(8))
			.build()
		) {

			for (MSKStar.State state : states) {

				// how should we define energies of conformations?
				state.confEcalc = new ConfEnergyCalculator.Builder(state.confSpace, ecalc)
					.setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(state.confSpace, ecalc)
						.build()
						.calcReferenceEnergies()
					)
					.build();

				// calc energy matrix
				EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(state.confEcalc)
					.build()
					.calcEnergyMatrix();
				state.fragmentEnergies = emat;

				// how should confs be ordered and searched?
				state.confTreeFactory = (rcs) -> new ConfAStarTree.Builder(emat, rcs)
					.setMaxNumNodes(boundedMemory ? 100000L : null)
					.setTraditional()
					.build();
			}
		}
	}

	private static void prepStates(Problem problem, Runnable block) {

		// make the ecalc from all the conf spaces
		List<SimpleConfSpace> confSpaces = problem.mskstar.states.stream()
			.map(state -> state.confSpace)
			.collect(Collectors.toList());
		try (EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpaces, problem.ffparams)
			.setParallelism(Parallelism.makeCpu(8))
			.build()) {

			// refresh the conf ecalcs
			for (MSKStar.State state : problem.mskstar.states) {
				state.confEcalc = new ConfEnergyCalculator(state.confEcalc, ecalc);
			}

			block.run();
		}
	}


	public static void main(String[] args) {
		bruteForce("2RL0 PPI", make2RL0PPI(false), 0.01);
		bruteForce("2RL0 only one mutant", make2RL0OnlyOneMutant(), 0.01);
		bruteForce("2RL0 space without wild type", make2RL0SpaceWithoutWildType(), 0.01);
	}

	private static void bruteForce(String name, Problem problem, double epsilon) {

		prepStates(problem, () -> {

			log("\n\n%s\n", name);

			// explicitly enumerate all the sequences
			List<Sequence> sequences = new ArrayList<>();
			if (problem.mskstar.seqSpace.containsWildTypeSequence()) {
				sequences.add(problem.mskstar.seqSpace.makeWildTypeSequence());
			}
			sequences.addAll(problem.mskstar.seqSpace.getMutants(problem.mskstar.maxSimultaneousMutations));

			for (Sequence sequence : sequences) {

				// calculate all the pfuncs
				Map<MSKStar.State,MathTools.DoubleBounds> stateBounds = new HashMap<>();
				for (MSKStar.State state : problem.mskstar.states) {

					// calculate the pfunc
					RCs rcs = sequence.makeRCs(state.confSpace);
					GradientDescentPfunc pfunc = new GradientDescentPfunc(state.confEcalc);
					pfunc.init(state.confTreeFactory.apply(rcs), rcs.getNumConformations(), epsilon);
					pfunc.compute();
					stateBounds.put(state, pfunc.makeResult().values.calcFreeEnergyBounds());
				}


				// build the assert statement
				// e.g.: assertSequence(sequences, "phe asp glu thr phe lys ile thr", -20.9375054147, -20.9429662431);
				MathTools.DoubleBounds objectiveBounds = problem.mskstar.objective.calc(stateBounds);
				log("assertSequence(sequences, \"%s\", %12.6f, %12.6f);",
					sequence.toString(Sequence.Renderer.ResTypeMutations),
					objectiveBounds.lower,
					objectiveBounds.upper
				);
			}
		});
	}

	@Test
	public void test2RL0() {
		Problem problem = make2RL0PPI(false);
		prepStates(problem, () -> check2RL0PPI(problem.mskstar));
	}

	@Test
	public void test2RL0BoundedMemory() {
		Problem problem = make2RL0PPI(true);
		prepStates(problem, () -> check2RL0PPI(problem.mskstar));
	}

	@Test
	public void test2RL0OnlyOneMutant() {
		Problem problem = make2RL0OnlyOneMutant();
		prepStates(problem, () -> check2RL0OnlyOneMutant(problem.mskstar));
	}

	@Test
	public void test2RL0SpaceWithoutWildType() {
	Problem problem = make2RL0SpaceWithoutWildType();
		prepStates(problem, () -> check2RL0WithoutWildType(problem.mskstar));
	}

	private static void assertSequence(List<MSKStar.SequenceInfo> sequences, String seqStr, double objectiveLower, double objectiveUpper) {

		// find the sequence
		MSKStar.SequenceInfo info = sequences.stream()
			.filter(seq -> seq.sequence.toString(Sequence.Renderer.ResTypeMutations).equals(seqStr))
			.findAny()
			.orElseThrow(() -> new NoSuchElementException("sequence " + seqStr + " was not found"));

		// make sure the objective bounds contain the expected bounds
		final double boundsEpsilon = 1e-6; // expected bounds only specified to 6 decimals
		assertThat(info.objective.lower, lessThanOrEqualTo(objectiveLower + boundsEpsilon));
		assertThat(info.objective.upper, greaterThanOrEqualTo(objectiveUpper - boundsEpsilon));
	}

	private static void assertSequenceOrder(List<MSKStar.SequenceInfo> sequences) {

		// make sure sequences come out in order of weakly increasing lower bounds
		for (int i=1; i<sequences.size(); i++) {
			assertThat(sequences.get(i-1).objective.lower, lessThanOrEqualTo(sequences.get(i).objective.lower));
		}
	}
}
