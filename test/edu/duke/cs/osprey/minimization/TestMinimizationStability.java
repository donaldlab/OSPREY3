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

package edu.duke.cs.osprey.minimization;

import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.ParameterizedMoleculeCopy;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.FFInterGen;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestMinimizationStability extends TestBase {
	
	private static SearchProblem search;
	private static List<ScoredConf> confs;
	private static ForcefieldParams ffparams;
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		
		// for these small problems, more than one thread is actually slower
		MultiTermEnergyFunction.setNumThreads(1);
		
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addMutable("39 43", "ALA");
		resFlex.addFlexible("40 41 42 44 45");
		boolean doMinimize = true;
		boolean addWt = true;
		boolean useEpic = false;
		boolean useTupleExpansion = false;
		boolean useEllipses = false;
		boolean useERef = false;
		boolean addResEntropy = false;
		boolean addWtRots = false;
		ArrayList<String[]> moveableStrands = new ArrayList<String[]>();
		ArrayList<String[]> freeBBZones = new ArrayList<String[]>();
		
		search = new SearchProblem(
			"test", "examples/1CC8/1CC8.ss.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
		
		ffparams = makeDefaultFFParams();
		
		// calc the energy matrix
		File ematFile = new File("/tmp/testMinimizationStability.emat.dat");
		if (ematFile.exists()) {
			search.emat = (EnergyMatrix)ObjectIO.readObject(ematFile.getAbsolutePath(), false);
		} else {
			search.emat = new SimpleEnergyMatrixCalculator.Cpu(2, ffparams, search.confSpace, search.shellResidues).calcEnergyMatrix();
			ObjectIO.writeObject(search.emat, ematFile.getAbsolutePath());
		}
		
		// don't prune anything
		search.pruneMat = new PruningMatrix(search.confSpace, 1000);
		
		// build A* tree
		ConfAStarTree tree = new ConfAStarTree.Builder(search.emat, search.pruneMat)
			.setMPLP(new ConfAStarTree.MPLPBuilder()
				.setNumIterations(1)
			).build();
		
		// get the confs
		final int NumConfs = 10;
		confs = new ArrayList<>();
		for (int i=0; i<NumConfs; i++) {
			confs.add(tree.nextConf());
		}
	}
	
	@Test
	public void sharedMol() {
		
		ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
		
		System.out.println("minimizing in ascending order...");
		double[] ascendingEnergies = new double[confs.size()];
		for (int i=0; i<confs.size(); i++) {
			ascendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		System.out.println("minimizing in descending order...");
		double[] descendingEnergies = new double[confs.size()];
		for (int i=confs.size() - 1; i>=0; i--) {
			descendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		// energies should be exactly the same
		checkResults(ascendingEnergies, descendingEnergies, 0);
	}
	
	@Test
	public void separateMols() {
		
		System.out.println("minimizing in ascending order...");
		double[] ascendingEnergies = new double[confs.size()];
		for (int i=0; i<confs.size(); i++) {
			ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
			ascendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		System.out.println("minimizing in descending order...");
		double[] descendingEnergies = new double[confs.size()];
		for (int i=confs.size() - 1; i>=0; i--) {
			ParameterizedMoleculeCopy pmol = new ParameterizedMoleculeCopy(search.confSpace);
			descendingEnergies[i] = minimize(confs.get(i), pmol);
		}
		
		// energies should be exactly the same
		checkResults(ascendingEnergies, descendingEnergies, 0);
	}
	
	private double minimize(ScoredConf conf, ParameterizedMoleculeCopy pmol) {
		CpuConfMinimizer minimizer = new CpuConfMinimizer.Builder(
			ffparams,
			(mol) -> FFInterGen.makeFullConf(search.confSpace, search.shellResidues, mol),
			search.confSpace
		).build();
		double energy = minimizer.minimize(conf).getEnergy();
		minimizer.cleanup();
		return energy;
	}
	
	private void checkResults(double[] ascendingEnergies, double[] descendingEnergies, double energyEpsilon) {
		
		double sum = 0;
		
		int n = ascendingEnergies.length;
		for (int i=0; i<n; i++) {
			
			double absErr = Math.abs(ascendingEnergies[i] - descendingEnergies[i]);
			sum += absErr*absErr;
			
			System.out.println(String.format("asc: %20.12f   desc: %20.12f   err: %e   %s",
				ascendingEnergies[i], descendingEnergies[i],
				absErr,
				absErr > energyEpsilon ? "   <-- too much err" : ""
			));
		}
		
		double rmsd = Math.sqrt(sum/n);
		System.out.println(String.format("RMS err: %12.12f", rmsd));
		
		// check the results
		for (int i=0; i<n; i++) {
			assertThat(ascendingEnergies[i], isAbsolutely(descendingEnergies[i], energyEpsilon));
		}
	}
	
}
