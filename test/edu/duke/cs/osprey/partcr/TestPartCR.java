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

package edu.duke.cs.osprey.partcr;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.ConfSearch;
import edu.duke.cs.osprey.confspace.ConfSearch.EnergiedConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.gmec.GMECFinder;
import edu.duke.cs.osprey.gmec.MinimizingConfEnergyCalculator;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestPartCR extends TestBase {
	
	@BeforeClass
	public static void before() {
		initDefaultEnvironment();
		MultiTermEnergyFunction.setNumThreads(1);
	}
	
	@Test
	public void testPartCr()
	throws Exception {
		
		SearchProblem search = makeSearch();
		
		// calc the energy matrix once
		search.emat = new SimpleEnergyMatrixCalculator.Cpu(2, EnvironmentVars.curEFcnGenerator.ffParams, search.confSpace, search.shellResidues).calcEnergyMatrix();
		
		List<EnergiedConf> expectedConfs = getConfs(search, false);
		List<EnergiedConf> partcrConfs = getConfs(search, true);
		
		// first make sure the expected confs are actually what we expect
		EnergiedConf minGMEC = expectedConfs.get(0);
		assertThat(minGMEC.getAssignments(), is(new int[] { 23, 23, 3, 5 }));
		assertThat(minGMEC.getScore(), isRelatively(-66.883873, 1e-6));
		assertThat(minGMEC.getEnergy(), isRelatively(-64.270196, 1e-6));
		assertThat(expectedConfs.size(), is(15));
		
		// then make sure PartCR did the right thing
		// the two conf lists should be identical, PartCR is only an optimization
		// except for maybe the minimized energies
		// the minimizer is slightly non-deterministic, or the initial conditions change over time
		final double scoreEpsilon = 1e-10;
		final double energyEpsilon = 1e-6;
		assertThat(expectedConfs.size(), is(partcrConfs.size()));
		for (int i=0; i<expectedConfs.size(); i++) {
			
			EnergiedConf expectedConf = expectedConfs.get(i);
			EnergiedConf partcrConf = partcrConfs.get(i);
			
			assertThat(expectedConf.getAssignments(), is(partcrConf.getAssignments()));
			assertThat(expectedConf.getScore(), isRelatively(partcrConf.getScore(), scoreEpsilon));
			assertThat(expectedConf.getEnergy(), isRelatively(partcrConf.getEnergy(), energyEpsilon));
		}
	}
	
	private SearchProblem makeSearch() {
		
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addMutable("39 43", "ALA VAL LEU ILE");
		resFlex.addFlexible("40 41");
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
		return new SearchProblem(
			"test", "examples/1CC8/1CC8.ss.pdb", 
			resFlex.flexResList, resFlex.allowedAAs, addWt, doMinimize, useEpic, new EPICSettings(), useTupleExpansion, new LUTESettings(),
			new DEEPerSettings(), moveableStrands, freeBBZones, useEllipses, useERef, addResEntropy, addWtRots, null,
			false, new ArrayList<>()
		);
	}
	
	private List<EnergiedConf> getConfs(SearchProblem search, boolean usePartCR)
	throws Exception {
		
		// configure DEE
		// don't worry about the pruning interval now, GMECFinder will configure it later
		double pruningInterval = 0;
		search.pruneMat = new PruningMatrix(search.confSpace, pruningInterval);
		boolean typeDep = false;
		double boundsThresh = 100;
		int algOption = 1;
		boolean useFlags = true;
		boolean useTriples = false;
		boolean preDACS = false;
		boolean useTupExp = false;
		double stericThresh = 100;
		PruningControl pruningControl = new PruningControl(
			search, pruningInterval, typeDep, boundsThresh, algOption,
			useFlags, useTriples, preDACS, search.useEPIC, useTupExp, stericThresh
		);
		pruningControl.setReportMode(PruningControl.ReportMode.Short);
		
		// configure A* search
		ConfSearchFactory astarFactory = new ConfSearchFactory() {
			@Override
			public ConfSearch make(EnergyMatrix emat, PruningMatrix pmat) {
				return new ConfAStarTree.Builder(search.emat, search.pruneMat)
					.setMPLP(new ConfAStarTree.MPLPBuilder()
						.setNumIterations(1)
					).setShowProgress(true)
					.build();
			}
		};
		
		// configure what energies to use
		ForcefieldParams ffparams = EnvironmentVars.curEFcnGenerator.ffParams;
		GMECConfEnergyCalculator.Async ecalc = MinimizingConfEnergyCalculator.make(ffparams, search);
		
		// configure the GMEC finder
		// NOTE: PartCR doesn't help as much with energy window designs
		// but we still want to test that it works correctly
		double I0 = 10;
		double Ew = 1;
		boolean useIMinDEE = true;
		boolean useContFlex = true;
		boolean useEPIC = false;
		boolean checkApproxE = true;
		boolean outputGMECStruct = false;
		boolean eFullConfOnly = false;
		File tmpFile = File.createTempFile("partcrTestConfs", ".txt");
		tmpFile.deleteOnExit();
		GMECFinder gmecFinder = new GMECFinder();
		gmecFinder.init(
			search, pruningControl, astarFactory, ecalc,
			Ew, useIMinDEE, I0, useContFlex, useTupExp, useEPIC,
			checkApproxE, outputGMECStruct, eFullConfOnly, tmpFile.getAbsolutePath(), stericThresh
		);
		gmecFinder.setLogConfsToConsole(false);
		
		// configure PartCR if needed
		if (usePartCR) {
			gmecFinder.setConfPruner(new PartCRConfPruner(search, Ew));
		}
		
		// GO GO GO!!
		return gmecFinder.calcGMEC();
	}
}
