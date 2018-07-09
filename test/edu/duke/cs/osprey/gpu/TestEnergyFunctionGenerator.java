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

package edu.duke.cs.osprey.gpu;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.TestBase;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.GpuEnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.energy.forcefield.GpuForcefieldEnergy;
import edu.duke.cs.osprey.gpu.opencl.GpuQueuePool;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tupexp.LUTESettings;

public class TestEnergyFunctionGenerator extends TestBase {
	
	@BeforeClass
	public static void before() {
		
		initDefaultEnvironment();
		
		// don't use energy function-level parallelism
		MultiTermEnergyFunction.setNumThreads(1);
	}
	
	private SearchProblem makeSearch(boolean doMinimize) {
	
		ResidueFlexibility resFlex = new ResidueFlexibility();
		resFlex.addMutable("39 43", "ALA");
		resFlex.addFlexible("40 41");
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
	
	private GpuEnergyFunctionGenerator makeGpuEgen() {
		return new GpuEnergyFunctionGenerator(makeDefaultFFParams(), new GpuQueuePool(1, 1));
	}
	
	private void testSingles(boolean doMinimize) {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		GpuEnergyFunctionGenerator gpuegen = makeGpuEgen();
		
		SearchProblem search = makeSearch(doMinimize);
		for (int pos=0; pos<search.confSpace.numPos; pos++) {
			Residue res = search.confSpace.posFlex.get(pos).res;
			
			double energy = egen.singleResEnergy(res).getEnergy();
			double gpuenergy = getGpuEnergy(gpuegen.singleResEnergy(res));
			
			assertThat(gpuenergy, isRelatively(energy));
		}
	}
	
	private double getGpuEnergy(GpuForcefieldEnergy efunc) {
		double energy = efunc.getEnergy();
		efunc.clean();
		return energy;
	}

	@Test
	public void testSinglesRigid() {
		testSingles(false);
	}
	
	@Test
	public void testSinglesContinuous() {
		testSingles(true);
	}
	
	private void testPairs(boolean doMinimize) {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		GpuEnergyFunctionGenerator gpuegen = makeGpuEgen();
		
		SearchProblem search = makeSearch(doMinimize);
		for (int pos1=0; pos1<search.confSpace.numPos; pos1++) {
			Residue res1 = search.confSpace.posFlex.get(pos1).res;
			
			for (int pos2=0; pos2<pos1; pos2++) {
				Residue res2 = search.confSpace.posFlex.get(pos2).res;
				
				double energy = egen.resPairEnergy(res1, res2).getEnergy();
				double gpuenergy = getGpuEnergy(gpuegen.resPairEnergy(res1, res2));
				
				assertThat(gpuenergy, isRelatively(energy));
			}
		}
	}
	
	@Test
	public void testPairsRigid() {
		testPairs(false);
	}
	
	@Test
	public void testPairsContinuous() {
		testPairs(true);
	}
	
	private void testIntraAndShell(boolean doMinimize) {
		
		EnergyFunctionGenerator egen = EnvironmentVars.curEFcnGenerator;
		GpuEnergyFunctionGenerator gpuegen = makeGpuEgen();
		
		SearchProblem search = makeSearch(doMinimize);
		for (int pos=0; pos<search.confSpace.numPos; pos++) {
			Residue res = search.confSpace.posFlex.get(pos).res;
			
			double energy = egen.intraAndShellEnergy(res, search.shellResidues).getEnergy();
			double gpuenergy = getGpuEnergy(gpuegen.intraAndShellEnergy(res, search.shellResidues));
			
			assertThat(gpuenergy, isRelatively(energy));
		}
	}
	
	@Test
	public void testIntraAndShellRigid() {
		testIntraAndShell(false);
	}
	
	@Test
	public void testIntraAndShellContinuous() {
		testIntraAndShell(true);
	}
}
