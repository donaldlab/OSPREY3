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

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.gmec.GMECConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.GMECFinder;
import edu.duke.cs.osprey.partcr.pickers.ConfPicker;
import edu.duke.cs.osprey.partcr.pickers.WalkingConfPicker;
import edu.duke.cs.osprey.partcr.scorers.RCScorer;
import edu.duke.cs.osprey.partcr.scorers.VolumeRCScorer;
import edu.duke.cs.osprey.partcr.splitters.NAryRCSplitter;
import edu.duke.cs.osprey.partcr.splitters.RCSplitter;

public class PartCRConfPruner implements GMECFinder.ConfPruner {
	
	private SearchProblem search;
	private double Ew;
	
	public PartCRConfPruner(SearchProblem search, double Ew) {
		this.search = search;
		this.Ew = Ew;
	}

	@Override
	public void prune(List<ScoredConf> confs, GMECConfEnergyCalculator confEcalc) {
				
		if (confs.isEmpty()) {
			return;
		}
		
		// IMPORTANT: make sure PartCR's ecalc matches the confEcalc!
		
		// this test will probably fail for EPIC/LUTE runs, which means PartCR will not work correctly there
		// because PartCR redistributes the terms of the energy function and assumes they're minimization energies
		// to avoid this problem, we need to refactor our energy calculation framework to be more modular
		// to allow this kind of redistribution, but be agnostic about what kinds of energies are being calculated
		// basically, we need a way to calculate PartCR.calcPosMinimizedEnergy() for every kind of energy ConfEnergyCalculator uses
		// this can't be fixed at runtime, so bail hard
		System.out.println("\nChecking to see if it's ok to use PartCR with this config...");
		double expectedMinimizedEnergy = confEcalc.calcEnergy(confs.get(0)).getEnergy();
		double partcrMinimizedEnergy = search.minimizedEnergy(confs.get(0).getAssignments());
		double absErr = Math.abs(expectedMinimizedEnergy - partcrMinimizedEnergy);
		// NOTE: ideally, these energies should be bitwise-identical, since search.minimizedEnergy() should be being used in both cases,
		// but I'm starting to suspect our minimizer is slightly non-deterministic, and/or the initial conditions are changing over time
		if (absErr > 1e-6) {
			throw new Error("This PartCR implmentation is not compatible with the current Osprey config!"
				+ "\nNeed to refactor energy calculation to make it work.\n"
				+ "\nOr the minimizer is just non-deterministic. The abs error is: " + absErr
			);
		}
		
		// TODO: also check energy decomposition
		
		System.out.println("Yup, we're good. Full speed ahead, captain!");
		
		// TODO: some enterprising student could try to optimize the PartCR configuration
		// i.e., see what settings/heuristics work best on a wide variety of design problems
		
		//ConfPicker picker = new FirstConfPicker();
		ConfPicker picker = new WalkingConfPicker();
		
		//RCScorer scorer = new NopRCScorer();
		//RCScorer scorer = new SplitsRCScorer();
		RCScorer scorer = new VolumeRCScorer();
		
		RCSplitter splitter = new NAryRCSplitter();
		//RCSplitter splitter = new BinaryRCSplitter();
		
		PartCR pcr = new PartCR(search, EnvironmentVars.curEFcnGenerator.ffParams, Ew, confs);
		pcr.setPicker(picker);
		pcr.setScorer(scorer);
		pcr.setSplitter(splitter);
		
		// run it!!
		int numStrikes = 3;
		pcr.autoIterate(numStrikes);
	}
}
