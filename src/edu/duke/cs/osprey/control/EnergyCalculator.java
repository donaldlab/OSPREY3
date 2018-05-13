/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.control;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.PositionConfSpace;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MoleculeModifierAndScorer;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Residue;

public class EnergyCalculator {

	public void run(ConfigFileParser cfp) {
		
		SearchProblem search = cfp.getSearchProblem();
		
		System.out.println();
		
		// TODO: refactor this so we don't use a search space
		// then we don't have to require all the params below
		
		// this calculator doesn't really work with perturbations
		// the whole point is to evaluate the input structure without changing it
		if (cfp.params.getBool("doPerturbations")) {
			throw new Error("DEEPer not supported by Energy Calculator");
		}
		
		// also, it doesn't make much sense to run on non-wild-type rotamers
		if (!cfp.params.getBool("addWTRots")) {
			throw new Error("Must turn on addWTRots to use Energy Calculator. Otherwise, you're just evaluating arbitrary conformations");
		}
		
		System.out.println("Setting wild-type rotamers...");
		for (int i=0; i<search.confSpace.numPos; i++) {
			PositionConfSpace pos = search.confSpace.posFlex.get(i);
			
			if (pos.wtRCs.size() > 1) {
				System.out.println(String.format("WARNING: %d wild-type rotamers found for residue at design index %d. Picking first one aribtrarily",
					pos.wtRCs.size(), pos.designIndex
				));
			}
			RC rc = pos.wtRCs.get(0);
			
			// AAO 2016: switch template assumes an amino acid. will throw an exception otherwise.
			Residue res = search.confSpace.m.getResByPDBResNumber( search.confSpace.flexibleRes.get(i) );
            if(HardCodedResidueInfo.hasAminoAcidBB(res) && !res.fullName.startsWith("FOL")) {
            	search.confSpace.mutDOFs.get(i).switchToTemplate(rc.template);
            }
		}
		
		double energy;
		if (cfp.params.getBool("doMinimize")) {
			
			// run CCD over the continuous degrees of freedom
			// (on existing structure, not any RCs)
			System.out.println("Building energy function...");
			MoleculeModifierAndScorer objFunc = new MoleculeModifierAndScorer(search.fullConfE, search.confSpace);
			System.out.println(String.format("Minimizing %d degrees of freedom...", objFunc.getNumDOFs()));
			DoubleMatrix1D optDOFVals = new CCDMinimizer(objFunc, false).minimize().dofValues;
			
			System.out.println("Calculating energy...");
			energy = objFunc.getValue(optDOFVals);
			
		} else {
			
			// just evaluate the energy function
			System.out.println("Calculating energy...");
			energy = search.fullConfE.getEnergy();
		}
		
		System.out.println(String.format("Energy: %f kcal/mol\n", energy));
	}
}
