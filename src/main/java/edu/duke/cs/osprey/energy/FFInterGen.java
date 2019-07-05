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

package edu.duke.cs.osprey.energy;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

/**
 * forcefield interactions generator
 * 
 * Use ResInterGen and ForcefieldInteractions(ResidueInteractions,Molecule) instead
 */
@Deprecated
public class FFInterGen {
	
	public static ForcefieldInteractions makeSingleRes(ConfSpace confSpace, int pos, Molecule mol) {
		return makeSingleRes(matchResidue(confSpace.posFlex.get(pos).res, mol));
	}

	public static ForcefieldInteractions makeSingleRes(SimpleConfSpace confSpace, int pos, Molecule mol) {
		return makeSingleRes(mol.getResByPDBResNumber(confSpace.positions.get(pos).resNum));
	}
		
	public static ForcefieldInteractions makeSingleRes(Residue res) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResidue(res);
		return interactions;
	}
	
	public static ForcefieldInteractions makeResPair(ConfSpace confSpace, int pos1, int pos2, Molecule mol) {
		return makeResPair(
			matchResidue(confSpace.posFlex.get(pos1).res, mol),
			matchResidue(confSpace.posFlex.get(pos2).res, mol)
		);
	}
	
	public static ForcefieldInteractions makeResPair(SimpleConfSpace confSpace, int pos1, int pos2, Molecule mol) {
		return makeResPair(
			mol.getResByPDBResNumber(confSpace.positions.get(pos1).resNum),
			mol.getResByPDBResNumber(confSpace.positions.get(pos2).resNum)
		);
	}

	public static ForcefieldInteractions makeResPair(Residue res1, Residue res2) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResiduePair(res1, res2);
		return interactions;
	}
	
	public static ForcefieldInteractions makeIntraAndShell(ConfSpace confSpace, int pos, List<Residue> shellResidues, Molecule mol) {
		return makeIntraAndShell(
			matchResidue(confSpace.posFlex.get(pos).res, mol),
			matchResidues(shellResidues, mol)
		);
	}
	
	public static ForcefieldInteractions makeIntraAndShell(SimpleConfSpace confSpace, int pos, Molecule mol) {
		return makeIntraAndShell(
			mol.getResByPDBResNumber(confSpace.positions.get(pos).resNum),
			mol.getResiduesByPDBResNumbers(confSpace.shellResNumbers)
		);
	}
	
	public static ForcefieldInteractions makeIntraAndShell(Residue res, List<Residue> shellResidues) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResidue(res);
		for (Residue shellRes : shellResidues) {
			interactions.addResiduePair(res, shellRes);
		}
		return interactions;
	}
	
	public static ForcefieldInteractions makeFullConf(ConfSpace confSpace, List<Residue> shellResidues) {
		return makeFullConf(confSpace, shellResidues, null);
	}
	
	public static ForcefieldInteractions makeFullConf(ConfSpace confSpace, List<Residue> shellResidues, Molecule mol) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (int pos1=0; pos1<confSpace.posFlex.size(); pos1++) {
			Residue res1 = confSpace.posFlex.get(pos1).res;
			res1 = matchResidue(res1, mol);
			
			// intra energy
			interactions.addResidue(res1);
			
			// pair energies
			for (int pos2=0; pos2<pos1; pos2++) {
				Residue res2 = confSpace.posFlex.get(pos2).res;
				res2 = matchResidue(res2, mol);
				interactions.addResiduePair(res1, res2);
			}
			
			// shell energies
			for (Residue shellRes : shellResidues) {
				shellRes = matchResidue(shellRes, mol);
				interactions.addResiduePair(res1, shellRes);
			}
		}
		
		return interactions;
	}
	
	public static ForcefieldInteractions makeFullConf(SimpleConfSpace confSpace, Molecule mol) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (int pos1=0; pos1<confSpace.positions.size(); pos1++) {
			Residue res1 = mol.getResByPDBResNumber(confSpace.positions.get(pos1).resNum);
			
			// intra energy
			interactions.addResidue(res1);
			
			// pair energies
			for (int pos2=0; pos2<pos1; pos2++) {
				Residue res2 = mol.getResByPDBResNumber(confSpace.positions.get(pos2).resNum);
				interactions.addResiduePair(res1, res2);
			}
			
			// shell energies
			for (String resNum : confSpace.shellResNumbers) {
				Residue shellRes = mol.getResByPDBResNumber(resNum);
				interactions.addResiduePair(res1, shellRes);
			}
		}
		
		return interactions;
	}
	
	public static ForcefieldInteractions makeFullMol(Molecule mol) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (Residue res1 : mol.residues) {
			
			// intra energy
			interactions.addResidue(res1);
			
			// pair energies
			for (Residue res2 : mol.residues) {
				if (res1 != res2) {
					interactions.addResiduePair(res1, res2);
				}
			}
		}
		
		return interactions;
	}
	
	private static Residue matchResidue(Residue res, Molecule mol) {
		if (mol == null) {
			return res;
		}
		return mol.residues.get(res.indexInMolecule);
	}
	
	private static List<Residue> matchResidues(List<Residue> residues, Molecule mol) {
		if (mol == null) {
			return residues;
		}
		List<Residue> out = new ArrayList<>();
		for (Residue res : residues) {
			out.add(matchResidue(res, mol));
		}
		return out;
	}
}
