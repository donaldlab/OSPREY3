package edu.duke.cs.osprey.energy;

import java.util.List;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldInteractions;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class ForcefieldInteractionsGenerator {
		
	public ForcefieldInteractions makeSingleRes(Residue res) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResidue(res);
		return interactions;
	}

	public ForcefieldInteractions makeResPair(Residue res1, Residue res2) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResiduePair(res1, res2);
		return interactions;
	}
	
	public ForcefieldInteractions makeIntraAndShell(Residue res, List<Residue> shellResidues) {
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		interactions.addResidue(res);
		for (Residue shellRes : shellResidues) {
			interactions.addResiduePair(res, shellRes);
		}
		return interactions;
	}
	
	public ForcefieldInteractions makeIntraAndDistributedShell(Residue res, List<Residue> shellResidues, int numPos, double singleWeight) {
		checkSingleWeight(singleWeight);
		return makeIntraAndShell(res, shellResidues);
	}

	public ForcefieldInteractions makeResPairAndDistributedShell(Residue res1, Residue res2, List<Residue> shellResidues, int numPos, double singleWeight) {
		checkSingleWeight(singleWeight);
		return makeResPair(res1, res2);
	}
	
	public ForcefieldInteractions makeFullConf(ConfSpace confSpace, List<Residue> shellResidues) {
		return makeFullConf(confSpace, shellResidues, null);
	}
	
	public ForcefieldInteractions makeFullConf(ConfSpace confSpace, List<Residue> shellResidues, Molecule mol) {
		
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
	
	public ForcefieldInteractions makeFullMol(Molecule mol) {
		
		ForcefieldInteractions interactions = new ForcefieldInteractions();
		
		for (Residue res1 : mol.residues) {
			
			// intra energy
			interactions.addResidue(res1);
			
			// pair energies
			for (Residue res2 : mol.residues) {
				interactions.addResiduePair(res1, res2);
			}
		}
		
		return interactions;
	}
	
	private Residue matchResidue(Residue res, Molecule mol) {
		if (mol == null) {
			return res;
		}
		return mol.residues.get(res.indexInMolecule);
	}
	
	private void checkSingleWeight(double singleWeight) {
		
		// only single weight of 1 supported
		if (singleWeight != 1) {
			throw new UnsupportedOperationException("only all-on-pairs shell energy distribution supported in GPU forcefields so far");
		}
	}
}
