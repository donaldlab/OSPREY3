package edu.duke.cs.osprey.confspace;

import java.util.ArrayList;
import java.util.List;

import edu.duke.cs.osprey.restypes.DAminoAcidHandler;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.restypes.HardCodedResidueInfo;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;

public class Strand {
	
	public final Molecule mol;
	
	public Strand(Molecule mol) {
		// by default, use the full-length molecule
		this(mol, mol.residues.get(0).getPDBResNumber(), mol.residues.get(mol.residues.size() - 1).getPDBResNumber());
	}
	
	public Strand(Molecule mol, String firstResNumber, String lastResNumber) {
		
		// build our molecule copy from the residue subset
		this.mol = new Molecule();
		for (Residue res : mol.getResRangeByPDBResNumber(firstResNumber, lastResNumber)) {
			this.mol.appendResidue(new Residue(res));
		}
	}
	
	public List<String> assignTemplatesAndMarkBonds(GenericResidueTemplateLibrary templateLib, boolean deleteNonTemplateResidues) {
		
		List<String> deletedResidueNames = new ArrayList<>();
		
		for (int resNum = mol.residues.size() - 1; resNum >= 0; resNum--) {
			// go through residues backwards so we can delete some if needed
			
			Residue res = mol.residues.get(resNum);
			
			// We accept D-amino acid named using the usual L names,
			// but must change them here so the right template name is used
			DAminoAcidHandler.tryRenamingAsD(res);
			
			// try to assign the template
			boolean templateAssigned = res.assignTemplate(templateLib);
			
			// delete unrecognized residues if desired
			if (!templateAssigned && deleteNonTemplateResidues) {
				deletedResidueNames.add(res.fullName);
				mol.deleteResidue(resNum);
			}
		}
		
		// assigning templates marks intra-res bonds; we can now mark inter-res too
		HardCodedResidueInfo.markInterResBonds(mol);
		
		return deletedResidueNames;
	}
}
