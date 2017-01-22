package edu.duke.cs.osprey.restypes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.Residue;

/***
 * Rotamer library which specifies specific rotamers at 
 * each design position. Used for generating alternate
 * conformations when input structure is sufficiently high 
 * resolution.
 * @author JJ
 *
 */
public class PositionSpecificRotamerLibrary extends ResidueTemplateLibrary {

	private Map<Integer,Map<String,List<ResidueTemplate>>> positionSpecificRotamers = new HashMap<>();
	
	public PositionSpecificRotamerLibrary(Strand strand) {
		
		for (int i=0; i<strand.mol.residues.size(); i++) {
			
			// gather all the residues for each position, even alts
			List<Residue> residues = new ArrayList<>();
			residues.add(strand.mol.residues.get(i));
			residues.addAll(strand.mol.getAlternates(i));
			
			addResidueTemplate(i, residues.get(0).template.name, ResidueTemplate.makeFromResidueConfs(residues));
		}
	}

	@Override
	public int numRotForResType(int pos, String resType, double phi, double psi) {
		return getTemplatesForDesignIndex(pos).get(resType).size();
	}

	@Override
	public double getDihedralForRotamer(int pos, String resType, double phi, double psi, int rotNum, int dihedralNum) {
		return getTemplatesForDesignIndex(pos).get(resType).get(pos).getRotamericDihedrals(phi, psi, rotNum, dihedralNum);
	}

	public Map<String,List<ResidueTemplate>> getTemplatesForDesignIndex(int pos) {
		if(!positionSpecificRotamers.containsKey(pos)) {
			positionSpecificRotamers.put(pos, new HashMap<>());
		}
		return positionSpecificRotamers.get(pos);
	}

	private void addResidueTemplate(int residueIndex, String resType, ResidueTemplate allowedConformations) {
		Map<String,List<ResidueTemplate>> templatesAtDesignIndex = getTemplatesForDesignIndex(residueIndex);
		if(!templatesAtDesignIndex.containsKey(resType)) {
			templatesAtDesignIndex.put(resType, new ArrayList<>());
		}
		templatesAtDesignIndex.get(resType).add(allowedConformations);
	}
}
