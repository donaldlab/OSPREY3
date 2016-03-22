package edu.duke.cs.osprey.restypes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBRotamerReader;
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

	private Map<Integer, Map<String, List<ResidueTemplate>>> positionSpecificRotamers = new HashMap<>();
	public static PositionSpecificRotamerLibrary generateLibraryFromPDB(String pdbFileName)
	{

		PositionSpecificRotamerLibrary library = new PositionSpecificRotamerLibrary();
		Molecule m = PDBFileReader.readPDBFile(pdbFileName);
		PDBRotamerReader.createTemplates(m, library, pdbFileName);


		for(Residue r: m.residues)
		{
			outputResidueTemplateInfo(r);
		}
		return library;
	}

	// This constructor left intentionally empty. Never create one in a non-static fashion.
	private void PositionSpecificRotamerLibrary(){}

	private static void outputResidueTemplateInfo(Residue r)
	{
	}


	@Override
	public int numRotForResType(int pos, String resType, double phi, double psi) {
		Map<String, List<ResidueTemplate>> templatesAtDesignIndex = getTemplatesForDesignIndex(pos);
		return templatesAtDesignIndex.get(resType).size();
	}

	@Override
	public double getDihedralForRotamer(int pos, String resType, double phi, double psi,
			int rotNum, int dihedralNum) {
		Map<String, List<ResidueTemplate>> templatesAtDesignIndex = getTemplatesForDesignIndex(pos);
		return templatesAtDesignIndex.get(resType).get(pos).getRotamericDihedrals(phi, psi, rotNum, dihedralNum);
	}

	private Map<String, List<ResidueTemplate>> getTemplatesForDesignIndex (int pos) {
		if(!positionSpecificRotamers.containsKey(pos))
			positionSpecificRotamers.put(pos, new HashMap<>());
		return positionSpecificRotamers.get(pos);
	}

	public void addResidueTemplate (int residueIndex, String resType, ResidueTemplate allowedConformations) {
		// TODO Auto-generated method stub
		Map<String, List<ResidueTemplate>> templatesAtDesignIndex = getTemplatesForDesignIndex(residueIndex);
		if(!templatesAtDesignIndex.containsKey(resType))
			templatesAtDesignIndex.put(resType, new ArrayList<ResidueTemplate>());
		templatesAtDesignIndex.get(resType).add(allowedConformations);
	}
}
