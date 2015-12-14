package edu.duke.cs.osprey.restypes;

import java.util.ArrayList;

public abstract class ResidueTemplateLibrary {

	public ArrayList<ResidueTemplate> templates = new ArrayList<>();

	public ResidueTemplateLibrary() {
	}

	public int numDihedralsForResType(String resType) {
		//number of free dihedrals (sidechain dihedrals for actual amino acids)
		return firstTemplate(resType).numDihedrals;
	}

	/**
	 * PGC 2015: 
	 * Returns the number of rotamers for the specified residue type, for backbone dependent or backbone independent rotamers.
	 * @param pos 
	 * @param resType in three letter amino acid type
	 * @param phi The backbone phi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @param psi The backbone psi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @return The number of rotamers.  
	 */
	public abstract int numRotForResType(int pos, String resType, double phi, double psi);

	/**
	 * PGC 2015:
	 * get ideal dihedral value for a particular rotamer of a particular residue type, for backbone dependent or backbone independent rotamers.
	 * 
	 * @param resType in three letter amino acid type
	 * @param phi The backbone phi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @param psi The backbone psi angle for backbone dependent rotamers; will be ignored if backbone dependent rotamer libraries are not used.
	 * @param rotNum The rotamer number within this residue type.
	 * @param dihedralNum The dihedral number within this rotamer.
	 * @return
	 */
	public double getDihedralForRotamer(int pos, String resType, double phi, double psi,
			int rotNum, int dihedralNum) {
		return firstTemplate(resType).getRotamericDihedrals(phi, psi, rotNum, dihedralNum);
	}
	public double getDihedralForRotamer(int pos, String resType, int rotNum, int dihedralNum) {
		return firstTemplate(resType).getRotamericDihedrals(0, 0, rotNum, dihedralNum);
	}
	public double getDihedralForRotamer(String resType, int rotNum, int dihedralNum) {
		return firstTemplate(resType).getRotamericDihedrals(0, 0, rotNum, dihedralNum);
	}

	protected ResidueTemplate firstTemplate(String resType){
		//get the first template with the given residue type 
		//templates with the same type will all have the same rotamers, etc.
		//so this is good for looking up things like that
		for(ResidueTemplate templ : templates){
			if( templ.name.equalsIgnoreCase(resType) )
				return templ;
		}
		throw new RuntimeException("ERROR: template not found for residue type "+resType);
	}

	protected void addResidueTemplate(ResidueTemplate template)
	{
		templates.add(template);
	}
}