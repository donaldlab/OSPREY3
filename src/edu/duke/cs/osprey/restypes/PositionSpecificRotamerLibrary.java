package edu.duke.cs.osprey.restypes;

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
    
    public static PositionSpecificRotamerLibrary generateLibraryFromPDB(String pdbFileName)
    {
        Molecule m = PDBFileReader.readPDBFile(pdbFileName);
        PDBRotamerReader.addAlternates(m, pdbFileName);
        
        PositionSpecificRotamerLibrary library = new PositionSpecificRotamerLibrary();
        
        for(Residue r: m.residues)
        {
            outputResidueTemplateInfo(r);
            library.addResidueTemplate(r.template);
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
		// TODO Auto-generated method stub
		return 0;
	}
}
