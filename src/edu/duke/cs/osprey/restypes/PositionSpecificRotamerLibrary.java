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
public class PositionSpecificRotamerLibrary {
    
    public static PositionSpecificRotamerLibrary generateLibraryFromPDB(String pdbFileName)
    {
        Molecule m = PDBFileReader.readPDBFile(pdbFileName);
        PDBRotamerReader.addAlternates(m, pdbFileName);
        
        for(Residue r: m.residues)
        {
            outputResidueTemplateInfo(r);
        }
        return new PositionSpecificRotamerLibrary();
    }
    
    private void PositionSpecificRotamerLibrary()
    {
        
    }
    
    private static void outputResidueTemplateInfo(Residue r)
    {
    }
}
