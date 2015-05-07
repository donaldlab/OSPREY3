/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

//testing that DOFs are applied correctly
//going to use the file 1CC8.ss.pdb so need to have it on hand

import edu.duke.cs.osprey.dof.FreeDihedral;
import edu.duke.cs.osprey.dof.ResidueTypeDOF;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.Protractor;

/**
 *
 * @author mhall44
 */
public class DOFTests {
    
    public static void testMutation(){
        Molecule m = PDBFileReader.readPDBFile("1CC8.ss.pdb");
        ResidueTypeDOF mutDOF = new ResidueTypeDOF(m.residues.get(37));//Ser 39 originally
        mutDOF.mutateTo("ALA");
        PDBFileWriter.writePDBFile(m, "1CC8.S39A.pdb");
        System.out.println("Wrote mutation test output: 1CC8.S39A.pdb");
    }
    
    public static void testDihedral(){
        
        Molecule m = PDBFileReader.readPDBFile("1CC8.ss.pdb");
        Residue res = m.residues.get(37);
        
        FreeDihedral chi1 = new FreeDihedral(res,0);//Ser 39
        FreeDihedral chi2 = new FreeDihedral(res,1);
        
        chi1.apply(45.);
        chi2.apply(-121);
        
        //measure dihedrals.  Start by collecting coordinates
        double N[] = res.getCoordsByAtomName("N");
        double CA[] = res.getCoordsByAtomName("CA");
        double CB[] = res.getCoordsByAtomName("CB");
        double OG[] = res.getCoordsByAtomName("OG");
        double HG[] = res.getCoordsByAtomName("HG");
        
        double chi1Measured = Protractor.measureDihedral(new double[][] {N,CA,CB,OG});
        double chi2Measured = Protractor.measureDihedral(new double[][] {CA,CB,OG,HG});
        
        System.out.println("chi1 applied as 45, measured as "+chi1Measured);
        System.out.println("chi2 applied as -121, measured as "+chi2Measured);
        System.out.println("Outputting dihedral-adjusted structure as 1CC8.dih.pdb");
        PDBFileWriter.writePDBFile(m, "1CC8.dih.pdb");
    }
    
}
