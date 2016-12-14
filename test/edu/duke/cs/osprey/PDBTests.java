/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import java.io.IOException;

/**
 *
 * @author mhall44
 */
public class PDBTests {
    
    public static void testPDBReadWrite(){
        
        EnvironmentVars.assignTemplatesToStruct = false;//this would create changes that would go to copy2
        Molecule m = PDBFileReader.readPDBFile("1CC8.copy.pdb", null);
        EnvironmentVars.assignTemplatesToStruct = true;
        PDBFileWriter.writePDBFile(m, "testResults/1CC8.copy2.pdb");
                
        try {
            int diffState = Runtime.getRuntime().exec("diff 1CC8.copy.pdb testResults/1CC8.copy2.pdb").getInputStream().read();
            if(diffState==-1)//indicates no difference
                System.out.println("PDB READ/WRITE TEST PASSED");
            else
                System.out.println("PDB READ/WRITE TEST FAILED");
        }
        catch(IOException e){System.err.println("IO EXCEPTION DURING DIFF");}


    }
}
