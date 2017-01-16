/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;


import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 *
 * @author mhall44
 */
public class TestPDBReadWrite extends TestBase {
    
    
    @BeforeClass
    public static void before() {
            initDefaultEnvironment();
    }
    
    
    @Test
    public void testPDBReadWrite(){
        
        String origFileName = "examples/1CC8/1CC8.copy.pdb";
        String copiedFileName = "examples/1CC8/testResults/1CC8.copy2.pdb";
        
        EnvironmentVars.assignTemplatesToStruct = false;//this would create changes that would go to copy2
        Molecule m = PDBFileReader.readPDBFile(origFileName, null);
        EnvironmentVars.assignTemplatesToStruct = true;//return to normal in case we do other stuff later
        PDBFileWriter.writePDBFile(m, copiedFileName);
        
        try {
            byte[] origBytes = Files.readAllBytes(Paths.get(origFileName));
            byte[] copiedBytes = Files.readAllBytes(Paths.get(copiedFileName));
            
            String origString = new String(origBytes, StandardCharsets.UTF_8);
            String copiedString = new String(copiedBytes, StandardCharsets.UTF_8);
            
            assertThat(origString, is(copiedString));
        }
        catch(IOException e){throw new RuntimeException("IO EXCEPTION COMPARING FILES");}
    }
    
    
}
