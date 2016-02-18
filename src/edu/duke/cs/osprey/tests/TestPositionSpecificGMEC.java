package edu.duke.cs.osprey.tests;

import org.junit.Test;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.control.GMECFinder;
import edu.duke.cs.osprey.restypes.PositionSpecificRotamerLibrary;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import junit.framework.TestCase;

public class TestPositionSpecificGMEC extends TestCase {

    ConfigFileParser cfp;
    String PDBFileLocation = "test/4NPD/4NPD.pdb";
    protected void setUp () throws Exception {
        super.setUp();
        
        EnvironmentVars.assignTemplatesToStruct = true;
        EnvironmentVars.resTemplates = null;
        
        String[] testArgs = new String[]{"-c", "test/4NPD/KStar.cfg", "Dummy command", "test/4NPD/DEE.cfg", "test/4NPD/System.cfg"};
        cfp = new ConfigFileParser(testArgs);//args 1, 3+ are configuration files
        cfp.loadData();
    }

    protected void tearDown () throws Exception {
        super.tearDown();
    }
    
    @Test
    public void testLoadPDBRotamers() throws Exception
    {
        ResidueTemplateLibrary library = PositionSpecificRotamerLibrary.generateLibraryFromPDB(PDBFileLocation);
    }
    
    @Test
    public void testFindAlternateGMEC() throws Exception
    {
        ResidueTemplateLibrary library = PositionSpecificRotamerLibrary.generateLibraryFromPDB(PDBFileLocation);
        GMECFinder gf = new GMECFinder(cfp);
        gf.calcGMEC();
    }


}
