package edu.duke.cs.osprey;

import org.junit.Test;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.control.GMECFinder;
import edu.duke.cs.osprey.restypes.PositionSpecificRotamerLibrary;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import junit.framework.TestCase;

public class TestPositionSpecificGMEC extends TestCase {

    ConfigFileParser cfp;
    String PDBFileLocation = "examples/4NPD/4NPD.pdb";
    protected void setUp () throws Exception {
        super.setUp();
        
        EnvironmentVars.assignTemplatesToStruct = true;
        EnvironmentVars.resTemplates = null;
        
        cfp = ConfigFileParser.makeFromFilePaths(
            "examples/4NPD/KStar.cfg",
            "examples/4NPD/DEE.cfg",
            "examples/4NPD/System.cfg"
        );
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
        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
        gf.calcGMEC();
    }


}
