package edu.duke.cs.osprey;

import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.control.GMECFinder;
import edu.duke.cs.osprey.restypes.PositionSpecificRotamerLibrary;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;

public class TestGenerateRotamerLibrary {
    
    public static void main(String[] args)
    {
        // Init Environment for test
        initEnvironmentVariables();
        
        String PDBFileLocation = "test/4NPD/4NPD.pdb";
        
        ResidueTemplateLibrary library = PositionSpecificRotamerLibrary.generateLibraryFromPDB(PDBFileLocation);
    }
    
    private static void initEnvironmentVariables()
    {
        EnvironmentVars.assignTemplatesToStruct = true;
        EnvironmentVars.resTemplates = null;
        
        ConfigFileParser cfp = ConfigFileParser.makeFromResources(
            "/examples/4NPD/KStar.cfg",
            "/examples/4NPD/DEE.cfg",
            "/examples/4NPD/System.cfg"
        );
        cfp.loadData();
        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
    }

}
