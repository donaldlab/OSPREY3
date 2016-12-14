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
        
        String[] testArgs = new String[]{"-c", "test/4NPD/KStar.cfg", "Dummy command", "test/4NPD/DEE.cfg", "test/4NPD/System.cfg"};
        ConfigFileParser cfp = new ConfigFileParser(testArgs);//args 1, 3+ are configuration files
        cfp.loadData();
        GMECFinder gf = new GMECFinder();
        gf.init(cfp);
    }

}
