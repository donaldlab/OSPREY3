/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;


@Deprecated
public class EnvironmentVars {
    // TODO: get rid of global state
	
	//key parameter sets to use throughout the program for energy-function and flexibility reference
        public static EnergyFunctionGenerator curEFcnGenerator;
        public static ResidueTemplateLibrary resTemplates;
        
        //Regulation of structure read-in/template assignment
        public static boolean assignTemplatesToStruct = true;//Assign templates when we read in a template.
        //We'll keep this true unless explicitly set otherwise since most calculations require templates
        public static boolean deleteNonTemplateResidues = true;//Delete residues for which we don't have a template
        
        public static boolean alwaysIdealizeSidechainsAfterMutation;
        
        public static boolean useMPI = false;//distribute things like energy matrix calculations, K* calculation for sequences using MPI
        
        public static double DUNBRACK_PROBABILTY_CUTOFF = 0.001;
        
        
        //data files directory
        private static String dataDir = "./";
        
        //Some people like to know if their residues are being deleted
        //so we write this in a special warning log
        //If you indicate "None" for the filename then this is null (no log will be made)
        public static SpecialWarningLog deletedResWarningLog = null; 
        
        
        public static String getDataDir() {
		return dataDir;
	}

	public static void setDataDir(String dd) {
		if(!dd.endsWith("/") || !dd.endsWith("\\")){
			dd = dd.concat("/");
		}
		EnvironmentVars.dataDir = dd;
	}
        
        
        static void openSpecialWarningLogs(ConfigFileParser cfp){
            //so far just the one deleted residues log
            String deletedResLogName = cfp.params.getValue("DeletedResWarningLog");
            if( ! deletedResLogName.equalsIgnoreCase("None") ){
                deletedResWarningLog = new SpecialWarningLog(deletedResLogName);
            }
        }
        
        static void closeSpecialWarningLogs(){
            if(deletedResWarningLog!=null)
                deletedResWarningLog.close();
        }
        
        //ALL THE BELOW ARE STILL IN QUESTION!!!
        
	/*
	        

	//static BigInteger maxKSconfs;
	//static boolean useMaxKSconfs = false;
	
	public static boolean useEntropy = false;
	private static double[] entropy = null;
	public static double entropyScaling = 0;
	
	public static String ksConfDir = "ksConfs";
	
	public static RotamerLibrary aaRotLib;

        public static boolean autoFix = true;//Should structures being read in be autofixed?
	
        
        
        public static boolean useMPLP = false;
        public static int MPLP_iterations = 100;

	
	public static void setForcefld(String frcefld) {
            //set up the current force field to be what the string indicates (AMBER or whatever)
            curForcefieldParams = new ForcefieldParams(frcefld);
	}
	

	public static double[] getEntropyTerms() {
		if(entropy == null){
			double[] entr = {0,0.5,0.75,0.75,0.58,0.99,0.97,1.14,1.53,1.19,1.12,2.21,2.13,0.99,0.99,0.99,0.61,1.65,0.81,2.02,0};
			entropy = entr;	
			for(int i=0;i<entropy.length;i++)
				entropy[i] *= entropyScaling;
		}
		
		return entropy;
	}
	
	public static void setEntropyScale(double es){
		if(es>0)
			useEntropy = true;
		entropyScaling = es;
	}
	
	public static void setAArotLib(String rl){
		aaRotLib = new RotamerLibrary(rl, true);
	}
        */

}

