/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.SearchSpace;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.pruning.PruningControl;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

/**
 *
 * @author mhall44
 */
public class ConfigFileParser {
    //An object that parses configuration files and uses them to initialize objects needed
    //in various calculations (i.e. for data file loading, conf space definition, pruning, etc.)
    
    ParamSet params = new ParamSet();
    
    public ConfigFileParser(String[] args){
        //parse all config files into params
        
        //check format of args
        if(!args[0].equalsIgnoreCase("-c"))
            throw new RuntimeException("ERROR: bad arguments (should start with -c)");
        
        params.addParamsFromFile(args[1]);//KStar.cfg file
        for(int argNum=3; argNum<args.length; argNum++)//System.cfg, etc.
            params.addParamsFromFile(args[argNum]);
    }
    
    
    //creation of objects needed in calculations like GMEC and K*
    SearchSpace getSearchSpace(){//this version is for a single search space...can modify for
        //additional states (unbound, etc.)
        
        String name = params.getValue("RUNNAME");
        
        ArrayList<String> flexRes = getFlexRes();
        ArrayList<ArrayList<String>> allowedAAs = getAllowedAAs();
        
        if(flexRes.size() != allowedAAs.size()){
            throw new RuntimeException("ERROR: Number of flexible positions different in flexible residue "
                    + "and allowed AA type parameters!");
        }
        
        return new SearchSpace( name, params.getValue("PDBNAME"), 
                flexRes, allowedAAs,
                params.getBool("AddWT",true), 
                params.getBool("doMinimize",false) );//CURRENTLY JUST SC MINIMIZATION...
    }
    
    
    
    ArrayList<String> getFlexRes(){
        //list of flexible residues.  PDB-based residue numbers
        //we'll include all flexible residues: for compatibility (MAY BE TEMPORARY),
        //all residues in a "StrandMut" record will be included here
        ArrayList<String> flexResList = new ArrayList<>();
        
        for(String param : params.searchParams("STRANDMUT")){//STRANDMUT0, etc.
            
            if(param.equalsIgnoreCase("STRANDMUTNUMS"))//means something different
                continue;
            
            String resListString = params.getValue(param);
            StringTokenizer tokenizer = new StringTokenizer(resListString);
            
            while(tokenizer.hasMoreTokens()){
                flexResList.add( tokenizer.nextToken() );
            }
        }
        
        return flexResList;
    }
    
    
    ArrayList<ArrayList<String>> getAllowedAAs(){
        //List allowed AA types for each flexible position
        //for compatibility (MAY BE TEMPORARY),
        //we'll go through all resAllowed records, ordering first by strand num
        //then pos num in strand
        
        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
        
        //handle better later (now assuming old-style numbering and no more than 10 strands)...
        for(int str=0; str<10; str++){
            ArrayList<String> resAllowedRecords = params.searchParams("RESALLOWED"+str);
            for(String param : resAllowedRecords){
                String allowedAAString = params.getValue(param);
                
                //parse AA types from allowedAAString
                ArrayList<String> resAllowedAAs = new ArrayList<>();
                StringTokenizer tokenizer = new StringTokenizer(allowedAAString);
            
                while(tokenizer.hasMoreTokens()){
                    resAllowedAAs.add( tokenizer.nextToken() );
                }
                
                allowedAAs.add(resAllowedAAs);
            }
        }
        
        return allowedAAs;
    }
    
    
    
    PruningControl setupPruning(SearchSpace searchSpace, double curEw){
        //setup pruning.  Conformations in searchSpace more than curEw over the GMEC are liable to pruning
        
        return new PruningControl(searchSpace, curEw, params.getBool("TYPEDEP",false), 
            params.getDouble("BOUNDSTHRESH",100), params.getInt("ALGOPTION",1), 
            params.getBool("USEFLAGS",true),
            params.getBool("USETRIPLES",false), false);//FOR NOW NO DCAS
    }
    
    
    
    
    
    //loading of data files
    //residue templates, rotamer libraries, forcefield parameters, and Ramachandran data
    void loadData(){
        
        EnvironmentVars.setDataDir(params.getValue("DataDir"));
        
        //a lot of this depends on forcefield type so figure that out first
        //general forcefield data loaded into the ForcefieldParams in EnvironmentVars
        ForcefieldParams curForcefieldParams = new ForcefieldParams(
                params.getValue("Forcefield","AMBER"),
                params.getBool("DISTDEPDIELECT",true),
		params.getDouble("DIELECTCONST",6.0),
                params.getDouble("VDWMULT",0.95),
		params.getBool("DOSOLVATIONE",true),
		params.getDouble("SOLVSCALE",0.5),
                params.getBool("HELECT",true),
                params.getBool("HVDW",true) );
        
        
        EnvironmentVars.curEFcnGenerator = new EnergyFunctionGenerator( 
                curForcefieldParams,
                params.getDouble("SHELLDISTCUTOFF",Double.POSITIVE_INFINITY) );
        
        String[] resTemplateFiles = getResidueTemplateFiles(curForcefieldParams.forcefld);
        
        ResidueTemplateLibrary resTemplates = new ResidueTemplateLibrary( resTemplateFiles, curForcefieldParams );
        
        //load template coordinates (necessary for all residues we might want to mutate to)
        //these will be matched to templates
        resTemplates.loadTemplateCoords("all_amino_coords.in");
        
        //load rotamer libraries; residues mentioned here will be matched to templates
        resTemplates.loadRotamerLibrary(params.getValue("ROTFILE","LovellRotamer.dat"));//see below; also gRotFile0 etc
        
        EnvironmentVars.resTemplates = resTemplates;
        
        
        //load RamachandranChecker here!
        
        
        /*
         * A lot of this is similar to this KSParser.setConfigPars code:
         * 
         hElect = (new Boolean((String)rParams.getValue("HELECT", "true"))).booleanValue();
		hVDW = (new Boolean((String)rParams.getValue("HVDW", "true"))).booleanValue();
		hSteric = (new Boolean((String)rParams.getValue("HSTERIC","false"))).booleanValue();
		distDepDielect = (new Boolean((String)rParams.getValue("DISTDEPDIELECT","true"))).booleanValue();
		dielectConst = (new Double((String)rParams.getValue("DIELECTCONST","6.0"))).doubleValue();
		doDihedE = (new Boolean((String)rParams.getValue("DODIHEDE","false"))).booleanValue();
		doSolvationE = (new Boolean((String)rParams.getValue("DOSOLVATIONE","true"))).booleanValue();
		solvScale = (new Double((String)rParams.getValue("SOLVSCALE","0.5"))).doubleValue();
		softvdwMultiplier = (new Double((String)rParams.getValue("VDWMULT","0.95"))).doubleValue();
		stericThresh = (new Double((String)rParams.getValue("STERICTHRESH","0.4"))).doubleValue();
		softStericThresh = (new Double((String)rParams.getValue("SOFTSTERICTHRESH","1.5"))).doubleValue();
		EnvironmentVars.setDataDir(rParams.getValue("DATADIR","./"));
		EnvironmentVars.setForcefld(rParams.getValue("FORCEFIELD","AMBER"));
		double entropyScale = (new Double((String)rParams.getValue("ENTROPYSCALE","0.0"))).doubleValue();
		EnvironmentVars.setEntropyScale(entropyScale);
		
		EnvironmentVars.setAArotLib(EnvironmentVars.getDataDir().concat(rParams.getValue("ROTFILE","LovellRotamer.dat")));

                EnvironmentVars.autoFix = new Boolean((String)rParams.getValue("AUTOFIX","true")).booleanValue();
                
                String ramaGlyFile = (String)rParams.getValue("RAMAGLYFILE","rama500-gly-sym.data");

                if( ! ramaGlyFile.equalsIgnoreCase("none") ){
                    String ramaFiles[] = { EnvironmentVars.dataDir + ramaGlyFile,
                    EnvironmentVars.dataDir + (String)rParams.getValue("RAMAPROFILE","rama500-pro.data"),
                    EnvironmentVars.dataDir + (String)rParams.getValue("RAMAGENFILE","rama500-general.data"),
                    EnvironmentVars.dataDir + (String)rParams.getValue("RAMAPREPROFILE","rama500-prepro.data")
                    };
                    RamachandranChecker.getInstance().readInputFiles( ramaFiles );
                }

         */
    }
    
    
    String[] getResidueTemplateFiles (ForcefieldParams.FORCEFIELD forcefld){
        //return names of residue template files
        
        //template file names are currently fixed
        String aaFilename=null, aaNTFilename=null, aaCTFilename=null, grFilename=null;
        
        switch(forcefld){
                case AMBER:
                        //KER: This is for the standard amber parameters
                        aaFilename =  "all_amino94.in";
                        aaNTFilename =  "all_aminont94.in";
                        aaCTFilename =  "all_aminoct94.in";
                        grFilename = "all_nuc94_and_gr.in";
                        break;
                case CHARMM22:
                        //KER: This if for using the charmm22 parameters:
                        aaFilename = "all_amino_charmm22.txt";
                        aaNTFilename = "all_amino_charmm22_nt.txt";
                        aaCTFilename = "all_amino_charmm22_ct.txt";
                        grFilename = "all_nuc_and_gr_charmm.in";
                        break;
                case CHARMM19NEUTRAL:
                        //KER: This is for CHARMM19 parameters:
                        aaFilename =  "all_amino_charmm19_neutral.in";
                        aaNTFilename =  "all_amino_charmm19_neutral_nt.in";
                        aaCTFilename =  "all_amino_charmm19_neutral_ct.in";
                        grFilename = "all_nuc_and_gr_charmm.in";
                        break;
                case CHARMM19:
                        aaFilename =  "all_amino_charmm19.in";
                        aaNTFilename =  "all_amino_charmm19_nt.in";
                        aaCTFilename =  "all_amino_charmm19_ct.in";
                        grFilename = "all_nuc_and_gr_charmm.in";
                        break;
                default:
                        System.out.println("FORCEFIELD not recognized...Exiting");
                        System.exit(0);
        }
        
        return new String[] {
            aaFilename, aaNTFilename, aaCTFilename, grFilename
        };
    }
    
    
    
}
