/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.astar.GMECMutSpace;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.RamachandranChecker;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.Forcefield;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams.SolvationForcefield;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.tools.StringParsing;
import edu.duke.cs.osprey.tupexp.LUTESettings;

/**
 *
 * @author mhall44
 */
public class ConfigFileParser {
    //An object that parses configuration files and uses them to initialize objects needed
    //in various calculations (i.e. for data file loading, conf space definition, pruning, etc.)
    
    public final ParamSet params;
    
    public static ConfigFileParser makeFromFilePaths(String ... paths) {
        return makeFromFilePaths(Arrays.asList(paths));
    }
    
    public static ConfigFileParser makeFromFilePaths(Iterable<String> paths) {
        ConfigFileParser cfp = new ConfigFileParser();
        for (String path : paths) {
            cfp.params.addParamsFromFile(new File(path));
        }
        return cfp;
    }
    
    public static ConfigFileParser makeFromFiles(Iterable<File> files) {
        ConfigFileParser cfp = new ConfigFileParser();
        for (File file : files) {
            cfp.params.addParamsFromFile(file);
        }
        return cfp;
    }
    
    public ConfigFileParser() {
        params = new ParamSet();
    }
    
    public ConfigFileParser(ConfigFileParser other) {
        params = new ParamSet(other.params);
    }
    
    protected DEEPerSettings setupDEEPer(){
        //Set up the DEEPerSettings object, including the PertSet (describes the perturbations)
        DEEPerSettings dset = new DEEPerSettings(
                params.getBool("doPerturbations"),
                params.getRunSpecificFileName("perturbationFile", ".pert"),
                params.getBool("selectPerturbations"),
                params.getValue("startingPerturbationFile"),
                params.getBool("onlyStartingPerturbations"),
                params.getDouble("maxShearParam"),
                params.getDouble("maxBackrubParam"),
                params.getBool("selectLCAs"),
                getFlexRes(), 
                params.getValue("PDBNAME"),
                params.getBool("DORAMACHECK"),
				EnvironmentVars.resTemplates
        );
        
        dset.loadPertFile(null);//load the PertSet from its file
        return dset;
    }
    
    
    private ArrayList<String[]> freeBBZoneTermini(){
        //Read the termini of the BBFreeBlocks, if any
        ArrayList<String[]> ans = new ArrayList<>();
        
        for(String rt : params.searchParams("BBFREEBLOCK")){
            //So for example BBFREEBLOCK0 120 125 would mean make a BBFreeBlock for res 120-125
            //lexical ordering for blocks is OK
            String strandLimitsString = params.getValue(rt);

            String[] termini = 
                { StringParsing.getToken(strandLimitsString, 1),
                    StringParsing.getToken(strandLimitsString, 2) };

            ans.add(termini);
        }
        
        return ans;
    }
    
    
    private ArrayList<String[]> moveableStrandTermini(){
        //Read the strands that are going to translate and rotate
        //Let's say they can do this regardless of what doMinimize says (that's for sidechains)
        ArrayList<String[]> ans = new ArrayList<>();
        
        for(String rt : params.searchParams("STRANDROTTRANS")){
            if(params.getBool(rt)){
                //So rt = STRANDROTTRANS0 here means strand 0 should translate & rotate
                //OK to go through these params in lexical ordering
                String strandNum = rt.substring(14);
                String strandLimitsString = params.getValue("STRAND"+strandNum);
                
                String[] termini = 
                    { StringParsing.getToken(strandLimitsString, 1),
                        StringParsing.getToken(strandLimitsString, 2) };
                
                ans.add(termini);
            }
        }
        
        return ans;
    }
    
    //creation of objects needed in calculations like GMEC and K*
    public SearchProblem getSearchProblem(){//this version is for a single search problem...can modify for
        //additional states (unbound, etc.)
        
        String name = params.getValue("RUNNAME");
        
        ArrayList<String> flexRes = getFlexRes();
        ArrayList<ArrayList<String>> allowedAAs = getAllowedAAs();
        
        if(flexRes.size() != allowedAAs.size()){
            throw new RuntimeException("ERROR: Number of flexible positions different in flexible residue "
                    + "and allowed AA type parameters!");
        }
        
        System.out.println("CREATING SEARCH PROBLEM.  NAME: "+name);
        
        ArrayList<String[]> moveableStrands = moveableStrandTermini();
        ArrayList<String[]> freeBBZones = freeBBZoneTermini();
        DEEPerSettings dset = setupDEEPer();
        
        SearchProblem search = new SearchProblem( name, params.getFile("PDBNAME").getAbsolutePath(), 
                flexRes, allowedAAs,
                params.getBool("AddWT"), 
                params.getBool("doMinimize"),
                params.getBool("UseEPIC"),
                new EPICSettings(params),
                params.getBool("UseTupExp"),
                new LUTESettings(params),
                dset, moveableStrands, freeBBZones,
                params.getBool("useEllipses"),
                params.getBool("useERef"),
                params.getBool("AddResEntropy"),
                params.getBool("addWTRots"),
                null,
                params.getBool("useVoxelG"),
                getWtRotOnlyRes()
        );
        
        search.numEmatThreads = params.getInt("EmatThreads");
        
        return search;
    }
    
    
    
    public ArrayList<String> getFlexRes(){
        //list of flexible residues.  PDB-based residue numbers
        //we'll include all flexible residues: for compatibility (MAY BE TEMPORARY),
        //all residues in a "StrandMut" record will be included here
        //"StrandMutNums" means something different and thus will be excluded
        ArrayList<String> flexResList = new ArrayList<>();
        
        int numStrands = params.searchParams("STRANDMUT").size() 
                - params.searchParams("STRANDMUTNUMS").size();
        
        //must go through the strands in order to get the right residues 
        for(int smNum=0; smNum<numStrands; smNum++){//must do these in order
            //so we get the right residues
            String param = "STRANDMUT"+smNum;//STRANDMUT0, etc
            
            String resListString = params.getValue(param);
            StringTokenizer tokenizer = new StringTokenizer(resListString);
            
            while(tokenizer.hasMoreTokens()){
                flexResList.add( tokenizer.nextToken() );
            }
        }
        
        return flexResList;
    }
    
    
    protected ArrayList<String> getWtRotOnlyRes(){
        //List of residues for which we'll only include the wild-type rotamer
        ArrayList<String> wtRotOnlyRes = new ArrayList<>();
        String val = params.getValue("WTRotOnlyRes");
        
        StringTokenizer tokenizer = new StringTokenizer(val);
        while(tokenizer.hasMoreTokens()){
            wtRotOnlyRes.add(tokenizer.nextToken());
        }
        
        return wtRotOnlyRes;
    }
    
    
    public ArrayList<ArrayList<String>> getAllowedAAs(){
        //List allowed AA types for each flexible position
        //We can accept either RESALLOWED0_0 (for flexible res 0 of strand 0)
        //or RESALLOWED255 (for residue with PDB number 255)
        //but which one is used should be consistent across all residues
        ArrayList<String> resAllowedRecords = params.searchParams("RESALLOWED");
        
        if(resAllowedRecords.isEmpty())//no flexible residues
            return new ArrayList<ArrayList<String>>();
        else {
            boolean usingStrandFormat = resAllowedRecords.get(0).contains("_");
            for(int rec=1; rec<resAllowedRecords.size(); rec++){
                if( usingStrandFormat != resAllowedRecords.get(rec).contains("_") ){
                    throw new RuntimeException("ERROR: Inconsistent formatting of resAllowed records"
                            + " (should be all by PDB residue number or all by strand)");
                }
            }
            
            if(usingStrandFormat)
                return getAllowedAAsByStrand();
            else
                return getAllowedAAsByPDBResNum();
        }
    }
        
    
    
    ArrayList<ArrayList<String>> getAllowedAAsByStrand(){
        //we'll go through all resAllowed records, ordering first by strand num
        //then pos num in strand
        
        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
        
        //handle better later (now assuming old-style numbering)...
        for(int str=0; true; str++){
            ArrayList<String> resAllowedRecords = params.searchParams("RESALLOWED"+str);
            int numRecordsInStrand = resAllowedRecords.size();
            
            //must go through residues in numerical order
            for(int recNum=0; recNum<numRecordsInStrand; recNum++){
                String param = "RESALLOWED" + str + "_" + recNum;
                    
                String allowedAAString = params.getValue(param);
                
                //parse AA types from allowedAAString
                ArrayList<String> resAllowedAAs = new ArrayList<>();
                StringTokenizer tokenizer = new StringTokenizer(allowedAAString);
            
                while(tokenizer.hasMoreTokens()){
                    resAllowedAAs.add( tokenizer.nextToken().toUpperCase() );
                }
                
                allowedAAs.add(resAllowedAAs);
            }
            
            if(numRecordsInStrand==0)//finished with strands that have flexible residues
                break;
        }
        
        return allowedAAs;
    }
    
    
    
    ArrayList<ArrayList<String>> getAllowedAAsByPDBResNum(){
        //we'll go through all resAllowed records, ordering first by strand num
        //then pos num in strand
        
        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
        ArrayList<String> flexRes = getFlexRes();
        
        for(int flexResNum=0; flexResNum<flexRes.size(); flexResNum++){
            String param = "RESALLOWED" + flexRes.get(flexResNum);
            String allowedAAString = params.getValue(param);

            //parse AA types from allowedAAString
            ArrayList<String> resAllowedAAs = new ArrayList<>();
            StringTokenizer tokenizer = new StringTokenizer(allowedAAString);

            while(tokenizer.hasMoreTokens()){
                resAllowedAAs.add( tokenizer.nextToken().toUpperCase() );
            }

            allowedAAs.add(resAllowedAAs);
        }
        
        return allowedAAs;
    }
    
    
    public PruningControl setupPruning(SearchProblem searchSpace, double pruningInterval, boolean useEPIC, boolean useTupExp){
        //setup pruning.  Conformations in searchSpace more than (Ew+Ival) over the GMEC are liable to pruning
        
        //initialize the pruning matrix for searchSpace, if not already initialized
        //or if pruningInterval lower (so it may have pruned tuples that shouldn't
        //be pruned with our new pruningInterval)
        boolean initPruneMat = false;
        if(searchSpace.pruneMat==null)
            initPruneMat = true;
        else if(searchSpace.pruneMat.getPruningInterval() < pruningInterval)
            initPruneMat = true;
        
        if(initPruneMat)
            searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, pruningInterval);
        
        return new PruningControl( searchSpace, pruningInterval, params.getBool("TYPEDEP"), 
            params.getDouble("BOUNDSTHRESH"), params.getInt("ALGOPTION"), 
            params.getBool("USEFLAGS"),
            params.getBool("USETRIPLES"), false, useEPIC, useTupExp,
            params.getDouble("STERICTHRESH") );//FOR NOW NO DACS
    }
    
    
    
    
    
    //loading of data files
    //residue templates, rotamer libraries, forcefield parameters, and Ramachandran data
    // TODO: update CFP members instead of EnvironmentVars, then replace loadData() with getters
    @Deprecated
    public void loadData(){
    	
    	Forcefield ff = Forcefield.get(params.getValue("Forcefield"));
        
        //a lot of this depends on forcefield type so figure that out first
        //general forcefield data loaded into the ForcefieldParams in EnvironmentVars
        ForcefieldParams ffparams = new ForcefieldParams(ff);
        ffparams.distDepDielect = params.getBool("DISTDEPDIELECT");
        ffparams.dielectric = params.getDouble("DIELECTCONST");
        ffparams.vdwMultiplier = params.getDouble("VDWMULT");
        ffparams.solvScale = params.getDouble("SOLVSCALE");
        ffparams.hElect = params.getBool("HELECT");
        ffparams.hVDW = params.getBool("HVDW");
        ffparams.shellDistCutoff = params.getDouble("SHELLDISTCUTOFF");
        
        if (params.getBool("USEPOISSONBOLTZMANN")) {
            ffparams.solvationForcefield = SolvationForcefield.PoissonBoltzmann;
        } else if (params.getBool("DOSOLVATIONE")) {
            ffparams.solvationForcefield = SolvationForcefield.EEF1;
        } else {
            ffparams.solvationForcefield = null;
        }
        
        EnvironmentVars.curEFcnGenerator = new EnergyFunctionGenerator(ffparams);
        
        // make the template library
        ResidueTemplateLibrary.Builder templateLibBuilder = new ResidueTemplateLibrary.Builder(ff)
			.clearRotamers()
            .addRotamers(params.readPath("ROTFILE"))
			.clearResidueEntropies()
            .addResidueEntropies(params.readPath("RESENTROPYFILE"));

        if (params.getBool("UseDunbrackRotamers")) {
			templateLibBuilder.addBackboneDependentRotamers(params.readPath("DUNBRACKROTFILE"));
		}

        // AAO 2016: load generic rotamer libraries
        for(String grotFile : params.searchParams("GROTFILE")) {
        	templateLibBuilder.addRotamers(params.readPath(grotFile));
        }

		EnvironmentVars.resTemplates = templateLibBuilder.build();

		// load rama data
        if (!params.getValue("RAMAGLYFILE").equalsIgnoreCase("none")) {
            RamachandranChecker.getInstance().readInputFiles(
                params.readPath("RAMAGLYFILE"),
                params.readPath("RAMAPROFILE"),
                params.readPath("RAMAGENFILE"),
                params.readPath("RAMAPREPROFILE")
            );
        }
        
        // TODO: pass this to the conf space somehow
        EnvironmentVars.alwaysIdealizeSidechainsAfterMutation = params.getBool("ALWAYSIDEALIZESIDECHAINSAFTERMUTATION");
    }

    //GMEC mut files
    public boolean hasGMECMutFile(){
        return ! params.getValue("GMECMutFile", "None").equalsIgnoreCase("None");
    }
    
    public GMECMutSpace parseGMECMutFile(ConfSpace confSpace){
        String mutFileName = params.getValue("GMECMutFile", "None");
        if(mutFileName.equalsIgnoreCase("None"))
            return null;
        else
            return new GMECMutSpace(mutFileName, confSpace);
    }
}
