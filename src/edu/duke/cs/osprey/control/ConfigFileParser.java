/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.dof.deeper.DEEPerSettings;
import edu.duke.cs.osprey.dof.deeper.RamachandranChecker;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.restypes.GenericResidueTemplateLibrary;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.pruning.PruningControl;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import edu.duke.cs.osprey.tools.StringParsing;
import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *
 * @author mhall44
 */
public class ConfigFileParser {
    //An object that parses configuration files and uses them to initialize objects needed
    //in various calculations (i.e. for data file loading, conf space definition, pruning, etc.)

    public ParamSet params = new ParamSet();

    boolean verbose = true;

    public ConfigFileParser(String[] args) {
        //parse all config files into params

        //check format of args
        if (!args[0].equalsIgnoreCase("-c")) {
            throw new RuntimeException("ERROR: bad arguments (should start with -c)");
        }

        params.addParamsFromFile(args[1]);//KStar.cfg file
        EnvironmentVars.setDataDir(params.getValue("DataDir"));

        for (int argNum = 3; argNum < args.length; argNum++)//System.cfg, etc.
        {
            params.addParamsFromFile(args[argNum]);
        }

        params.addDefaultParams();//We'll look for this in DataDir
    }

    DEEPerSettings setupDEEPer() {
        //Set up the DEEPerSettings object, including the PertSet (describes the perturbations)
        String runName = params.getValue("runName");

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
                params.getValue("PDBNAME")
        );

        dset.loadPertFile();//load the PertSet from its file
        return dset;
    }

    //HMN: create bound/unbound specific version 
    DEEPerSettings setupDEEPer(ArrayList<String> flexRes) {
        //Set up the DEEPerSettings object, including the PertSet (describes the perturbations)
        String runName = params.getValue("runName");

        DEEPerSettings dset = new DEEPerSettings(
                params.getBool("doPerturbations"),
                params.getValue("perturbationFile"),
                params.getBool("selectPerturbations"),
                params.getValue("startingPerturbationFile"),
                params.getBool("onlyStartingPerturbations"),
                params.getDouble("maxShearParam"),
                params.getDouble("maxBackrubParam"),
                params.getBool("selectLCAs"),
                flexRes,
                params.getValue("PDBNAME")
        );

        dset.loadPertFile();//load the PertSet from its file
        return dset;
    }

    private ArrayList<String[]> freeBBZoneTermini() {
        //Read the termini of the BBFreeBlocks, if any
        ArrayList<String[]> ans = new ArrayList<>();

        for (String rt : params.searchParams("BBFREEBLOCK")) {
            //So for example BBFREEBLOCK0 120 125 would mean make a BBFreeBlock for res 120-125
            //lexical ordering for blocks is OK
            String strandLimitsString = params.getValue(rt);

            String[] termini
                    = {StringParsing.getToken(strandLimitsString, 1),
                        StringParsing.getToken(strandLimitsString, 2)};

            ans.add(termini);
        }

        return ans;
    }

    //HMN: create  bound/unbound specific versions of freeBBZoneTermini
    private ArrayList<String[]> freeBBZoneTerminiBound() {
        return freeBBZoneTermini();
    }

    //HMN: freeBBZoneTermini for unbound 0 strand
    private ArrayList<String[]> freeBBZoneTerminiUnbound0() {
        //Read the termini of the BBFreeBlocks, if any
        ArrayList<String[]> ans = new ArrayList<>();

        for (String rt : params.searchParams("BBFREEBLOCK0")) {
            //So for example BBFREEBLOCK0 120 125 would mean make a BBFreeBlock for res 120-125
            //lexical ordering for blocks is OK
            String strandLimitsString = params.getValue(rt);

            String[] termini
                    = {StringParsing.getToken(strandLimitsString, 1),
                        StringParsing.getToken(strandLimitsString, 2)};

            ans.add(termini);
        }

        return ans;
    }

    //HMN: freeBBZoneTermini for unbound 1 strand
    private ArrayList<String[]> freeBBZoneTerminiUnbound1() {
        //Read the termini of the BBFreeBlocks, if any
        ArrayList<String[]> ans = new ArrayList<>();

        for (String rt : params.searchParams("BBFREEBLOCK1")) {
            //So for example BBFREEBLOCK0 120 125 would mean make a BBFreeBlock for res 120-125
            //lexical ordering for blocks is OK
            String strandLimitsString = params.getValue(rt);

            String[] termini
                    = {StringParsing.getToken(strandLimitsString, 1),
                        StringParsing.getToken(strandLimitsString, 2)};

            ans.add(termini);
        }

        return ans;
    }

    private ArrayList<String[]> moveableStrandTermini() {
        //Read the strands that are going to translate and rotate
        //Let's say they can do this regardless of what doMinimize says (that's for sidechains)
        ArrayList<String[]> ans = new ArrayList<>();

        for (String rt : params.searchParams("STRANDROTTRANS")) {
            if (params.getBool(rt)) {

                //So rt = STRANDROTTRANS0 here means strand 0 should translate & rotate
                //OK to go through these params in lexical ordering
                String strandNum = rt.substring(14);
                String strandLimitsString = params.getValue("STRAND" + strandNum);

                String[] termini
                        = {StringParsing.getToken(strandLimitsString, 1),
                            StringParsing.getToken(strandLimitsString, 2)};

                ans.add(termini);
            }
        }

        return ans;
    }

    //HMN: Make sure we only get moveable strands for Bound or Unbound states that 
    ///we are interested in
    private ArrayList<String[]> moveableStrandTerminiBound() {
        return moveableStrandTermini();
    }

    private ArrayList<String[]> moveableStrandTerminiUnboud0() {
        //Read the strands that are going to translate and rotate
        //Let's say they can do this regardless of what doMinimize says (that's for sidechains)
        ArrayList<String[]> ans = new ArrayList<>();

        for (String rt : params.searchParams("STRANDROTTRANS0")) {
            if (params.getBool(rt)) {
                //So rt = STRANDROTTRANS0 here means strand 0 should translate & rotate
                //OK to go through these params in lexical ordering
                String strandNum = rt.substring(14);
                String strandLimitsString = params.getValue("STRAND" + strandNum);

                String[] termini
                        = {StringParsing.getToken(strandLimitsString, 1),
                            StringParsing.getToken(strandLimitsString, 2)};

                ans.add(termini);
            }
        }
        return ans;
    }

    private ArrayList<String[]> moveableStrandTerminiUnboud1() {
        //Read the strands that are going to translate and rotate
        //Let's say they can do this regardless of what doMinimize says (that's for sidechains)
        ArrayList<String[]> ans = new ArrayList<>();

        for (String rt : params.searchParams("STRANDROTTRANS1")) {
            if (params.getBool(rt)) {
                //So rt = STRANDROTTRANS0 here means strand 0 should translate & rotate
                //OK to go through these params in lexical ordering
                String strandNum = rt.substring(14);
                String strandLimitsString = params.getValue("STRAND" + strandNum);

                String[] termini
                        = {StringParsing.getToken(strandLimitsString, 1),
                            StringParsing.getToken(strandLimitsString, 2)};

                ans.add(termini);
            }
        }
        return ans;
    }

    //creation of objects needed in calculations like GMEC and K*
    public SearchProblem getSearchProblem() {//this version is for a single search problem...can modify for
        //additional states (unbound, etc.)

        String name = params.getValue("RUNNAME");

        ArrayList<String> flexRes = getFlexRes();
        ArrayList<ArrayList<String>> allowedAAs = getAllowedAAs();

        if (flexRes.size() != allowedAAs.size()) {
            throw new RuntimeException("ERROR: Number of flexible positions different in flexible residue "
                    + "and allowed AA type parameters!");
        }
        if (verbose) {
            System.out.println("CREATING SEARCH PROBLEM.  NAME: " + name);
        }
        ArrayList<String[]> moveableStrands = moveableStrandTermini();
        ArrayList<String[]> freeBBZones = freeBBZoneTermini();
        DEEPerSettings dset = setupDEEPer();

        return new SearchProblem(name, params.getValue("PDBNAME"),
                flexRes, allowedAAs,
                params.getBool("AddWT"),
                params.getBool("AddWTRots"),
                params.getBool("doMinimize"),
                params.getBool("UseEPIC"),
                new EPICSettings(params),
                params.getBool("UseTupExp"),
                dset, moveableStrands, freeBBZones,
                params.getBool("useEllipses"),
                params.getBool("useERef"),
                params.getBool("AddResEntropy"));

    }


    public ArrayList<String> getFlexRes() {
        //list of flexible residues.  PDB-based residue numbers
        //we'll include all flexible residues: for compatibility (MAY BE TEMPORARY),
        //all residues in a "StrandMut" record will be included here
        //"StrandMutNums" means something different and thus will be excluded
        ArrayList<String> flexResList = new ArrayList<>();

        int numStrands = params.searchParams("STRANDMUT").size()
                - params.searchParams("STRANDMUTNUMS").size();

        //must go through the strands in order to get the right residues 
        for (int smNum = 0; smNum < numStrands; smNum++) {//must do these in order
            //so we get the right residues
            String param = "STRANDMUT" + smNum;//STRANDMUT0, etc

            String resListString = params.getValue(param);
            StringTokenizer tokenizer = new StringTokenizer(resListString);

            while (tokenizer.hasMoreTokens()) {
                flexResList.add(tokenizer.nextToken());
            }
        }

        return flexResList;
    }

    public int getMutableStrand() {
        int mutStrand = -1;
        boolean mutStrandSet = false;
        //handle better later (now assuming old-style numbering)...
        for (int str = 0; str < 2; str++) {
            ArrayList<String> resAllowedRecords = params.searchParams("RESALLOWED" + str);
            int numRecordsInStrand = resAllowedRecords.size();
            //must go through residues in numerical order
            for (int recNum = 0; recNum < numRecordsInStrand; recNum++) {
                String param = "RESALLOWED" + str + "_" + recNum;
                String allowedAAString = params.getValue(param);
                if ((allowedAAString.trim().split(" ").length > 1)
                        || ((allowedAAString.trim().split(" ").length == 1) && (!allowedAAString.trim().split(" ")[0].isEmpty()))) {
                    if (mutStrandSet && (str != mutStrand)) {
                        throw new RuntimeException("ERROR: DEE.cfg must have mutations only on strand0");
                    } else {
                        mutStrand = str;
                        mutStrandSet = true;
                    }
                }
            }
        }
        return mutStrand;
    }

    public ArrayList<ArrayList<String>> getAllowedAAs() {
        //List allowed AA types for each flexible position
        //We can accept either RESALLOWED0_0 (for flexible res 0 of strand 0)
        //or RESALLOWED255 (for residue with PDB number 255)
        //but which one is used should be consistent across all residues
        ArrayList<String> resAllowedRecords = params.searchParams("RESALLOWED");

        if (resAllowedRecords.isEmpty())//no flexible residues
        {
            return new ArrayList<ArrayList<String>>();
        } else {
            boolean usingStrandFormat = resAllowedRecords.get(0).contains("_");
            for (int rec = 1; rec < resAllowedRecords.size(); rec++) {
                if (usingStrandFormat != resAllowedRecords.get(rec).contains("_")) {
                    throw new RuntimeException("ERROR: Inconsistent formatting of resAllowed records"
                            + " (should be all by PDB residue number or all by strand)");
                }
            }

            if (usingStrandFormat) {
                return getAllowedAAsByStrand();
            } else {
                return getAllowedAAsByPDBResNum();
            }
        }
    }

    ArrayList<ArrayList<String>> getAllowedAAsByStrand() {
        //we'll go through all resAllowed records, ordering first by strand num
        //then pos num in strand

        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();

        //handle better later (now assuming old-style numbering)...
        for (int str = 0; true; str++) {
            ArrayList<String> resAllowedRecords = params.searchParams("RESALLOWED" + str);
            int numRecordsInStrand = resAllowedRecords.size();

            //must go through residues in numerical order
            for (int recNum = 0; recNum < numRecordsInStrand; recNum++) {
                String param = "RESALLOWED" + str + "_" + recNum;

                String allowedAAString = params.getValue(param);

                //parse AA types from allowedAAString
                ArrayList<String> resAllowedAAs = new ArrayList<>();
                StringTokenizer tokenizer = new StringTokenizer(allowedAAString);

                while (tokenizer.hasMoreTokens()) {
                    resAllowedAAs.add(tokenizer.nextToken());
                }

                allowedAAs.add(resAllowedAAs);
            }

            if (numRecordsInStrand == 0)//finished with strands that have flexible residues
            {
                break;
            }
        }

        return allowedAAs;
    }

    ArrayList<ArrayList<String>> getAllowedAAsByPDBResNum() {
        //we'll go through all resAllowed records, ordering first by strand num
        //then pos num in strand

        ArrayList<ArrayList<String>> allowedAAs = new ArrayList<>();
        ArrayList<String> flexRes = getFlexRes();

        for (int flexResNum = 0; flexResNum < flexRes.size(); flexResNum++) {
            String param = "RESALLOWED" + flexRes.get(flexResNum);
            String allowedAAString = params.getValue(param);

            //parse AA types from allowedAAString
            ArrayList<String> resAllowedAAs = new ArrayList<>();
            StringTokenizer tokenizer = new StringTokenizer(allowedAAString);

            while (tokenizer.hasMoreTokens()) {
                resAllowedAAs.add(tokenizer.nextToken());
            }

            allowedAAs.add(resAllowedAAs);
        }

        return allowedAAs;
    }

    public PruningControl setupPruning(SearchProblem searchSpace, double pruningInterval, boolean useEPIC, boolean useTupExp) {

        //setup pruning.  Conformations in searchSpace more than (Ew+Ival) over the GMEC are liable to pruning
        //initialize the pruning matrix for searchSpace, if not already initialized
        //or if pruningInterval lower (so it may have pruned tuples that shouldn't
        //be pruned with our new pruningInterval)
        boolean initPruneMat = false;
        if (searchSpace.pruneMat == null) {
            initPruneMat = true;
        } else if (searchSpace.pruneMat.getPruningInterval() < pruningInterval) {
            initPruneMat = true;
        }

        if (initPruneMat) {
            searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, pruningInterval);
        }

        return new PruningControl(searchSpace, pruningInterval, params.getBool("TYPEDEP"),
                params.getDouble("BOUNDSTHRESH"), params.getInt("ALGOPTION"),
                params.getBool("USEFLAGS"),
                params.getBool("USETRIPLES"), false, useEPIC, useTupExp,
                params.getDouble("STERICTHRESH"));//FOR NOW NO DACS

    }

    //loading of data files
    //residue templates, rotamer libraries, forcefield parameters, and Ramachandran data
    public void loadData() {

        boolean usePoissonBoltzmann = params.getBool("USEPOISSONBOLTZMANN");
        boolean useEEF1 = params.getBool("DOSOLVATIONE") && (!usePoissonBoltzmann);

        //a lot of this depends on forcefield type so figure that out first
        //general forcefield data loaded into the ForcefieldParams in EnvironmentVars
        ForcefieldParams curForcefieldParams = new ForcefieldParams(
                params.getValue("Forcefield"),
                params.getBool("DISTDEPDIELECT"),
                params.getDouble("DIELECTCONST"),
                params.getDouble("VDWMULT"),
                useEEF1,//Only EEF1 solvation is part of the forcefield (P-B calls Delphi)
                params.getDouble("SOLVSCALE"),
                params.getBool("HELECT"),
                params.getBool("HVDW"));

        EnvironmentVars.curEFcnGenerator = new EnergyFunctionGenerator(
                curForcefieldParams,
                params.getDouble("SHELLDISTCUTOFF"),
                usePoissonBoltzmann);

        String[] resTemplateFiles = getResidueTemplateFiles(curForcefieldParams.forcefld);

        GenericResidueTemplateLibrary resTemplates = new GenericResidueTemplateLibrary(resTemplateFiles, curForcefieldParams);

        //load template coordinates (necessary for all residues we might want to mutate to)
        //these will be matched to templates
        resTemplates.loadTemplateCoords("all_amino_coords.in");

        //load rotamer libraries; the names of residues as they appear in the rotamer library file will be matched to templates
        boolean dunbrackRots = params.getBool("UseDunbrackRotamers");
        // PGC 2015: Always load the Lovell Rotamer Library.
        resTemplates.loadRotamerLibrary(params.getValue("ROTFILE"), false);//see below; also gRotFile0 etc
        if (dunbrackRots) { // Use the dunbrack rotamer library
            resTemplates.loadRotamerLibrary(params.getValue("DUNBRACKROTFILE"), true);//see below; also gRotFile0 etc
        }

//        resTemplates.loadResEntropy(params.getValue("RESENTROPYFILE"));

        //let's make D-amino acid templates by inverting the L-amino acid templates 
        resTemplates.makeDAminoAcidTemplates();

        EnvironmentVars.resTemplates = resTemplates;

        String ramaGlyFile = params.getValue("RAMAGLYFILE");

        if (!ramaGlyFile.equalsIgnoreCase("none")) {
            String ramaFiles[] = {EnvironmentVars.getDataDir() + ramaGlyFile,
                EnvironmentVars.getDataDir() + params.getValue("RAMAPROFILE"),
                EnvironmentVars.getDataDir() + params.getValue("RAMAGENFILE"),
                EnvironmentVars.getDataDir() + params.getValue("RAMAPREPROFILE")

            };
            RamachandranChecker.getInstance().readInputFiles(ramaFiles);
        }

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

    String[] getResidueTemplateFiles(ForcefieldParams.FORCEFIELD forcefld) {
        //return names of residue template files

        //template file names are currently fixed
        String aaFilename = null, aaNTFilename = null, aaCTFilename = null, grFilename = null;

        switch (forcefld) {
            case AMBER:
                //KER: This is for the standard amber parameters
                aaFilename = "all_amino94.in";
                aaNTFilename = "all_aminont94.in";
                aaCTFilename = "all_aminoct94.in";
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
                aaFilename = "all_amino_charmm19_neutral.in";
                aaNTFilename = "all_amino_charmm19_neutral_nt.in";
                aaCTFilename = "all_amino_charmm19_neutral_ct.in";
                grFilename = "all_nuc_and_gr_charmm.in";
                break;
            case CHARMM19:
                aaFilename = "all_amino_charmm19.in";
                aaNTFilename = "all_amino_charmm19_nt.in";
                aaCTFilename = "all_amino_charmm19_ct.in";
                grFilename = "all_nuc_and_gr_charmm.in";
                break;
            default:
                System.out.println("FORCEFIELD not recognized...Exiting");
                System.exit(0);
        }

        return new String[]{
            aaFilename, aaNTFilename, aaCTFilename, grFilename
        };
    }

    //HMN: get search problems for classic multistate bound vs. unbound states
    ///Bound, unbound1, unbound2
    public SearchProblem[] getMSDSearchProblems() {
        String name = params.getValue("RUNNAME");

        //make sure we have two strands
        if (!(params.getInt("NUMOFSTRANDS") == 2)) {
            throw new RuntimeException("Cannot get 3 Search Problems since numOfStrans != 2");
        }

        int numMut0 = Integer.parseInt(params.getValue("STRANDMUTNUMS").split(" ")[0]);
        int numMut1 = Integer.parseInt(params.getValue("STRANDMUTNUMS").split(" ")[1]);

        ArrayList<String> flexResBound = getFlexRes();
        ArrayList<ArrayList<String>> allowedAAsBound = getAllowedAAs();

        if (!(flexResBound.size() == numMut0 + numMut1)) {
            throw new RuntimeException("ERROR: strandMutNums does not equal size of strandMut0 and strandMut1");
        }

        if (flexResBound.size() != allowedAAsBound.size()) {
            throw new RuntimeException("ERROR: Number of flexible positions different in flexible residue "
                    + "and allowed AA type parameters!");
        }

        //Get flexibleRes and allowedAAs for each unbound strand
        ArrayList<String> flexResUnbound0 = new ArrayList<>();
        ArrayList<ArrayList<String>> allowedAAsUnbound0 = new ArrayList<>();

        ArrayList<String> flexResUnbound1 = new ArrayList<>();;
        ArrayList<ArrayList<String>> allowedAAsUnbound1 = new ArrayList<>();

        int iter = 0;
        for (int i = 0; i < numMut0; i++) {
            flexResUnbound0.add(flexResBound.get(i));
            allowedAAsUnbound0.add(allowedAAsBound.get(i));
            iter++;
        }
        for (int j = iter; j < iter + numMut1; j++) {
            flexResUnbound1.add(flexResBound.get(j));
            allowedAAsUnbound1.add(allowedAAsBound.get(j));
        }

        if ((flexResUnbound0.size() != allowedAAsUnbound0.size()) || (flexResUnbound1.size() != allowedAAsUnbound1.size())) {
            throw new RuntimeException("ERROR: Number of flexible positions at a strand is different "
                    + "in flexible residue and allowed AA type parameters!");
        }

        ArrayList<String[]> moveableStrandsBound = moveableStrandTerminiBound();
        ArrayList<String[]> moveableStrandsUnbound0 = moveableStrandTerminiUnboud0();
        ArrayList<String[]> moveableStrandsUnbound1 = moveableStrandTerminiUnboud1();

        ArrayList<String[]> freeBBZonesBound = freeBBZoneTerminiBound();
        ArrayList<String[]> freeBBZonesUnbound0 = freeBBZoneTerminiUnbound0();
        ArrayList<String[]> freeBBZonesUnbound1 = freeBBZoneTerminiUnbound1();

        DEEPerSettings dsetBound = setupDEEPer(flexResBound);
        DEEPerSettings dsetUnbound0 = setupDEEPer(flexResUnbound0);
        DEEPerSettings dsetUnbound1 = setupDEEPer(flexResUnbound1);

        int[] strand0Termini = new int[2];
        strand0Termini[0] = Integer.parseInt(params.getValue("STRAND0").split(" ")[0]);
        strand0Termini[1] = Integer.parseInt(params.getValue("STRAND0").split(" ")[1]);
        int[] strand1Termini = new int[2];
        strand1Termini[0] = Integer.parseInt(params.getValue("STRAND1").split(" ")[0]);
        strand1Termini[1] = Integer.parseInt(params.getValue("STRAND1").split(" ")[1]);

        SearchProblem[] searchProblems = new SearchProblem[3];
        String pdbFile = params.getValue("PDBNAME");
        String pdbName = pdbFile.split(".pdb")[0];

        Molecule molComplex = PDBFileReader.readPDBFileBetweenTermini(pdbFile, strand0Termini[0], strand0Termini[1], strand1Termini[0], strand1Termini[1]);

        //get current working directory and write out file
        File cwdFile = new File("");
        String cwd = cwdFile.getAbsolutePath() + "/";
        String pdbComplex = pdbName + "_complex.pdb";
        PDBFileWriter.writePDBFile(molComplex, pdbComplex);

        searchProblems[0] = new SearchProblem(name + "_Bound", pdbComplex,
                flexResBound, allowedAAsBound,
                params.getBool("AddWT"),
                params.getBool("AddWTRots"),
                params.getBool("doMinimize"),
                params.getBool("UseEPIC"),
                new EPICSettings(params),
                params.getBool("UseTupExp"),
                dsetBound, moveableStrandsBound, freeBBZonesBound,
                params.getBool("useEllipses"),
                params.getBool("useERef"),
                params.getBool("AddResEntropy")
        );
        //Delete file 
        File pdbComplexFile = new File(cwd + pdbComplex);
        try {
            pdbComplexFile.delete();
        } catch (Exception e) {
            System.out.println("File cannot be deleted");
        }

        Molecule molLigand = PDBFileReader.readPDBFileBetweenTermini(pdbFile, strand0Termini[0], strand0Termini[1]);
        String pdbLigand = pdbName + "_strand1.pdb";
        PDBFileWriter.writePDBFile(molLigand, pdbLigand);

        searchProblems[1] = new SearchProblem(name + "_strand1", pdbLigand,
                flexResUnbound0, allowedAAsUnbound0,
                params.getBool("AddWT"),
                params.getBool("AddWTRots"),
                params.getBool("doMinimize"),
                params.getBool("UseEPIC"),
                new EPICSettings(params),
                params.getBool("UseTupExp"),
                dsetUnbound0, moveableStrandsUnbound0, freeBBZonesUnbound0,
                params.getBool("useEllipses"),
                params.getBool("useERef"),
                params.getBool("AddResEntropy")
        );
        //Delete file 
        File pdbLigandFile = new File(cwd + pdbLigand);
        try {
            pdbLigandFile.delete();
        } catch (Exception e) {
            System.out.println("File cannot be deleted");
        }

        Molecule molProtein = PDBFileReader.readPDBFileBetweenTermini(pdbFile, strand1Termini[0], strand1Termini[1]);
        String pdbProtein = pdbName + "_strand2.pdb";
        PDBFileWriter.writePDBFile(molProtein, pdbProtein);

        searchProblems[2] = new SearchProblem(name + "_strand2", pdbProtein,
                flexResUnbound1, allowedAAsUnbound1,
                params.getBool("AddWT"),
                params.getBool("AddWTRots"),
                params.getBool("doMinimize"),
                params.getBool("UseEPIC"),
                new EPICSettings(params),
                params.getBool("UseTupExp"),
                dsetUnbound1, moveableStrandsUnbound1, freeBBZonesUnbound1,
                params.getBool("useEllipses"),
                params.getBool("useERef"),
                params.getBool("AddResEntropy")
        );
        //Delete file 
        File pdbProteinFile = new File(cwd + pdbProtein);
        try {
            pdbProteinFile.delete();
        } catch (Exception e) {
            System.out.println("File cannot be deleted");
        }

        return searchProblems;
    }

    public int[] getNumMutPerStrand() {
        int numMut0 = Integer.parseInt(params.getValue("STRANDMUTNUMS").split(" ")[0]);
        int numMut1 = Integer.parseInt(params.getValue("STRANDMUTNUMS").split(" ")[1]);
        int[] numMutPerStrand = new int[2];
        numMutPerStrand[0] = numMut0;
        numMutPerStrand[1] = numMut1;
        return numMutPerStrand;
    }

// Getter function for the params.
    public ParamSet getParams() {
        return this.params;
    }

}
