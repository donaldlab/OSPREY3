/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.astar.comets.COMETSTree;
import edu.duke.cs.osprey.astar.comets.LME;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.StringParsing;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *
 * @author mhall44
 */
class COMETSDoer {
    
    COMETSTree tree;//The tree used for the COMETS search
    int numSeqsWanted;//How many sequences to enumerate
    
    LME objFcn;//objective function for the COMETS search
    LME[] constraints;//constraints for the COMETS search
    int numStates;//number of states considered
    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null;//AA types allowed at each mutable position
    
    ArrayList<ArrayList<Integer>> mutable2StatePosNums = new ArrayList<>();
    //For each state, a list of which flexible positions are mutable
    //these will be listed directly in Multistate.cfg under "STATEMUTRES0" etc.
    
    String stateArgs[][];//"GMEC search" arguments for each state
    
    public COMETSDoer (String args[]){
        //fill in all the settings
        //each state will have its own config file parser
        
        System.out.println();
        System.out.println("Performing multistate A*");
        System.out.println();
        
        //check format of args
        if(!args[0].equalsIgnoreCase("-c"))
            throw new RuntimeException("ERROR: bad arguments (should start with -c)");
        
        
        ParamSet sParams = new ParamSet();
        
        sParams.addParamsFromFile(args[1]);//read KStar.cfg parameters
        sParams.addParamsFromFile(args[3]);//read multistate parameters
        sParams.addDefaultParams();
        
        String defaultCFGName = args[1];//Each state could have its own KStar.cfg, but this is the default
        
        numStates = sParams.getInt("NUMSTATES");
        numTreeLevels = sParams.getInt("NUMMUTRES");
        int numConstr = sParams.getInt("NUMCONSTR");
        
        stateArgs = new String[numStates][];
        
        objFcn = new LME( sParams.getValue("OBJFCN"), numStates );
        
        constraints = new LME[numConstr];
        for(int constr=0; constr<numConstr; constr++)
            constraints[constr] = new LME( sParams.getValue("CONSTR"+constr), numStates );
        
        
        SearchProblem[] stateSP = new SearchProblem[numStates];    

        
        System.out.println();
        System.out.println("Preparing matrices and search problems for multistate A*");
        System.out.println();
        
        for(int state=0; state<numStates; state++){
            mutable2StatePosNums.add( stateMutablePos(state,sParams,numTreeLevels) );
            stateSP[state] = makeStateSearchProblem(state, sParams, defaultCFGName);
            
            System.out.println();
            System.out.println("State "+state+" matrices ready.");
            System.out.println();
        }
        
        //we can collect info on allowed sequences from any state (we choose state 0)
        //but we'll check that they all agree

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)
        int numMaxMut = sParams.getInt("NUMMAXMUT");
        String[] wtSeq = null;
        if(numMaxMut>-1){
            wtSeq = parseWTSeq( sParams.getValue("WTSEQ"), numTreeLevels );
        }

        numSeqsWanted = sParams.getInt("NUMSEQS");
        boolean outputGMECStructs = sParams.getBool("OutputStateGMECStructs");
        
        tree = new COMETSTree(numTreeLevels, objFcn, constraints, 
            AATypeOptions, numMaxMut, wtSeq, numStates, stateSP, 
            mutable2StatePosNums, outputGMECStructs);
    }
    
    
    private String[] parseWTSeq(String seq, int numTreeLevels){
        //WT sequence (for purposes of capping # mutations) specified like ARG LYS...,
        StringTokenizer st = new StringTokenizer(seq);
        if(st.countTokens()!=numTreeLevels)
            throw new RuntimeException("ERROR: wrong number of residues in WT seq: "+seq);
        
        String wt[] = new String[numTreeLevels];
        
        for(int level=0; level<numTreeLevels; level++)
            wt[level] = st.nextToken();
        
        return wt;
    }
    
    
    private SearchProblem makeStateSearchProblem(int state, ParamSet sParams, String defaultCFGName){
        //read state-specific configuration files and create a search problem object
        //defaultCFGName is for KStar.cfg
        String cfgName = sParams.getValue("STATEKSFILE"+state);
        if(cfgName.equalsIgnoreCase("DEFAULT"))
            cfgName = defaultCFGName;

        //We expect input like
        //STATECFGFILES0 System0.cfg DEE0.cfg
        String stateConfigFiles = sParams.getValue("STATECFGFILES"+state);
        String stateSysFile = StringParsing.getToken(stateConfigFiles, 1);
        String stateDEEFile = StringParsing.getToken(stateConfigFiles, 2);

        stateArgs[state] = new String[] {"-c", cfgName, "n/a", stateSysFile, stateDEEFile};

        ConfigFileParser stateCFGP = new ConfigFileParser(stateArgs[state]);
        stateCFGP.loadData();
        SearchProblem searchProb = stateCFGP.getSearchProblem();

        if ( stateCFGP.params.getBool("doMinimize") && (!searchProb.useTupExpForSearch) ) {
            throw new RuntimeException("ERROR: COMETS requires LUTE to do continuous flexibility");
        }
        
        precomputeMatrices(state, stateCFGP, searchProb);
        
        handleAATypeOptions(state, stateCFGP);
                
        return searchProb;
    }
    
    
    private void handleAATypeOptions(int state, ConfigFileParser cfgP){
        //Given the config file parser for a state, make sure AATypeOptions
        //matches the allowed AA types for this state
        
        ArrayList<ArrayList<String>> stateAAOptions = cfgP.getAllowedAAs();
        
        Molecule wtMolec = PDBFileReader.readPDBFile( cfgP.params.getValue("PDBName") );
        ArrayList<String> flexRes = cfgP.getFlexRes();
        
        
        if(AATypeOptions == null)//still null...fill in based on this state
            AATypeOptions = new ArrayList<>();
        
        for(int mutPos=0; mutPos<numTreeLevels; mutPos++){

            int flexPos = mutable2StatePosNums.get(state).get(mutPos);
            
            ArrayList<String> statePosOptions = stateAAOptions.get(flexPos);
            if(cfgP.params.getBool("AddWT")){
                Residue res = wtMolec.getResByPDBResNumber( flexRes.get(flexPos) );
                if( ! StringParsing.containsIgnoreCase(statePosOptions, res.template.name) )
                    statePosOptions.add(res.template.name);
            }
            
            if(AATypeOptions.size()<=mutPos){//need to fill in based on this state
                AATypeOptions.add(statePosOptions);
            }
            else {//check to make sure they match
                ArrayList<String> posOptions = AATypeOptions.get(mutPos);

                if(posOptions.size() != statePosOptions.size()){
                    throw new RuntimeException("ERROR: Current state has "+
                            statePosOptions.size()+" AA types allowed for position "+mutPos
                            +" compared to "+posOptions.size()+" for previous states");
                }

                for(int a=0; a<posOptions.size(); a++){
                    String aa1 = posOptions.get(a);
                    String aa2 = statePosOptions.get(a);

                    if( ! aa1.equalsIgnoreCase(aa2) ){
                        throw new RuntimeException("ERROR: Current state has AA type "+
                            aa2+" where previous states have AA type "+aa1+", at position "+mutPos);
                    }
                }
            }
        }
    }
    
    
    private ArrayList<Integer> stateMutablePos(int state, ParamSet sParams, int numTreeLevels){
        //read the list of which of this state's flexible positions are mutable
        ArrayList<Integer> m2s = new ArrayList<>();
        String stateMutRes = (String)sParams.getValue("STATEMUTRES"+state);
        StringTokenizer st = new StringTokenizer(stateMutRes);

        while(st.hasMoreTokens())
            m2s.add(Integer.valueOf(st.nextToken()));

        if(m2s.size()!=numTreeLevels){
            throw new RuntimeException("ERROR: SeqTree has "+numTreeLevels+" mutable positions "
                    +" but "+m2s.size()+" are listed for state "+state);
        }
            
        return m2s;
    }
    
    
    public void calcBestSequences(){
                    
        System.out.println("Performing multistate A*");
        
            
        //how many sequences to enumerate

        long startAStarTime = System.currentTimeMillis();


        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
            int seq[] = tree.nextConf();//this will find the best sequence and print it
            if(seq == null)//empty sequence...indicates no more sequence possibilities
                break;
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));
        
        
        //DEBUG!!!  Checking by exhaustive search
        //exhaustiveMultistateSearch();
    }
        
        
    void precomputeMatrices(int state, ConfigFileParser cfgP, SearchProblem sp){
        //Precompute energy and pruning matrices for a state
        //as if computing a GMEC
        
        cfgP.params.setValue("TYPEDEP", "TRUE");//we need to do type-dependent pruning
        
        GMECFinder gf = new GMECFinder(cfgP); 
        gf.searchSpace = sp;//Using our search problem here instead of full GMEC calculation
        
        double ival = chooseIVal(state,cfgP,gf);
        
        gf.precomputeMatrices(ival);//no Ew...just looking at GMECs
        //and for any valid sequence the GMEC shouldn't be pruned (unless ival manually set lower)
    }
    
    
    double chooseIVal(int state, ConfigFileParser cfgP, GMECFinder gf){
        //choose the iMinDEE interval we want to use for a state
        
        if( ! gf.useTupExp )//discrete flexibility & pairwise E-fcn...iVal can be 0
            return 0;
        
        double iVal = cfgP.params.getDouble("COMETSIVAL");
        //If iVal provided manually, use that
        
        if(Double.isNaN(iVal)){//need to choose ival
            
            
            if(gf.EFullConfOnly)//do only steric pruning
                iVal = Double.POSITIVE_INFINITY;
            else {
            
                //If there's a stability constr for this state, compute lowest pairwise bound
                //and then we know we can use I = max allowed seq GMEC - lowest pairwise bound
                //Else, use only steric pruning (ival = infinity)
                
                double maxStateGMEC = Double.POSITIVE_INFINITY;//max allowed seq GMEC for state

                for(LME constr : constraints){
                    double coeffs[] = constr.getCoeffs();

                    boolean isStabilityConstr = true;//is this a stability constraint for the current state?

                    for(int state2=0; state2<numStates; state2++){
                        if(state2==state && Math.abs(coeffs[state2]-1) > 1e-8)
                            isStabilityConstr = false;
                        else if (state2!=state && Math.abs(coeffs[state2]) > 1e-8)
                            isStabilityConstr = false;
                    }

                    if(isStabilityConstr){
                        maxStateGMEC = -constr.getConstTerm();
                        break;
                    }
                }

                
                if(maxStateGMEC==Double.POSITIVE_INFINITY)//no stability constr...can't bound iVal
                    iVal = Double.POSITIVE_INFINITY;
                else{

                    //compute pairwise lowest bound

                    gf.searchSpace.loadEnergyMatrix();

                    //allocate pruning matrices, but don't prune since ival unknown
                    cfgP.setupPruning(gf.searchSpace, Double.POSITIVE_INFINITY, false, false);
                    double lowestBound = gf.lowestPairwiseBound(gf.searchSpace);
                    iVal = maxStateGMEC - lowestBound;
                }
            }
        }
        
        System.out.println("State "+state+" is using ival "+iVal);
        return iVal;
    }
        
    
    
    //For quality control, it's good to be able to check COMETS results by exhaustive search
    void exhaustiveMultistateSearch(){
        
        System.out.println();
        System.out.println("CHECKING COMETS RESULT BY EXHAUSTIVE SEARCH");
        System.out.println();
                
        
        ArrayList<String[]> seqList = listAllSeqs();
        int numSeqs = seqList.size();
        double stateGMECs[][] = new double[numSeqs][numStates];
        
        for(int state=0; state<numStates; state++){
            for(int seqNum=0; seqNum<numSeqs; seqNum++){
                stateGMECs[seqNum][state] = calcStateGMEC(state, seqList.get(seqNum));
            }
        }
        
        //now find the best sequence and obj fcn value
        int topSeqNum = -1;
        double bestSeqScore = Double.POSITIVE_INFINITY;
        
        for(int seqNum=0; seqNum<numSeqs; seqNum++){
            boolean constrSatisfied = true;
            for(LME constr : constraints){
                if(constr.eval(stateGMECs[seqNum]) > 0)
                    constrSatisfied = false;
            }
            
            if(constrSatisfied){
                double seqScore = objFcn.eval(stateGMECs[seqNum]);
                if(seqScore < bestSeqScore){
                    bestSeqScore = seqScore;
                    topSeqNum = seqNum;
                }
            }
        }
        
        System.out.println();
        if(topSeqNum == -1)
            System.out.println("EXHAUSTIVE SEARCH FINDS NO CONSTR-SATISFYING SEQUENCES");
        else{
            System.out.println("EXHAUSTIVE MULTISTATE BEST SCORE: "+bestSeqScore+" SEQUENCE: ");
            for(String aa : seqList.get(topSeqNum))
                System.out.println(aa);
        }
        System.out.println();
    }
    
    
    private ArrayList<String[]> listAllSeqs(){
        //List all possible sequence for the mutable residues,
        //based on AATypeOptions
        return listAllSeqsHelper(0);
    }
    
    
    private ArrayList<String[]> listAllSeqsHelper(int mutPos){
        //List all partial sequences for the subset of mutable positions
        //starting at mutPos and going to the last mutable position
        ArrayList<String[]> ans = new ArrayList<>();
        
        if(mutPos == numTreeLevels){
            ans.add(new String[0]);
        }
        else {
            ArrayList<String[]> subList = listAllSeqsHelper(mutPos+1);
            for(String AAType : AATypeOptions.get(mutPos)){
                for(String[] subSeq : subList){
                    String seq[] = new String[numTreeLevels-mutPos];
                    System.arraycopy(subSeq, 0, seq, 1, numTreeLevels-mutPos-1);
                    seq[0] = AAType;
                    ans.add(seq);
                }
            }
        }
        
        return ans;
    }
    
    
    private double calcStateGMEC(int state, String[] AATypes){
        //Calculate the GMEC for the specified state for sequence AATypes
        ConfigFileParser cfp = new ConfigFileParser(stateArgs[state]);
        cfp.loadData();
        
        //set up for particular AA types
        cfp.params.setValue("RUNNAME", "EXHAUSTIVE_SEQ_"+System.currentTimeMillis());
        
        //also want to restrict to only this seq: addWT = false
        cfp.params.setValue("ADDWT", "false");        
        
        //DEBUG!!!
        //cfp.params.setValue("USETUPEXP", "false");
        
        String seq[] = new String[cfp.getFlexRes().size()];
        for(int mutPos=0; mutPos<numTreeLevels; mutPos++)
            seq[mutable2StatePosNums.get(state).get(mutPos)] = AATypes[mutPos];
        
        
        //matching format to ConfigFileParser.getAllowedAAs
        int posCount = 0;
        
        for(int str=0; str<10; str++){
            ArrayList<String> resAllowedRecords = cfp.params.searchParams("RESALLOWED"+str);
            int numRecordsInStrand = resAllowedRecords.size();
            
            //must go through residues in numerical order
            for(int recNum=0; recNum<numRecordsInStrand; recNum++){
                
                if(seq[posCount] != null){
                    String param = "RESALLOWED" + str + "_" + recNum;
                    cfp.params.setValue(param, seq[posCount]);
                }
                
                posCount++;
            }
        }
        
        GMECFinder gf = new GMECFinder(cfp);
        return gf.calcGMEC();
    }
    
}
