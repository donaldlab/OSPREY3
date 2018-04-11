/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.util.ArrayList;
import java.util.StringTokenizer;

import edu.duke.cs.osprey.astar.comets.COMETSTree;
import edu.duke.cs.osprey.astar.comets.LME;
import edu.duke.cs.osprey.astar.ewakstar.EWAKLME;
import edu.duke.cs.osprey.astar.ewakstar.EWAKStarTree;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.ewakstar.EWAKStarLowEnergyFinder;
import edu.duke.cs.osprey.gmec.GMECFinder;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;
import edu.duke.cs.osprey.tools.FileTools.PathRoot;
import edu.duke.cs.osprey.tools.StringParsing;

/**
 *
 * @author lowegard (adapted from COMETSDoer.java by mhall44)
 */
public class EWAKStarDoer {

    EWAKStarTree tree;//The tree used for the COMETS search
    int numSeqsWanted;//How many sequences to enumerate

    EWAKLME objFcn;//objective function for the COMETS search
    int numStates;//number of states considered
    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null;//AA types allowed at each mutable position

    ArrayList<ArrayList<Integer>> mutable2StatePosNums = new ArrayList<>();
    //For each state, a list of which flexible positions are mutable
    //these will be listed directly in Multistate.cfg under "STATEMUTRES0" etc.

    private ConfigFileParser cfp;

    public EWAKStarDoer (ConfigFileParser cfp) {

        this.cfp = cfp;

        //fill in all the settings
        //each state will have its own config file parser

        System.out.println();
        System.out.println("Performing EWAK*");
        System.out.println();

        numStates = 3; //for EWAKStar we have PL, P, and L

        //NUMMUTRES should be the number of mutable and flexible residues
        numTreeLevels = cfp.params.getInt("NUMMUTRES");

        // original: objFcn = new LME( cfp.params.getValue("OBJFCN"), 1 );
        objFcn = new EWAKLME( "1", 1 );

        SearchProblem[] stateSP = new SearchProblem[numStates];

        System.out.println();
        System.out.println("Preparing matrix and search problem for the bound complex.");
        System.out.println();

        // create a search problem for the bound complex, PL
        mutable2StatePosNums.add( stateMutablePos(0,cfp.params,numTreeLevels) );
        stateSP[0] = makeStateSearchProblem(0);

        System.out.println();
        System.out.println("Bound complex matrix ready.");
        System.out.println();

        //we can collect info on allowed sequences from any state (we choose state 0)
        //but we'll check that they all agree

        //needed for COMETs tree - we'll just use the highest total number of mutations
        int numMaxMut = cfp.params.getInt("NUMMAXMUT");

        //wt seq
        String[] wtSeq = null;
        wtSeq = parseWTSeq( cfp.params.getValue("WTSEQ"), numTreeLevels );

        numSeqsWanted = cfp.params.getInt("NUMSEQS");

        tree = new EWAKStarTree(numTreeLevels, objFcn,
                AATypeOptions, numMaxMut, wtSeq, numStates, stateSP,
                mutable2StatePosNums, false);
    }


    private String[] parseWTSeq(String seq, int numTreeLevels){
        //WT sequence (for purposes of capping # mutations) specified like ARG LYS...,
        StringTokenizer st = new StringTokenizer(seq);
        if(st.countTokens()!=numTreeLevels)
            throw new RuntimeException("ERROR: wrong number of residues in WT seq: "+seq);

        String wt[] = new String[numTreeLevels];

        for(int level=0; level<numTreeLevels; level++)
            wt[level] = st.nextToken().toUpperCase();

        return wt;
    }

    private ConfigFileParser makeStateConfig(int state) {

        //read state-specific configuration files and create a search problem object
        ConfigFileParser stateCFGP;

        if (cfp.params.getValue("STATEKSFILE"+ state).equalsIgnoreCase("DEFAULT")) {

            // state has no config file, inherit from the main CFP
            stateCFGP = new ConfigFileParser(cfp);

        } else {

            // otherwise, start with that config (and no defaults)
            stateCFGP = new ConfigFileParser();
            stateCFGP.params.addParamsFromFile(cfp.params.getFile("STATEKSFILE" + state));
        }

        //We expect input like
        //STATECFGFILES0 System0.cfg DEE0.cfg
        String stateConfigFiles = cfp.params.getValue("STATECFGFILES"+state);
        String stateSysFile = StringParsing.getToken(stateConfigFiles, 1);
        String stateDEEFile = StringParsing.getToken(stateConfigFiles, 2);

        PathRoot root = cfp.params.getRoot("STATECFGFILES" + state);
        stateCFGP.params.addParams(root, stateSysFile);
        stateCFGP.params.addParams(root, stateDEEFile);
        stateCFGP.loadData();
        return stateCFGP;
    }


    private SearchProblem makeStateSearchProblem(int state){

        ConfigFileParser stateCFGP = makeStateConfig(state);
        SearchProblem searchProb = stateCFGP.getSearchProblem();

        precomputeMatrices(state, stateCFGP, searchProb);

        handleAATypeOptions(state, stateCFGP);

        return searchProb;
    }


    private void handleAATypeOptions(int state, ConfigFileParser cfgP){
        //Given the config file parser for a state, make sure AATypeOptions
        //matches the allowed AA types for this state

        ArrayList<ArrayList<String>> stateAAOptions = cfgP.getAllowedAAs();

        Molecule wtMolec = new Strand.Builder(PDBIO.readFile(cfgP.params.getFile("PDBName"))).build().mol;
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


    public ArrayList<String> calcBestSequences(){

        System.out.println("Performing multistate A*");

        //how many sequences to enumerate

        long startAStarTime = System.currentTimeMillis();

        ArrayList<String> bestSequences = new ArrayList<>();

        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
            //this will find the best sequence and print it
            ScoredConf conf = tree.nextConf();
            if (conf == null) {
                //empty sequence...indicates no more sequence possibilities
                break;
            } else {
                bestSequences.add(tree.seqAsString(conf.getAssignments()));
            }
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));


        //DEBUG!!!  Checking by exhaustive search
        //exhaustiveMultistateSearch();
        return bestSequences;
    }


    void precomputeMatrices(int state, ConfigFileParser cfgP, SearchProblem sp){
        //Precompute energy and pruning matrices for a state
        //as if computing a GMEC

        GMECFinder gf = new GMECFinder();
        gf.init(cfgP, sp); //Using our search problem here instead of full GMEC calculation

        double ival = cfgP.params.getDouble("IVAL");

        gf.precomputeMatrices(ival);//no Ew...just looking at GMECs
        //and for any valid sequence the GMEC shouldn't be pruned (unless ival manually set lower)
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


    private double calcLowEnergy(int state, String[] AATypes){
        //Calculate a low energy conformation for the specified state for sequence AATypes
        ConfigFileParser stateCfp = makeStateConfig(state);

        //set up for particular AA types
        stateCfp.params.setValue("RUNNAME", "SEQ_"+System.currentTimeMillis());

        //also want to restrict to only this seq: addWT = false, however, include WT rotamers
        stateCfp.params.setValue("ADDWT", "false");
        stateCfp.params.setValue("ADDWTROTS", "true");

        String seq[] = new String[stateCfp.getFlexRes().size()];
        for(int mutPos=0; mutPos<numTreeLevels; mutPos++) {
            seq[mutable2StatePosNums.get(state).get(mutPos)] = AATypes[mutPos];
        }

        //matching format to ConfigFileParser.getAllowedAAs
        int posCount = 0;

        for(int str=0; str<10; str++){
            ArrayList<String> resAllowedRecords = stateCfp.params.searchParams("RESALLOWED"+str);
            int numRecordsInStrand = resAllowedRecords.size();

            //must go through residues in numerical order
            for(int recNum=0; recNum<numRecordsInStrand; recNum++){

                if(seq[posCount] != null){
                    String param = "RESALLOWED" + str + "_" + recNum;
                    stateCfp.params.setValue(param, seq[posCount]);
                }

                posCount++;
            }
        }

        EWAKStarLowEnergyFinder gf = new EWAKStarLowEnergyFinder();
        gf.init(stateCfp);
        return gf.calcSequences().getEnergy();
    }

}
