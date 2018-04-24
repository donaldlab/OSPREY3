package edu.duke.cs.osprey.astar.ewakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ewakstar.EWAKStarSequenceAnalyzer;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.gmec.PruningSettings;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import org.apache.commons.collections4.map.LinkedMap;

import java.io.File;
import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;

public class NewEWAKStarDoer {


    NewEWAKStarTree treePL;//The tree used for the EWAKStar PL search - based on NewCOMETSTree
    NewEWAKStarTree treeL; //The tree used for the EWAKStar L search
    int numSeqsWanted;//How many sequences to enumerate for bound complex - defaulted to 10000

    Sequence wtSeq;

    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null; //AA types allowed at each mutable position

    ArrayList<Integer> mutablePosNums;
    ArrayList<Integer> mutablePosNumsL;

    SimpleConfSpace confSpace;
    SimpleConfSpace confSpaceL;
    PrecomputedMatrices precompMat;
    ConfEnergyCalculator confECalc;

    LinkedMap<String,Double> ewakstarPLConfs = new LinkedMap<>();
    LinkedMap<String,Double> ewakstarLConfs = new LinkedMap<>();

    double unboundEw;
    double boundEw;
    String wtSeqEWAKStar;
    String startResL;
    String endResL;
    Molecule mol;
    String[] mutResNums;
    double Ival;
    PruningSettings pruningSettings;
    String LmatrixName;

    HashMap<Sequence,javafx.util.Pair<BigDecimal,Double>> ligandPFs = new HashMap<>();

    public NewEWAKStarDoer (SimpleConfSpace confSpace, SimpleConfSpace confSpaceL, PrecomputedMatrices precompMat,
                            ArrayList<Integer> mutablePosNums, ArrayList<Integer> mutablePosNumsL, ArrayList<ArrayList<String>> AATypeOptions,
                            int numSeqsWanted, ConfEnergyCalculator confECalc, double unboundEw, double boundEw,
                            String startResL, String endResL, Molecule mol, String[] mutResNums, double Ival,
                            PruningSettings pruningSettings, String LmatrixName) {

        //fill in all the settings
        //each state will have its own config file parser

        this.mutablePosNums = mutablePosNums;
        this.mutablePosNumsL = mutablePosNumsL;
        this.AATypeOptions = AATypeOptions;
        numTreeLevels = AATypeOptions.size();
        this.numSeqsWanted = numSeqsWanted;
        this.confSpace = confSpace;
        this.confSpaceL = confSpaceL;
        this.precompMat = precompMat;
        this.confECalc = confECalc;
        this.unboundEw = unboundEw;
        this.boundEw = boundEw;
        this.startResL = startResL;
        this.endResL = endResL;
        this.mol = mol;
        this.mutResNums = mutResNums;
        this.Ival = Ival;
        this.pruningSettings = pruningSettings;
        this.LmatrixName = LmatrixName;

        this.wtSeq = Sequence.makeWildType(confSpaceL);
        //get wild-type sequence for the unbound complex, L
        this.wtSeqEWAKStar = Sequence.makeWildTypeEWAKStar(wtSeq);

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)

        treePL = new NewEWAKStarTree(numTreeLevels, AATypeOptions, wtSeq, confSpace, precompMat,
                mutablePosNums, confECalc);

    }


    public ArrayList<String> calcBestSequences(){

        System.out.println("Performing EWAK*");
        ArrayList<Sequence> bestSeqs = extractSeqByLB();

        ArrayList<String> keyArray = new ArrayList<>();
        HashMap<String, ArrayList<String>> newAAOptionsDict = new HashMap<>();
        for(Sequence s: bestSeqs){
            String[] myAA = s.toString().split(" ");
            for (String set: myAA){
                String[] keyAndValue = set.split("=");
                String myKey = keyAndValue[0];
                String myVal = keyAndValue[1];
                if (!newAAOptionsDict.containsKey(myKey)){
                    keyArray.add(myKey);
                    ArrayList<String> aaTypes = new ArrayList<>();
                    aaTypes.add(myVal);
                    newAAOptionsDict.put(myKey,aaTypes);
                }
                else if(!newAAOptionsDict.get(myKey).contains(myVal)){
                    newAAOptionsDict.get(myKey).add(myVal);
                }
            }
        }

        ArrayList<ArrayList<String>> newAAOptions = new ArrayList<>();
        for(String myKey: keyArray){
            newAAOptions.add(newAAOptionsDict.get(myKey));
        }

        //make new strand and search space using the unbound ligand
        Strand strandL = new Strand.Builder(mol).setResidues(startResL, endResL).build();

        for(int mutPos=0; mutPos<newAAOptions.size(); mutPos++)
            strandL.flexibility.get(mutResNums[mutPos]).setLibraryRotamers(newAAOptions.get(mutPos)).setContinuous();

        SimpleConfSpace newConfSpace = new SimpleConfSpace.Builder().addStrand(strandL).build();

        ForcefieldParams ffparams = new ForcefieldParams();
        EnergyCalculator ecalc = new EnergyCalculator.Builder(newConfSpace, ffparams).build();

        ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(newConfSpace, ecalc);

        //use reference energies
        SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(newConfSpace,ecalc).build();
        confEcalcBuilder.setReferenceEnergies(eref);

        ConfEnergyCalculator newConfECalc = confEcalcBuilder.build();

        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(newConfECalc)
                .setCacheFile(new File(LmatrixName))
                .build()
                .calcEnergyMatrix();

        PrecomputedMatrices newPrecompMat = new PrecomputedMatrices(Ival, unboundEw, "L", emat,
                newConfSpace, ecalc, newConfECalc, new EPICSettings(), new LUTESettings(),
                pruningSettings);

        treeL = new NewEWAKStarTree(numTreeLevels, newAAOptions, wtSeq, newConfSpace, newPrecompMat,
                mutablePosNumsL, newConfECalc);


        ArrayList<Sequence> bestLseqs = extractLSeqsByLB(bestSeqs);
        System.out.println(bestLseqs);

//        ArrayList<EWAKStarSequenceAnalyzer.Analysis> unboundEnsembles = calcUnboundEnsembles();

//        for (EWAKStarSequenceAnalyzer.Analysis ens : unboundEnsembles){
//            ligandPFs.put(ens.sequence, calcUnboundPF(ens));
//        }
//
//        System.out.println(ligandPFs);
        return null;
    }

    public ArrayList<Sequence> extractSeqByLB(){

        ArrayList<Sequence> bestSequences = new ArrayList<>();

        long startAStarTime = System.currentTimeMillis();

        Boolean wtFound = false;
        double wtEnergy = 0.0;
        Boolean didEW = false;

        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
            //this will find the best sequence and print it
            ScoredConf conf = treePL.nextConf();
            if (conf == null) {
                //empty sequence...indicates no more sequence possibilities
                break;
            } else {
                String curSeq = treePL.seqAsString(conf.getAssignments());
                Sequence newSequence = Sequence.makeFromEWAKStar(treePL.seqAsString(conf.getAssignments()), wtSeq, confSpaceL);
                if(!wtFound) {
                    if (curSeq.equals(wtSeqEWAKStar)) {
                        wtFound = true;
                        //the WT score is the minimized energy (not the lower-bound)
                        wtEnergy = treePL.wtMinimizedEnergy;
                        System.out.println("Found wild-type sequence in complex after " + seqNum + " sequences.");
                        System.out.println("Finding all sequences within " + boundEw + " kcal of WT minimized energy: "
                                + wtEnergy+" or until "+numSeqsWanted+" sequences are enumerated.");

                    }

                    ewakstarPLConfs.put(curSeq, conf.getScore());
                    bestSequences.add(newSequence);
                }
                else {
                    double curScore = conf.getScore();
                    if(curScore > wtEnergy+boundEw){
                        didEW = true;
                        ewakstarPLConfs.put(curSeq, conf.getScore());
                        bestSequences.add(newSequence);
                        break;
                    }
                    else{
                        ewakstarPLConfs.put(curSeq, conf.getScore());
                        bestSequences.add(newSequence);

                    }
                }
            }
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        if (didEW)
            System.out.println("\nFound all within the energy window of "+boundEw+" kcal. Number of sequences enumerated: "+Integer.toString(bestSequences.size()));

        else
            System.out.println("\nEnergy window not completed - hit sequence enumeration limit of "+numSeqsWanted+" sequences.");

        System.out.println("\nLower bounds on energies range from "+ewakstarPLConfs.get(ewakstarPLConfs.firstKey())+" to "+ewakstarPLConfs.get(ewakstarPLConfs.lastKey())+"\n");

        if (!wtFound){
            System.out.println("\nWARNING: The wild type sequence was not found in your search!\n");
        }

        //returns the bound complex's sequences enumerated in order of increasing lower bound
        return bestSequences;

    }

    public ArrayList<Sequence> extractLSeqsByLB(ArrayList<Sequence> plSeqs){

        ArrayList<Sequence> bestSequences = new ArrayList<>();

        long startAStarTime = System.currentTimeMillis();

        Boolean wtFound = false;
        Boolean didEW = false;
        double wtEnergy = 0.0;
        int numSeqs = 0;
        while(ewakstarLConfs.size()<plSeqs.size()){
            //this will find the best sequence and print it
            ScoredConf conf = treeL.nextConf();
            if (conf == null) {
                //empty sequence...indicates no more sequence possibilities
                break;

            } else {
                String curSeq = treeL.seqAsString(conf.getAssignments());
                Sequence newSequence = Sequence.makeFromEWAKStar(treeL.seqAsString(conf.getAssignments()), wtSeq, confSpaceL);
                if(!wtFound) {
                    if (curSeq.equals(wtSeqEWAKStar)) {
                        wtFound = true;
                        //the WT score is the minimized energy (not the lower-bound)
                        wtEnergy = treeL.wtMinimizedEnergy;
                        System.out.println("Found wild-type sequence in unbound after " + numSeqs + " sequences.");
                        System.out.println("Finding all sequences within " + unboundEw + " kcal of WT minimized energy: "
                                + wtEnergy +" or until all PL complex sequences are enumerated.");
                        ewakstarLConfs.put(curSeq, conf.getScore());
                        bestSequences.add(newSequence);
                        numSeqs++;

                    }

                    else if (plSeqs.contains(newSequence)) {
                        numSeqs++;
                        ewakstarLConfs.put(curSeq, conf.getScore());
                        bestSequences.add(newSequence);
                    }
                }

                else {
                    double curScore = conf.getScore();
                    //we don't want to stop unless we've found all of the PL sequences.
                    if(numSeqs == plSeqs.size() && curScore > wtEnergy+unboundEw){
                        didEW = true;
                        if(plSeqs.contains(newSequence)){
                            ewakstarLConfs.put(curSeq, conf.getScore());
                            bestSequences.add(newSequence);
                            numSeqs++;
                        }
                        break;
                    }

                    else if(numSeqs == plSeqs.size()){
                        if(plSeqs.contains(newSequence)){
                            ewakstarLConfs.put(curSeq, conf.getScore());
                            bestSequences.add(newSequence);
                            numSeqs++;
                        }
                        break;
                    }

                    else if(curScore > wtEnergy+unboundEw){
                        didEW = true;
                        if(plSeqs.contains(newSequence)){
                            ewakstarLConfs.put(curSeq, conf.getScore());
                            bestSequences.add(newSequence);
                            numSeqs++;
                        }
                    }

                    else{
                        if (plSeqs.contains(newSequence)) {
                            numSeqs++;
                            ewakstarLConfs.put(curSeq, conf.getScore());
                            bestSequences.add(newSequence);
                        }
                    }
                }
            }
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        double ewDiff = Math.abs(wtEnergy - ewakstarLConfs.get(ewakstarLConfs.lastKey()));

        if (didEW)
            System.out.println("\nEnumerated the energy window of "+unboundEw+" kcal. Continued up to "+ewDiff+" kcal difference. Number of sequences enumerated: "+Integer.toString(bestSequences.size()));

        else
            System.out.println("\nEnergy window not completed - found all "+numSeqs+" sequences.");


        System.out.println("\nLower bounds on energies range from "+ewakstarLConfs.get(ewakstarLConfs.firstKey())+" to "+ewakstarLConfs.get(ewakstarLConfs.lastKey())+"\n");

        if (!wtFound){
            System.out.println("WARNING: The wild type sequence was not found in your search!");
        }

        //returns the unbound complex's sequences enumerated in order of increasing lower bound
        return bestSequences;

    }



    public EWAKStarSequenceAnalyzer.Analysis calcEnsemble(Sequence seq, SimpleConfSpace myConfSpace, ConfEnergyCalculator
            myConfECalc, PrecomputedMatrices myPrecompMat){


        KStar.ConfSearchFactory confSearchFactory = (confSpace, rcs) -> new ConfAStarTree.Builder(confSpace, rcs)
                .setTraditional()
                .build();
        EWAKStarSequenceAnalyzer sa = new EWAKStarSequenceAnalyzer(myConfSpace, myConfECalc, myPrecompMat, confSearchFactory);

        ArrayList<EWAKStarSequenceAnalyzer.Analysis> seqAnalysis = new ArrayList<>();
        EWAKStarSequenceAnalyzer.Analysis analysis = sa.analyze(seq ,unboundEw);

        return analysis;
    }



    public javafx.util.Pair<BigDecimal,Double> calcUnboundPF(EWAKStarSequenceAnalyzer.Analysis ens){

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        BigDecimal pf = BigDecimal.ZERO;

        for (int i=0; i<ens.ensemble.analyses.size(); i++) {
            BigDecimal bwEnergy = bc.calc(ens.ensemble.analyses.get(i).epmol.energy);
            pf = pf.add(bwEnergy);
        }

        javafx.util.Pair<BigDecimal, Double> pfs = new javafx.util.Pair<BigDecimal, Double>(pf, Math.log10(pf.doubleValue()));

        return pfs;

    }

}
