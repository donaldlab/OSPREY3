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
import edu.duke.cs.osprey.ewakstar.EWAKStar;
import edu.duke.cs.osprey.ewakstar.EWAKStarBBKStar;
import edu.duke.cs.osprey.ewakstar.EWAKStarSequenceAnalyzer;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.gmec.PruningSettings;
import edu.duke.cs.osprey.kstar.BBKStar;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import org.apache.commons.collections4.map.LinkedMap;
import sun.java2d.pipe.SpanShapeRenderer;

import java.io.File;
import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.atomic.AtomicReference;

/** author: lowegard **/

public class NewEWAKStarDoer {

    public static class ConfSpaces {
        public ForcefieldParams ffparams;
        public SimpleConfSpace protein;
        public SimpleConfSpace ligand;
        public SimpleConfSpace complex;
    }

    public static class Results {
        public EWAKStarBBKStar bbkstar;
        public List<EWAKStar.ScoredSequence> sequences;
    }


    NewEWAKStarTree treePL;//The tree used for the EWAKStar PL search - based on NewCOMETSTree
    NewEWAKStarTree treeL; //The tree used for the EWAKStar L search
    NewEWAKStarTree treeP; //The tree used for the EWAKStar L search
    int numSeqsWanted;//How many sequences to enumerate for bound complex - defaulted to 10000

    Sequence fullWtSeq;
    Sequence wtSeqL;
    Sequence wtSeqP;

    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null; //AA types allowed at each mutable position

    ArrayList<Integer> mutablePosNums;
    ArrayList<Integer> mutablePosNumsL;
    ArrayList<Integer> mutablePosNumsP;

    ConfSpaces confSpaces = new ConfSpaces();

    PrecomputedMatrices precompMat;
    ConfEnergyCalculator confECalc;

    LinkedMap<String,Double> ewakstarPLConfs = new LinkedMap<>();
    LinkedMap<String,Double> ewakstarLConfs = new LinkedMap<>();
    LinkedMap<String,Double> ewakstarPConfs = new LinkedMap<>();

    double unboundEw;
    double boundEw;
    double ewakstarEw;
    double orderOfMag;
    String wtSeqEWAKStar;
    String startResL;
    String endResL;
    String startResP;
    String endResP;
    Molecule mol;
    String[] resNumsPL;
    String[] resNumsL;
    String[] resNumsP;
    double Ival;
    PruningSettings pruningSettings;
    String LmatrixName;
    String PmatrixName;

    HashMap<Sequence,javafx.util.Pair<BigDecimal,Double>> ligandPFs = new HashMap<>();

    public NewEWAKStarDoer (SimpleConfSpace confSpace, SimpleConfSpace confSpaceL, SimpleConfSpace confSpaceP, PrecomputedMatrices precompMat,
                            ArrayList<Integer> mutablePosNums, ArrayList<Integer> mutablePosNumsL, ArrayList<Integer> mutablePosNumsP,
                            ArrayList<ArrayList<String>> AATypeOptions, int numSeqsWanted, ConfEnergyCalculator confECalc,
                            double orderOfMag, double unboundEw, double boundEw, double ewakstarEw, String startResL, String endResL, String startResP,
                            String endResP, Molecule mol, String[] resNumsPL, String[] resNumsL, String[] resNumsP, double Ival, PruningSettings pruningSettings, String LmatrixName, String PmatrixName,
                            ForcefieldParams ffparams) {

        //fill in all the settings
        //each state will have its own config file parser

        this.mutablePosNums = mutablePosNums;
        this.mutablePosNumsL = mutablePosNumsL;
        this.mutablePosNumsP = mutablePosNumsP;
        this.AATypeOptions = AATypeOptions;
        numTreeLevels = AATypeOptions.size();
        this.numSeqsWanted = numSeqsWanted;
        this.confSpaces.complex = confSpace;
        this.confSpaces.protein = confSpaceP;
        this.confSpaces.ligand = confSpaceL;
        this.confSpaces.ffparams = ffparams;
        this.precompMat = precompMat;
        this.confECalc = confECalc;
        this.orderOfMag = orderOfMag;
        this.unboundEw = unboundEw;
        this.boundEw = boundEw;
        this.ewakstarEw = ewakstarEw;
        this.startResL = startResL;
        this.endResL = endResL;
        this.startResP = startResP;
        this.endResP = endResP;
        this.mol = mol;
        this.resNumsPL = resNumsPL;
        this.resNumsL = resNumsL;
        this.resNumsP = resNumsP;
        this.Ival = Ival;
        this.pruningSettings = pruningSettings;
        this.LmatrixName = LmatrixName;
        this.PmatrixName = PmatrixName;

        this.fullWtSeq = Sequence.makeWildType(confSpace);
        this.wtSeqL = Sequence.makeWildType(confSpaceL);
        this.wtSeqP = Sequence.makeWildType(confSpaceP);
        //get wild-type sequence for the unbound complex, L
        this.wtSeqEWAKStar = Sequence.makeWildTypeEWAKStar(fullWtSeq);

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)

        treePL = new NewEWAKStarTree(numTreeLevels, AATypeOptions, fullWtSeq, confSpace, precompMat,
                mutablePosNums, confECalc);

    }


    public ArrayList<String> calcBestSequences(){

        System.out.println("Performing EWAK*");

        ArrayList<Sequence> bestPLseqs = extractPLSeqsByLB();

        ArrayList<Sequence> filteredSeqsL = filterSequences(bestPLseqs, "L");
        ArrayList<Sequence> filteredSeqsP = filterSequences(bestPLseqs, "P");

        ArrayList<ArrayList<String>> newAAOptionsL = updateAminoAcids(filteredSeqsL);
        ArrayList<ArrayList<String>> newAAOptionsP = updateAminoAcids(filteredSeqsP);

        treeL = buildUnboundTree(newAAOptionsL, "L");
        treeP = buildUnboundTree(newAAOptionsP, "P");

        Sequence newSeq;
        ArrayList<Sequence> bestLseqs = new ArrayList<>();
        for(Sequence s:filteredSeqsL){
            newSeq = new Sequence(s, confSpaces.ligand);
            bestLseqs.add(newSeq);
        }
        this.wtSeqL = new Sequence(wtSeqL, confSpaces.ligand);

        ArrayList<Sequence> bestPseqs = new ArrayList<>();
        for(Sequence s:filteredSeqsP){
            newSeq = new Sequence(s, confSpaces.protein);
            bestPseqs.add(newSeq);
        }
        this.wtSeqP = new Sequence(wtSeqP, confSpaces.protein);

        ArrayList<Sequence> newLseqs = extractSeqsByLB(bestLseqs, "L");
        ArrayList<Sequence> newPseqs = extractSeqsByLB(bestPseqs, "P");

        ArrayList<Sequence> fullSeqs = makeFullSeqs(newPseqs, newLseqs, bestPLseqs);

        System.out.println(fullSeqs.size());

        Results ewakstarResults = runEWAKStarBBKStar(fullSeqs);

        System.out.println(ewakstarResults.toString());

        return null;
    }

    private ArrayList<Sequence> filterSequences (ArrayList<Sequence> bestPLseqs, String type){

        boolean foundStart = false;
        boolean foundEnd = false;
        ArrayList<Sequence> filteredSeqs = new ArrayList<>();
        ArrayList<String> resNumbers;
        Sequence wildtype;
        SimpleConfSpace thisConfSpace;

        if (type.equals("L")){
            resNumbers = new ArrayList<>(Arrays.asList(resNumsL));
            wildtype = wtSeqL;
            thisConfSpace = confSpaces.ligand;
        } else{
            resNumbers = new ArrayList<>(Arrays.asList(resNumsP));
            wildtype = wtSeqP;
            thisConfSpace = confSpaces.protein;
        }

        for (Sequence s: bestPLseqs) {
            foundEnd = false;
            foundStart = false;
            String newSeq = "";
            String[] seq = s.toString().split(" ");
            for (String str : seq) {
                String[] aaS = str.split("=");
                if (aaS[0].equals(resNumbers.get(0))) {
                    foundStart = true;
                    newSeq += aaS[1] + "_";
                } else if (foundStart && !foundEnd) {
                    if (aaS[0].equals(resNumbers.get(resNumbers.size()-1))) {
                        foundEnd = true;
                        newSeq += aaS[1] + "_";
                    } else {
                        newSeq += aaS[1] + "_";
                    }
                    if(foundEnd && aaS[0].equals(resNumbers.get(resNumbers.size()-1))){
                        filteredSeqs.add(Sequence.makeFromEWAKStar(newSeq, wildtype, thisConfSpace));
                    }
                } else if (foundEnd && !filteredSeqs.contains(Sequence.makeFromEWAKStar(newSeq, wildtype, thisConfSpace))){
                    filteredSeqs.add(Sequence.makeFromEWAKStar(newSeq, wildtype, thisConfSpace));
                }
            }
        }
        return filteredSeqs;
    }

    private Results runEWAKStarBBKStar(ArrayList<Sequence> inputSeqs){

        AtomicReference<Results> resultsRef = new AtomicReference<>(null);

        Parallelism parallelism = Parallelism.makeCpu(4);
        //Parallelism parallelism = Parallelism.make(4, 1, 8);

        // how should we compute energies of molecules?
        new EnergyCalculator.Builder(confSpaces.complex, confSpaces.ffparams)
                .setParallelism(parallelism)
                .use((minimizingEcalc) -> {

                    // how should we define energies of conformations?
                    EWAKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
                        return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, ecalcArg)
                                        .build()
                                        .calcReferenceEnergies()
                                ).build();
                    };

                    // how should confs be ordered and searched?
                    EWAKStar.ConfSearchFactory confSearchFactory = (emat, rcs) -> {
                        return new ConfAStarTree.Builder(emat, rcs)
                                .setTraditional()
                                .build();
                    };

                    // make a rigid energy calculator too
                    EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(minimizingEcalc)
                            .setIsMinimizing(false)
                            .build();

                    // run BBK*
                    EWAKStar.Settings kstarSettings = new EWAKStar.Settings.Builder()
                            .setEw(ewakstarEw)
                            .setStabilityThreshold(null)
                            .addScoreConsoleWriter()
                            .build();
                    EWAKStarBBKStar.Settings bbkstarSettings = new EWAKStarBBKStar.Settings.Builder()
                            .setNumConfsPerBatch(8)
                            .build();
                    Results results = new Results();
                    results.bbkstar = new EWAKStarBBKStar(confSpaces.protein, confSpaces.ligand, confSpaces.complex, rigidEcalc, minimizingEcalc, confEcalcFactory, confSearchFactory, kstarSettings, bbkstarSettings);
                    results.sequences = results.bbkstar.run(inputSeqs);

                    // pass back the ref
                    resultsRef.set(results);
                });

        return resultsRef.get();
    }

    private ArrayList<ArrayList<String>> updateAminoAcids(ArrayList<Sequence> bestSeqs){

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

        return newAAOptions;
    }

    private NewEWAKStarTree buildUnboundTree(ArrayList<ArrayList<String>> newAAOptions, String type){

        Strand strandUnbound;
        ArrayList<Integer> mutablePos;
        String[] residueNums;
        String matrixName;
        Sequence wildType;

        switch(type){
            case "P":
                strandUnbound = new Strand.Builder(mol).setResidues(startResP, endResP).build();
                mutablePos = mutablePosNumsP;
                residueNums = resNumsP;
                matrixName = PmatrixName;
                wildType = wtSeqP;
                break;
            case "L":
                strandUnbound = new Strand.Builder(mol).setResidues(startResL, endResL).build();
                mutablePos = mutablePosNumsL;
                residueNums = resNumsL;
                matrixName = LmatrixName;
                wildType = wtSeqL;
                break;
            default:
                throw new RuntimeException("Unrecognized state, cannot update amino acid types.");
        }


        for(Integer p: mutablePos)
            strandUnbound.flexibility.get(residueNums[p]).setLibraryRotamers(newAAOptions.get(p)).setContinuous();

        SimpleConfSpace newConfSpace = new SimpleConfSpace.Builder().addStrand(strandUnbound).build();

        if (type.equals("L"))
            this.confSpaces.ligand = newConfSpace;
        else
            this.confSpaces.protein = newConfSpace;

        ForcefieldParams ffparams = new ForcefieldParams();
        EnergyCalculator ecalc = new EnergyCalculator.Builder(newConfSpace, ffparams).build();

        ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(newConfSpace, ecalc);

        //use reference energies
        SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(newConfSpace,ecalc).build();
        confEcalcBuilder.setReferenceEnergies(eref);

        ConfEnergyCalculator newConfECalc = confEcalcBuilder.build();

        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(newConfECalc)
                .setCacheFile(new File(matrixName))
                .build()
                .calcEnergyMatrix();

        PrecomputedMatrices newPrecompMat = new PrecomputedMatrices(Ival, unboundEw, "L", emat,
                newConfSpace, ecalc, newConfECalc, new EPICSettings(), new LUTESettings(),
                pruningSettings);

        return new NewEWAKStarTree(newConfSpace.positions.size(), newAAOptions, wildType, newConfSpace, newPrecompMat,
                mutablePos, newConfECalc);

    }

    public ArrayList<Sequence> extractPLSeqsByLB(){

        ArrayList<Sequence> bestSequences = new ArrayList<>();

        long startAStarTime = System.currentTimeMillis();

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        double wtPfUB = Double.POSITIVE_INFINITY;
        double pfUB;

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
                Sequence newSequence = Sequence.makeFromEWAKStar(treePL.seqAsString(conf.getAssignments()), fullWtSeq, confSpaces.complex);
                pfUB = Math.log10(bc.calc(conf.getScore()).multiply(new BigDecimal (getNumConfsForSeq(newSequence))).doubleValue());
                if(!wtFound) {
                    if (curSeq.equals(wtSeqEWAKStar)) {
                        wtFound = true;
                        wtEnergy = treePL.wtMinimizedEnergy;
                        wtPfUB = pfUB;
                        System.out.println("Found wild-type sequence in complex after " + seqNum + " sequences.");
                        System.out.println("Finding all sequences within " + boundEw + " kcal of WT minimized energy: "
                                + wtEnergy+" or until "+numSeqsWanted+" sequences are enumerated.");

                    }
                    ewakstarPLConfs.put(curSeq, pfUB);
                    bestSequences.add(newSequence);
                }
                //only keep sequences that will provably have a partition function within orderOfMag of wild-type
                //also keeps any sequence where pfUB is infinity
                else{
                    if(conf.getScore() > wtEnergy+boundEw){
                        didEW = true;
                        if (pfUB > wtPfUB-orderOfMag) {
                            ewakstarPLConfs.put(curSeq, pfUB);
                            bestSequences.add(newSequence);
                        }
                        break;
                    }
                    else if (pfUB > wtPfUB-orderOfMag){
                        ewakstarPLConfs.put(curSeq, pfUB);
                        bestSequences.add(newSequence);
                    }
                }
            }
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        if (didEW)
            System.out.println("\nFound all within the energy window of "+boundEw+" kcal. \n"+Integer.toString(bestSequences.size())+" sequences enumerated with partition function upper bounds within "+orderOfMag+ " orders of magnitude of wild-type.");

        else
            System.out.println("\nEnergy window not completed - hit sequence enumeration limit of "+numSeqsWanted+" sequences. Keeping sequences within "+orderOfMag+" orders of magnitude of the wild-type partition function, we have "+Integer.toString(bestSequences.size())+" sequences.");

        System.out.println("\nUpper bounds on partition functions range from "+ewakstarPLConfs.get(ewakstarPLConfs.firstKey())+" to "+ewakstarPLConfs.get(ewakstarPLConfs.lastKey())+"\nThe upper bound on the partition function for the wild-type is "+wtPfUB);

        if (!wtFound){
            System.out.println("\nWARNING: The wild type sequence was not found in your search!\n");
        }

        //returns the bound complex's sequences enumerated in order of increasing lower bound
        return bestSequences;

    }

    public ArrayList<Sequence> extractSeqsByLB(ArrayList<Sequence> seqsToUse, String type){

        ArrayList<Sequence> bestSequences = new ArrayList<>();
        SimpleConfSpace curConfSpace;
        NewEWAKStarTree curTree;
        LinkedMap<String, Double> ewakstarConfs;
        Sequence curWT;
        String startRes;
        String endRes;


        if (type.equals("L")) {
            curConfSpace = confSpaces.ligand;
            curTree = treeL;
            ewakstarConfs = ewakstarLConfs;
            curWT = wtSeqL;
        }
        else{
            curConfSpace = confSpaces.protein;
            curTree = treeP;
            ewakstarConfs = ewakstarPConfs;
            curWT = wtSeqP;
        }

        long startAStarTime = System.currentTimeMillis();

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        Boolean wtFound = false;
        Boolean didEW = false;
        double wtEnergy = 0.0;
        double wtPfUB = Double.POSITIVE_INFINITY;
        int numSeqs = 0;
        double pfUB;

        while(ewakstarConfs.size()<seqsToUse.size()){
            //this will find the best sequence and print it
            ScoredConf conf = curTree.nextConf();
            if (conf == null) {
                //empty sequence...indicates no more sequence possibilities
                break;

            } else {
                String curSeq = curTree.seqAsString(conf.getAssignments());
                Sequence newSequence = Sequence.makeFromEWAKStar(curTree.seqAsString(conf.getAssignments()), curWT, curConfSpace);
                BigDecimal numConfs = new BigDecimal (getNumConfsForSeq(newSequence));
                pfUB = Math.log10(bc.calc(conf.getScore()).multiply(numConfs).doubleValue());
                if(!wtFound) {
                    if (newSequence.toString().equals(curWT.toString())) {
                        wtFound = true;
                        wtEnergy = curTree.wtMinimizedEnergy;
                        wtPfUB = pfUB;
                        System.out.println("Found wild-type sequence in unbound after " + numSeqs + " sequences.");
                        System.out.println("Finding all sequences within " + unboundEw + " kcal of WT minimized energy: "
                                + wtEnergy +" or until all PL complex sequences are enumerated.");

                        ewakstarConfs.put(curSeq, pfUB);
                        bestSequences.add(newSequence);
                        numSeqs++;

                    }

                    else if (seqsToUse.contains(newSequence)) {
                        numSeqs++;
                        ewakstarConfs.put(curSeq, pfUB);
                        bestSequences.add(newSequence);
                    }
                }

                else {
                    //we want to stop once we're outside of the eW
                    if(conf.getScore() > wtEnergy+unboundEw){
                        didEW = true;
                        if(seqsToUse.contains(newSequence) && pfUB > wtPfUB-orderOfMag){
                            ewakstarConfs.put(curSeq, pfUB);
                            bestSequences.add(newSequence);
                            numSeqs++;
                        }
                        break;
                    }

                    //we also want to stop if we've found all of the PL sequences even if we're not outside the window yet.
                    else if(numSeqs == seqsToUse.size()){
                        if(seqsToUse.contains(newSequence) && pfUB > wtPfUB-orderOfMag){
                            ewakstarConfs.put(curSeq, pfUB);
                            bestSequences.add(newSequence);
                            numSeqs++;
                        }
                        break;
                    }

                    else{
                        if (seqsToUse.contains(newSequence) && pfUB > wtPfUB-orderOfMag) {
                            numSeqs++;
                            ewakstarConfs.put(curSeq, pfUB);
                            bestSequences.add(newSequence);
                        }
                    }
                }
            }
        }


        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        //double ewDiff = Math.abs(wtEnergy - (ewakstarLConfs.get(ewakstarLConfs.lastKey()).doubleValue());

        if (didEW)
            System.out.println("\nEnumerated the energy window of "+unboundEw+" kcal. \nEnumerated: "+Integer.toString(bestSequences.size())+" sequences within "+orderOfMag+" orders of magnitude of the wild-type unbound partition function.");

        else
            System.out.println("\nEnergy window not completed - found all "+numSeqs+" sequences within "+orderOfMag+" orders of magnitude of the wild-type unbound partition function.");


        System.out.println("\nUpper bounds on partition functions range from "+ewakstarConfs.get(ewakstarConfs.firstKey())+" to "+ewakstarConfs.get(ewakstarConfs.lastKey())+"\nThe upper bound on the partition function for the wild-type sequence is "+wtPfUB);

        if (!wtFound){
            System.out.println("WARNING: The wild type sequence was not found in your search!");
        }

        //returns the unbound complex's sequences enumerated in order of increasing lower bound
        return bestSequences;

    }

    public BigInteger getNumConfsForSeq(Sequence seq){
        List<SimpleConfSpace.Position> posConfs = seq.confSpace.positions;
        BigInteger numConfs = BigInteger.ONE;

        for (SimpleConfSpace.Position p : posConfs){
            String curRes = seq.get(p);
            BigInteger count = BigInteger.ZERO;
            for (SimpleConfSpace.ResidueConf rc : p.resConfs){
                if (rc.template.name.equals(curRes)){
                    count = count.add(BigInteger.ONE);
                }
            }
            if (count.compareTo(BigInteger.ZERO) == 0){
                System.out.println("OH NO! IT IS ZERO!");
            }
            numConfs = numConfs.multiply(count);
        }
        return numConfs;
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

    private ArrayList<Sequence> makeFullSeqs(ArrayList<Sequence> mutSeqsP, ArrayList<Sequence> mutSeqsL, ArrayList<Sequence> mutSeqsPL){

        ArrayList<String> fullSeq = new ArrayList<>(Arrays.asList(fullWtSeq.toString().split(" ")));
        ArrayList<Sequence> newFullSeqs = new ArrayList<>();
        if (resNumsPL[0].equals(resNumsP[0])) {
            for (Sequence p : mutSeqsP) {
                ArrayList<String> newFullSeq = new ArrayList<>(Arrays.asList(p.toString().split(" ")));
                for (Sequence l : mutSeqsL) {
                    newFullSeq.addAll(Arrays.asList(l.toString().split(" ")));
                    int j = 0;
                    Sequence newSeq = Sequence.makeUnassigned(confSpaces.complex);
                    for (int i = 0; i < fullSeq.size(); i++) {
                        String[] wtA = fullSeq.get(i).split("=");
                        String[] newA = newFullSeq.get(j).split("=");
                        if (wtA[0].equals(newA[0])) {
                            newSeq.set(newA[0], newA[1]);
                        } else {
                            newSeq.set(wtA[0], wtA[1]);
                        }
                    }
                    if(mutSeqsPL.contains(newSeq))
                        newFullSeqs.add(newSeq);
                }
            }
        } else {
            for (Sequence l : mutSeqsL) {
                List<String> newFullSeq = Arrays.asList(l.toString().split(" "));
                for (Sequence p : mutSeqsP) {
                    newFullSeq.addAll(Arrays.asList(p.toString().split(" ")));
                    int j = 0;
                    Sequence newSeq = Sequence.makeUnassigned(confSpaces.complex);
                    for (int i = 0; i < fullSeq.size(); i++) {
                        String[] wtA = fullSeq.get(i).split("=");
                        String[] newA = newFullSeq.get(j).split("=");
                        if (wtA[0].equals(newA[0])) {
                            newSeq.set(newA[0], newA[1]);
                        } else {
                            newSeq.set(wtA[0], wtA[1]);
                        }
                    }
                    if(mutSeqsPL.contains(newSeq))
                        newFullSeqs.add(newSeq);
                }
            }

        }

        return newFullSeqs;

    }
}
