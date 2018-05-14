package edu.duke.cs.osprey.astar.ewakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.ewakstar.EWAKStar;
import edu.duke.cs.osprey.ewakstar.EWAKStarBBKStar;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.gmec.PruningSettings;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.tupexp.LUTESettings;
import org.apache.commons.collections4.map.LinkedMap;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;

/** author: lowegard **/

public class NewEWAKStarDoer {

    public static class ConfSpaces {
        public ForcefieldParams ffparams;
        public SimpleConfSpace protein;
        public SimpleConfSpace ligand;
        public SimpleConfSpace complex;
    }

    public static class Results {
        private EWAKStarBBKStar bbkstar;
        public List<EWAKStar.ScoredSequence> sequences;
    }


    private NewEWAKStarTree treePL;//The tree used for the EWAKStar PL search - based on NewCOMETSTree
    private NewEWAKStarTree treeL; //The tree used for the EWAKStar L search
    private NewEWAKStarTree treeP; //The tree used for the EWAKStar L search
    private int numSeqsWanted;//How many sequences to enumerate for bound complex - defaulted to 10000

    private Sequence fullWtSeq;
    private Sequence wtSeqL;
    private Sequence wtSeqP;

    private ArrayList<Integer> mutablePosNumsL;
    private ArrayList<Integer> mutablePosNumsP;

    private ConfSpaces confSpaces = new ConfSpaces();

    private LinkedMap<String,Double> ewakstarPLConfs = new LinkedMap<>();
    private LinkedMap<String,Double> ewakstarLConfs = new LinkedMap<>();
    private LinkedMap<String,Double> ewakstarPConfs = new LinkedMap<>();

    private double unboundEw;
    private double boundEw;
    private double epsilon;
    private double ewakstarEw;
    private int maxPFConfs;
    private double orderOfMag;
    private String wtSeqEWAKStar;
    private String startResL;
    private String endResL;
    private String startResP;
    private String endResP;
    private Molecule mol;
    private String[] resNumsPL;
    private String[] resNumsL;
    private String[] resNumsP;
    private double Ival;
    private String PLmatrixName;
    private EnergyCalculator minimizingEcalc;
    private EnergyMatrix ematPL;
    private EnergyMatrix ematL;
    private EnergyMatrix ematP;
    private int maxNumSeqs;

    private ConfEnergyCalculator confEnergyCalcL;
    private ConfEnergyCalculator confEnergyCalcP;
    private ConfEnergyCalculator confEnergyCalcPL;
    private ConfEnergyCalculator confRigidEnergyCalcL;
    private ConfEnergyCalculator confRigidEnergyCalcP;
    private ConfEnergyCalculator confRigidEnergyCalcPL;

    private boolean seqFilterOnly;

    private PruningSettings pruningSettings = new PruningSettings();

    public NewEWAKStarDoer (boolean seqFilterOnly, int maxNumSeqs, int maxPFConfs, double epsilon,
                            ConfEnergyCalculator confRigidEnergyCalcPL, ConfEnergyCalculator confEnergyCalcPL,
                            EnergyMatrix ematPL, EnergyCalculator minimizingEcalc, SimpleConfSpace confSpace,
                            SimpleConfSpace confSpaceL, SimpleConfSpace confSpaceP, Integer[] pos, Integer[] posL,
                            Integer[] posP, ArrayList<ArrayList<String>> AATypeOptions, int numSeqsWanted,
                            double orderOfMag, double unboundEw, double boundEw, double ewakstarEw, String startResL,
                            String endResL, String startResP, String endResP, Molecule mol, String[] resNumsPL,
                            String[] resNumsL, String[] resNumsP, double Ival, String PLmatrixName, ForcefieldParams ffparams) {

        //fill in all the settings
        //each state will have its own config file parser
        this.seqFilterOnly = seqFilterOnly;
        this.maxNumSeqs = maxNumSeqs;
        this.maxPFConfs = maxPFConfs;
        this.epsilon = epsilon;
        this.ematPL = ematPL;
        this.confEnergyCalcPL = confEnergyCalcPL;
        this.confRigidEnergyCalcPL = confRigidEnergyCalcPL;
        this.minimizingEcalc = minimizingEcalc;
        this.mutablePosNumsL = new ArrayList<>(Arrays.asList(posL));
        this.mutablePosNumsP = new ArrayList<>(Arrays.asList(posP));
        this.numSeqsWanted = numSeqsWanted;
        this.confSpaces.complex = confSpace;
        this.confSpaces.protein = confSpaceP;
        this.confSpaces.ligand = confSpaceL;
        this.confSpaces.ffparams = ffparams;
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
        this.PLmatrixName = PLmatrixName;

        this.fullWtSeq = Sequence.makeWildType(confSpace);
        this.wtSeqL = Sequence.makeWildType(confSpaceL);
        this.wtSeqP = Sequence.makeWildType(confSpaceP);
        //get wild-type sequence for the unbound complex, L
        this.wtSeqEWAKStar = Sequence.makeWildTypeEWAKStar(fullWtSeq);

        this.pruningSettings.typedep = true;

        PrecomputedMatrices newPrecompMat = new PrecomputedMatrices(Ival, boundEw, "PL", ematPL,
                confSpace, minimizingEcalc, confEnergyCalcPL, new EPICSettings(), new LUTESettings(),
                pruningSettings);

        treePL = new NewEWAKStarTree(AATypeOptions.size(), AATypeOptions, fullWtSeq, confSpace, newPrecompMat,
                new ArrayList<>(Arrays.asList(pos)), confEnergyCalcPL);

    }


    public ArrayList<Sequence> run(){

        System.out.println("Performing EWAK*");

        ArrayList<Sequence> bestPLseqs = extractPLSeqsByLB();

        ArrayList<Sequence> filteredSeqsL = filterSequences(bestPLseqs, "L");
        ArrayList<Sequence> filteredSeqsP = filterSequences(bestPLseqs, "P");

        ArrayList<ArrayList<String>> newAAOptions = updateAminoAcids(bestPLseqs);
        ArrayList<ArrayList<String>> newAAOptionsP = updateAminoAcids(filteredSeqsP);
        ArrayList<ArrayList<String>> newAAOptionsL = updateAminoAcids(filteredSeqsL);

        treeL = buildUnboundTree(newAAOptions, newAAOptionsL,"L");
        treeP = buildUnboundTree(newAAOptions, newAAOptionsP, "P");

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
        updatePLSettings();

        ArrayList<Sequence> bestFullseqs = new ArrayList<>();
        for(Sequence s:fullSeqs){
            newSeq = new Sequence(s, confSpaces.complex);
            bestFullseqs.add(newSeq);
        }
        this.fullWtSeq = new Sequence(fullWtSeq, confSpaces.complex);

        if (!seqFilterOnly) {
            runEWAKStarBBKStar(bestFullseqs);
        }

        if(seqFilterOnly) {
            writeSeqsToFile(bestFullseqs);
            return bestFullseqs;
        } else
            return null;
    }

    private void writeSeqsToFile(ArrayList<Sequence> seqs){

        File newFile = new File("ewakStar.filteredSeqs.txt");

        boolean append = true;
        boolean started = false;

        for (Sequence s: seqs){
            String curSeq = s.toString();
            if (!started) {
                started = true;
                append = false;
            }
            else{
                append = true;
            }
            try (FileWriter out = new FileWriter(newFile, append)) {
                out.write(curSeq);
                out.write("\n");
            } catch (IOException ex) {
                System.err.println("writing to file failed: " + newFile);
                ex.printStackTrace(System.err);
                System.err.println(curSeq);
            }
        }
    }

    private ArrayList<Sequence> filterSequences (ArrayList<Sequence> bestPLseqs, String type){

        boolean foundStart;
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
            foundStart = false;
            String newSeq = "";
            String[] seq = s.toString().split(" ");
            for (String str : seq) {
                String[] aaS = str.split("=");
                if (aaS[0].equals(resNumbers.get(0))) {
                    foundStart = true;
                }
                if(foundStart){
                    newSeq += aaS[1] + "_";
                    if(aaS[0].equals(resNumbers.get(resNumbers.size()-1))) {
                        Sequence curSeq = Sequence.makeFromEWAKStar(newSeq, wildtype, thisConfSpace);
                        if(!filteredSeqs.contains(curSeq)) {
                            filteredSeqs.add(curSeq);
                        }
                        break;
                    }
                }
            }
        }
        return filteredSeqs;
    }

    private void runEWAKStarBBKStar(ArrayList<Sequence> bestFullSeqs){

        // how should confs be ordered and searched?
        EWAKStar.ConfSearchFactory confSearchFactory = (emat, rcs) -> {
            return new ConfAStarTree.Builder(emat, rcs)
                    .setTraditional()
                    .build();
        };

        // run BBK*
        EWAKStar.Settings ewakstarSettings = new EWAKStar.Settings.Builder()
                .setEw(ewakstarEw)
                .setEpsilon(epsilon)
                .setMaxNumConfs(maxPFConfs)
                .setStabilityThreshold(null)
                .addScoreConsoleWriter()
                .setEnergyMatrixCachePattern(PLmatrixName)
                .addScoreFileWriter(new File("ewakStar.results.txt"))
                .build();

        EWAKStarBBKStar.Settings bbkstarSettings = new EWAKStarBBKStar.Settings.Builder()
                .setAllowedSeqs(bestFullSeqs).build();

        Results results = new Results();
        results.bbkstar = new EWAKStarBBKStar(maxNumSeqs, confSpaces.protein, confSpaces.ligand, confSpaces.complex,
                minimizingEcalc, confEnergyCalcPL, confEnergyCalcP, confEnergyCalcL, confRigidEnergyCalcPL,
                confRigidEnergyCalcP, confRigidEnergyCalcL, confSearchFactory, bbkstarSettings, ewakstarSettings,
                ematPL, ematP, ematL);
        results.sequences = results.bbkstar.run();

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


    private void updatePLSettings(){

        Strand strandP = confSpaces.protein.strands.get(0);
        Strand strandL = confSpaces.ligand.strands.get(0);

        this.confSpaces.complex = new SimpleConfSpace.Builder().addStrands(strandP, strandL).build();

        ForcefieldParams ffparams = new ForcefieldParams();
        EnergyCalculator ecalc = new EnergyCalculator.Builder(this.confSpaces.complex, ffparams).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(ecalc).setIsMinimizing(false).build();

        ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(this.confSpaces.complex, ecalc);
        ConfEnergyCalculator.Builder confRigidEcalcBuilder = new ConfEnergyCalculator.Builder(this.confSpaces.complex, rigidEcalc);

        //use reference energies
        SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(this.confSpaces.complex, ecalc).build();
        confEcalcBuilder.setReferenceEnergies(eref);
        confRigidEcalcBuilder.setReferenceEnergies(eref);

        this.confEnergyCalcPL = confEcalcBuilder.build();
        this.confRigidEnergyCalcPL = confRigidEcalcBuilder.build();

        // the pattern has a * right?
        if (PLmatrixName.indexOf('*') < 0) {
            throw new IllegalArgumentException("energyMatrixCachePattern (which is '" + PLmatrixName + "') has no wildcard character (which is *)");
        }

        String ematName = PLmatrixName.replace("*", "PL.emat");
        this.ematPL = new SimplerEnergyMatrixCalculator.Builder(this.confEnergyCalcPL)
                .setCacheFile(new File(ematName+".updated"))
                .build()
                .calcEnergyMatrix();

    }

    private NewEWAKStarTree buildUnboundTree(ArrayList<ArrayList<String>> newAAOptions, ArrayList<ArrayList<String>> newAAOptionsUB, String type) {

        Strand strandUnbound;
        ArrayList<Integer> mutablePos;
        String matrixName;
        Sequence wildType;

        switch (type) {
            case "P":
                strandUnbound = new Strand.Builder(mol).setResidues(startResP, endResP).build();
                mutablePos = mutablePosNumsP;
                matrixName = "ewakstar.P*";
                wildType = wtSeqP;
                break;
            case "L":
                strandUnbound = new Strand.Builder(mol).setResidues(startResL, endResL).build();
                mutablePos = mutablePosNumsL;
                matrixName = "ewakstar.L*";
                wildType = wtSeqL;
                break;
            default:
                throw new RuntimeException("Unrecognized state, cannot update amino acid types.");
        }

        for (Integer p : mutablePos)
            strandUnbound.flexibility.get(resNumsPL[p]).setLibraryRotamers(newAAOptions.get(p)).addWildTypeRotamers().setContinuous();

        for(int i=0; i<mutablePos.size();i++){
            mutablePos.set(i,i);
        }

        SimpleConfSpace newConfSpace = new SimpleConfSpace.Builder().addStrand(strandUnbound).build();

        if (type.equals("L"))
            this.confSpaces.ligand = newConfSpace;
        else if (type.equals("P"))
            this.confSpaces.protein = newConfSpace;

        ForcefieldParams ffparams = new ForcefieldParams();
        EnergyCalculator ecalc = new EnergyCalculator.Builder(newConfSpace, ffparams).build();
        EnergyCalculator rigidEcalc = new EnergyCalculator.SharedBuilder(ecalc).setIsMinimizing(false).build();

        ConfEnergyCalculator.Builder confEcalcBuilder = new ConfEnergyCalculator.Builder(newConfSpace, ecalc);
        ConfEnergyCalculator.Builder confRigidEcalcBuilder = new ConfEnergyCalculator.Builder(newConfSpace, rigidEcalc);

        //use reference energies
        SimpleReferenceEnergies eref = new SimpleReferenceEnergies.Builder(newConfSpace, ecalc).build();
        confEcalcBuilder.setReferenceEnergies(eref);
        confRigidEcalcBuilder.setReferenceEnergies(eref);

        ConfEnergyCalculator newConfECalc = confEcalcBuilder.build();
        ConfEnergyCalculator newRigidConfECalc = confRigidEcalcBuilder.build();

        if (type.equals("L")) {
            this.confEnergyCalcL = newConfECalc;
            this.confRigidEnergyCalcL = newRigidConfECalc;
        } else if (type.equals("P")) {
            this.confEnergyCalcP = newConfECalc;
            this.confRigidEnergyCalcP = newRigidConfECalc;
        }
        // the pattern has a * right?
        if (matrixName.indexOf('*') < 0) {
            throw new IllegalArgumentException("energyMatrixCachePattern (which is '" + matrixName + "') has no wildcard character (which is *)");
        }

        String ematName = matrixName.replace("*", ".emat");
        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(newConfECalc)
                .setCacheFile(new File(ematName))
                .build()
                .calcEnergyMatrix();

        if (type.equals("L"))
            this.ematL = emat;
        else if (type.equals("P"))
            this.ematP = emat;

        PrecomputedMatrices newPrecompMat = new PrecomputedMatrices(Ival, unboundEw, matrixName, emat,
                newConfSpace, ecalc, newConfECalc, new EPICSettings(), new LUTESettings(),
                pruningSettings);

        return new NewEWAKStarTree(newConfSpace.positions.size(), newAAOptionsUB, wildType, newConfSpace,
                newPrecompMat, mutablePos, newConfECalc);

    }

    public ArrayList<Sequence> extractPLSeqsByLB(){

        ArrayList<Sequence> bestSequences = new ArrayList<>();

        long startAStarTime = System.currentTimeMillis();

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        double wtPfUB = Double.POSITIVE_INFINITY;
        double pfUB;
        boolean didSeqMax = false;

        Boolean wtFound = false;
        double wtEnergy = 0.0;
        Boolean didEW = false;

        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
            //this will find the best sequence and print it
            ScoredConf conf = treePL.nextConf();
            if (conf == null) {
                //empty sequence...indicates no more sequence possibilities
                System.out.println("Ran out of sequences to enumerate.");
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
                        if(seqNum == numSeqsWanted)
                            didSeqMax = true;
                        break;
                    }
                    else if (pfUB > wtPfUB-orderOfMag){
                        ewakstarPLConfs.put(curSeq, pfUB);
                        bestSequences.add(newSequence);
                    }
                }
            }
            if(seqNum == numSeqsWanted)
                didSeqMax = true;
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        if (didEW)
            System.out.println("\nFound all within the energy window of "+boundEw+" kcal. \n"+Integer.toString(bestSequences.size())+" sequences enumerated with partition function upper bounds within "+orderOfMag+ " orders of magnitude of wild-type.");

        else if(didSeqMax)
            System.out.println("\nEnergy window not completed - hit sequence enumeration limit of "+numSeqsWanted+" sequences. Keeping sequences within "+orderOfMag+" orders of magnitude of the wild-type partition function, we have "+Integer.toString(bestSequences.size())+" sequences.");

        else
            System.out.println("\nFound "+Integer.toString(bestSequences.size())+" sequences within "+orderOfMag+" orders of magnitude of wild-type partition function.");

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
        boolean didSeqMax = false;


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
                System.out.println("Ran out of sequences to enumerate.");
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

                        if(seqsToUse.contains(newSequence) || newSequence.equals(curWT)) {
                            ewakstarConfs.put(curSeq, pfUB);
                            bestSequences.add(newSequence);
                            numSeqs++;
                        }

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

                        if(ewakstarConfs.size()==seqsToUse.size()){
                            didSeqMax = true;
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

                        if(ewakstarConfs.size()==seqsToUse.size()){
                            didSeqMax = true;
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
            if(ewakstarConfs.size()==seqsToUse.size()){
                didSeqMax = true;
            }
        }


        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        //double ewDiff = Math.abs(wtEnergy - (ewakstarLConfs.get(ewakstarLConfs.lastKey()).doubleValue());

        if (didEW)
            System.out.println("\nEnumerated the energy window of "+unboundEw+" kcal. \nEnumerated: "+Integer.toString(bestSequences.size())+" sequences within "+orderOfMag+" orders of magnitude of the wild-type unbound partition function.");

        else if(didSeqMax)
            System.out.println("\nEnergy window not completed - hit sequence enumeration limit of "+ seqsToUse.size()+" sequences. Keeping sequences within "+orderOfMag+" orders of magnitude of the wild-type partition function, we have "+Integer.toString(bestSequences.size())+" sequences.");

        else
            System.out.println("\nFound "+Integer.toString(bestSequences.size())+" sequences within "+orderOfMag+" orders of magnitude of wild-type partition function.");


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


    private ArrayList<Sequence> makeFullSeqs(ArrayList<Sequence> mutSeqsP, ArrayList<Sequence> mutSeqsL, ArrayList<Sequence> mutSeqsPL){

        ArrayList<String> fullSeq = new ArrayList<>(Arrays.asList(fullWtSeq.toString().split(" ")));
        ArrayList<Sequence> newFullSeqs = new ArrayList<>();
        if (resNumsPL[0].equals(resNumsP[0])) {
            for (Sequence p : mutSeqsP) {
                ArrayList<String> newFullSeq = new ArrayList<>(Arrays.asList(p.toString().split(" ")));
                for (Sequence l : mutSeqsL) {
                    newFullSeq.addAll(Arrays.asList(l.toString().split(" ")));
                    Sequence newSeq = Sequence.makeUnassigned(confSpaces.complex);
                    for (int i = 0; i < fullSeq.size(); i++) {
                        String[] newA = newFullSeq.get(i).split("=");
                        newSeq.set(newA[0], newA[1]);
                    }
                    if(mutSeqsPL.contains(newSeq))
                        newFullSeqs.add(newSeq);
                    newFullSeq = new ArrayList<>(Arrays.asList(p.toString().split(" ")));
                }
            }
        } else {
            for (Sequence l : mutSeqsL) {
                ArrayList<String> newFullSeq = new ArrayList<>(Arrays.asList(l.toString().split(" ")));
                for (Sequence p : mutSeqsP) {
                    newFullSeq.addAll(Arrays.asList(p.toString().split(" ")));
                    Sequence newSeq = Sequence.makeUnassigned(confSpaces.complex);
                    for (int i = 0; i < fullSeq.size(); i++) {
                        String[] newA = newFullSeq.get(i).split("=");
                        newSeq.set(newA[0], newA[1]);
                    }
                    if(mutSeqsPL.contains(newSeq))
                        newFullSeqs.add(newSeq);
                }
            }

        }

        return newFullSeqs;

    }
}
