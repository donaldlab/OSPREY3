package edu.duke.cs.osprey.astar.ewakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.ewakstar.EWAKStarSequenceAnalyzer;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;

public class NewEWAKStarDoer {


    NewEWAKStarTree treePL;//The tree used for the EWAKStar PL search - based on NewCOMETSTree
    NewEWAKStarTree treeL; //The tree used for the EWAKStar L search
    int numSeqsWanted;//How many sequences to enumerate for bound complex

    Sequence wtSeq;

    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null; //AA types allowed at each mutable position

    ArrayList<Integer> mutablePosNums;

    SimpleConfSpace[] confSpaces;
    PrecomputedMatrices[] precompMats;
    ConfEnergyCalculator[] confECalcs;
    ArrayList<Sequence> bestSequences = new ArrayList<>();
    ArrayList<String> ewakstarSeqs = new ArrayList<>();
    Double unboundEw;
    Double boundEw;

    HashMap<Sequence,javafx.util.Pair<BigDecimal,Double>> ligandPFs = new HashMap<>();

    public NewEWAKStarDoer (SimpleConfSpace[] confSpaces, PrecomputedMatrices[] precompMats,
                            ArrayList<Integer> mutablePosNums, ArrayList<ArrayList<String>> AATypeOptions,
                             int numSeqsWanted, ConfEnergyCalculator[] confECalcs, Double unboundEw, Double boundEw) {

        //fill in all the settings
        //each state will have its own config file parser

        this.mutablePosNums = mutablePosNums;
        this.AATypeOptions = AATypeOptions;
        numTreeLevels = AATypeOptions.size();
        this.numSeqsWanted = numSeqsWanted;
        this.confSpaces = confSpaces;
        this.precompMats = precompMats;
        this.confECalcs = confECalcs;
        this.unboundEw = unboundEw;
        this.boundEw = boundEw;

        this.wtSeq = Sequence.makeWildType(confSpaces[1]);

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)

        treePL = new NewEWAKStarTree(numTreeLevels, AATypeOptions, wtSeq, confSpaces[0], precompMats[0],
                mutablePosNums, confECalcs[0]);

        treeL = new NewEWAKStarTree(numTreeLevels, AATypeOptions, wtSeq, confSpaces[1], precompMats[1],
                mutablePosNums, confECalcs[1]);


    }


    public ArrayList<String> calcBestSequences(){

        System.out.println("Performing EWAK*");
        ArrayList<Sequence> bestSeqs = extractSeqByLB();
        System.out.println(ewakstarSeqs);
//        ArrayList<Sequence> bestLseqs = extractLSeqsByLB();
//        ArrayList<EWAKStarSequenceAnalyzer.Analysis> unboundEnsembles = calcUnboundEnsembles();

//        for (EWAKStarSequenceAnalyzer.Analysis ens : unboundEnsembles){
//            ligandPFs.put(ens.sequence, calcUnboundPF(ens));
//        }
//
//        System.out.println(ligandPFs);
        return null;
    }

    public ArrayList<Sequence> extractSeqByLB(){

        long startAStarTime = System.currentTimeMillis();

        //get wild-type sequence for the unbound complex, L
        Sequence WT = Sequence.makeWildType(confSpaces[1]);
        String wtSeq = Sequence.makeWildTypeEWAKStar(WT);
        Boolean wtFound = false;
        Double wtEnergy = 0.0;

        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
            //this will find the best sequence and print it
            ScoredConf conf = treePL.nextConf();
            if (conf == null) {
                //empty sequence...indicates no more sequence possibilities
                break;
            } else {
                if(!wtFound) {
                    String curSeq = treePL.seqAsString(conf.getAssignments());
                    if (curSeq.equals(wtSeq)) {
                        wtFound = true;
                        //the WT score is the minimized energy (not the lower-bound)
                        wtEnergy = conf.getScore();
                        System.out.println("Found wild-type sequence in complex after " + seqNum + " sequences.");
                        System.out.println("Finding all sequences within " + boundEw + "of WT minimized energy: " + wtEnergy);
                    }
                    ewakstarSeqs.add(curSeq);
                    Sequence newSequence = Sequence.makeFromEWAKStar(treePL.seqAsString(conf.getAssignments()), WT, confSpaces[1]);
                    bestSequences.add(newSequence);
                }
                else {
                    Double curScore = conf.getScore();
                    if(curScore > wtEnergy+boundEw){
                        break;
                    }
                    else{
                        String curSeq = treePL.seqAsString(conf.getAssignments());
                        ewakstarSeqs.add(curSeq);
                        Sequence newSequence = Sequence.makeFromEWAKStar(treePL.seqAsString(conf.getAssignments()), WT, confSpaces[1]);
                        bestSequences.add(newSequence);

                    }
                }
            }
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        System.out.println("Number of sequences enumerated: "+Integer.toString(bestSequences.size()));

        if (!wtFound){
            System.out.println("WARNING: The wild type sequence was not found in your search!");
        }

        //returns the bound complex's sequences enumerated in order of increasing lower bound
        return bestSequences;

    }

//    public ArrayList<Sequence> extractLSeqsByLB(){
//
//        long startAStarTime = System.currentTimeMillis();
//
//        //get wild-type sequence for the unbound complex, L
//        Sequence WT = Sequence.makeWildType(confSpaces[1]);
//        System.out.println(WT);
//
//        while()
//        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
//            //this will find the best sequence and print it
//            ScoredConf conf = treePL.nextConf();
//            if (conf == null) {
//                //empty sequence...indicates no more sequence possibilities
//                break;
//            } else {
//                Double confScore = conf.getScore();
//                ewakstarSeqs.add(treePL.seqAsString(conf.getAssignments()));
//                Sequence newSequence = Sequence.makeFromEWAKStar(treePL.seqAsString(conf.getAssignments()), WT, confSpaces[1]);
//                bestSequences.add(newSequence);
//            }
//        }
//
//        long stopTime = System.currentTimeMillis();
//        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));
//
//        System.out.println("Number of sequences enumerated: "+Integer.toString(bestSequences.size()));
//
//        //returns the bound complex's sequences enumerated in order of increasing lower bound
//        return bestSequences;
//
//    }

    //are these sequences we want to actually calculate scores for? Let's check.
//    public ArrayList<Sequence> checkSequences(ArrayList<Sequence> bestSeqs){
//
//        KStar.ConfSearchFactory confSearchFactory = (confSpace, rcs) -> new ConfAStarTree.Builder(confSpace, rcs)
//                .setTraditional()
//                .build();
//
//        EWAKStarSequenceAnalyzer sa = new EWAKStarSequenceAnalyzer(confSpaces[1], confECalcs[1], precompMats[1], confSearchFactory);
//
//        for (Sequence s: bestSeqs){
//            EWAKStarSequenceAnalyzer.Analysis analysis = sa.analyze(s,eW);
//        }
//
//
//    }
    public ArrayList<EWAKStarSequenceAnalyzer.Analysis> calcUnboundEnsembles(){


        KStar.ConfSearchFactory confSearchFactory = (confSpace, rcs) -> new ConfAStarTree.Builder(confSpace, rcs)
                .setTraditional()
                .build();
        EWAKStarSequenceAnalyzer sa = new EWAKStarSequenceAnalyzer(confSpaces[1], confECalcs[1], precompMats[1], confSearchFactory);
        //collect all of the analysis for each sequence
        ArrayList<EWAKStarSequenceAnalyzer.Analysis> seqAnalysis = new ArrayList<>();
//        EWAKStarSequenceAnalyzer seqa = new EWAKStarSequenceAnalyzer(confSpaces[1], confECalcs[1], precompMats[1], confSearchFactory);

        for(Sequence s: bestSequences){
            EWAKStarSequenceAnalyzer.Analysis analysis = sa.analyze(s,unboundEw);
            seqAnalysis.add(analysis);
        }

//        for (EWAKStarSequenceAnalyzer.Analysis a : seqAnalysis){
//            System.out.println(a.toString());
//        }

        return seqAnalysis;
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
