package edu.duke.cs.osprey.astar.ewakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.ewakstar.EWAKStarSequenceAnalyzer;
import edu.duke.cs.osprey.gmec.ConfSearchFactory;
import edu.duke.cs.osprey.gmec.EWAKStarConfAnalyzer;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;
import edu.duke.cs.osprey.kstar.KStar;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tests.EWAKStar;

import java.lang.reflect.Array;
import java.math.BigDecimal;
import java.util.ArrayList;

public class NewEWAKStarDoer {


    NewEWAKStarTree tree;//The tree used for the EWAKStar search - based on NewCOMETSTree
    int numSeqsWanted;//How many sequences to enumerate for bound complex

    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null; //AA types allowed at each mutable position

    ArrayList<Integer> mutablePosNums;

    SimpleConfSpace[] confSpaces;
    PrecomputedMatrices[] precompMats;
    ConfEnergyCalculator[] confECalcs;
    ArrayList<Sequence> bestSequences = new ArrayList<Sequence>();
    double eW;

    ArrayList<BigDecimal> ligandPFs = new ArrayList<BigDecimal>();

    public NewEWAKStarDoer (SimpleConfSpace[] confSpaces, PrecomputedMatrices[] precompMats,
                            ArrayList<Integer> mutablePosNums, ArrayList<ArrayList<String>> AATypeOptions,
                            String wtSeq[], int numSeqsWanted, ConfEnergyCalculator[] confECalcs, double eW) {

        //fill in all the settings
        //each state will have its own config file parser

        this.mutablePosNums = mutablePosNums;
        this.AATypeOptions = AATypeOptions;
        numTreeLevels = AATypeOptions.size();
        this.numSeqsWanted = numSeqsWanted;
        this.confSpaces = confSpaces;
        this.precompMats = precompMats;
        this.confECalcs = confECalcs;
        this.eW = eW;

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)

        tree = new NewEWAKStarTree(numTreeLevels, AATypeOptions, wtSeq, confSpaces[0], precompMats[0],
                mutablePosNums, confECalcs[0]);


    }


    public ArrayList<String> calcBestSequences(){

        System.out.println("Performing EWAK*");
        ArrayList<Sequence> bestSeqs = extractSeqByLB();
        ArrayList<EWAKStarSequenceAnalyzer.Analysis> unboundEnsembles = calcUnboundEnsembles(bestSeqs);

        for (EWAKStarSequenceAnalyzer.Analysis ens : unboundEnsembles){
            ligandPFs.add(calcUnboundPF(ens));
        }

        System.out.println(ligandPFs);
        return null;
    }

    public ArrayList<Sequence> extractSeqByLB(){

        long startAStarTime = System.currentTimeMillis();

        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
            //this will find the best sequence and print it
            ScoredConf conf = tree.nextConf();
            if (conf == null) {
                //empty sequence...indicates no more sequence possibilities
                break;
            } else {
                Sequence mySeq = Sequence.makeFromAssignments(confSpaces[1],conf.getAssignments());
                bestSequences.add(mySeq);
            }
        }

        long stopTime = System.currentTimeMillis();
        System.out.println("Sequence enumeration time: "+((stopTime-startAStarTime)/(60.0*1000.0)));

        System.out.println("Number of sequences enumerated: "+Integer.toString(bestSequences.size()));

        //returns the bound complex's sequences enumerated in order of increasing lower bound
        return bestSequences;

    }
    public ArrayList<EWAKStarSequenceAnalyzer.Analysis> calcUnboundEnsembles(ArrayList<Sequence> bestSeqs){

        KStar.ConfSearchFactory confSearchFactory = (confSpace, rcs) -> new ConfAStarTree.Builder(confSpace, rcs)
                .setTraditional()
                .build();

        EWAKStarSequenceAnalyzer sa = new EWAKStarSequenceAnalyzer(confSpaces[1], confECalcs[1], precompMats[1], confSearchFactory);

        //collect all of the analysis for each sequence
        ArrayList<EWAKStarSequenceAnalyzer.Analysis> seqAnalysis = new ArrayList<EWAKStarSequenceAnalyzer.Analysis>();


        for(Sequence s: bestSeqs){
            EWAKStarSequenceAnalyzer.Analysis analysis = sa.analyze(s,eW);
            seqAnalysis.add(analysis);
        }

//        for (EWAKStarSequenceAnalyzer.Analysis a : seqAnalysis){
//            System.out.println(a.toString());
//        }

        return seqAnalysis;
    }

    public BigDecimal calcUnboundPF(EWAKStarSequenceAnalyzer.Analysis ens){

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);
        BigDecimal pf = BigDecimal.ZERO;

        for (int i=0; i<ens.ensemble.analyses.size(); i++) {
            BigDecimal bwEnergy = bc.calc(ens.ensemble.analyses.get(i).epmol.energy);
            pf = pf.add(bwEnergy);
        }

        return pf;

    }
}
