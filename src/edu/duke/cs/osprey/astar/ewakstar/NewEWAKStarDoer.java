package edu.duke.cs.osprey.astar.ewakstar;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.gmec.PrecomputedMatrices;
import edu.duke.cs.osprey.confspace.ConfSearch.ScoredConf;

import java.util.ArrayList;

public class NewEWAKStarDoer {


    NewEWAKStarTree tree;//The tree used for the EWAKStar search - based on NewCOMETSTree
    int numSeqsWanted;//How many sequences to enumerate for bound complex

    int numTreeLevels;//number of mutable positions
    ArrayList<ArrayList<String>> AATypeOptions = null; //AA types allowed at each mutable position

    ArrayList<Integer> mutablePosNums;

    public NewEWAKStarDoer (SimpleConfSpace confSpace, PrecomputedMatrices precompMats,
                          ArrayList<Integer> mutablePosNums, ArrayList<ArrayList<String>> AATypeOptions,
                           String wtSeq[], int numSeqsWanted, ConfEnergyCalculator confECalc) {

        //fill in all the settings
        //each state will have its own config file parser

        this.mutablePosNums = mutablePosNums;
        this.AATypeOptions = AATypeOptions;
        numTreeLevels = AATypeOptions.size();
        this.numSeqsWanted = numSeqsWanted;

        //we can have a parameter numMaxMut to cap the number of deviations from the specified
        //wt seq (specified explicitly in case there is variation in wt between states...)

        tree = new NewEWAKStarTree(numTreeLevels, AATypeOptions, wtSeq, confSpace, precompMats,
                mutablePosNums, confECalc);
    }



    public ArrayList<String> calcBestSequences(){

        System.out.println("Performing EWAK*");


        //how many sequences to enumerate

        long startAStarTime = System.currentTimeMillis();

        ArrayList<String> bestSequences = new ArrayList<>();

        for(int seqNum=0; seqNum<numSeqsWanted; seqNum++){
            //this will find the best sequence and print it
            System.out.println("Sequence "+seqNum);
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

        return bestSequences;
    }



}
