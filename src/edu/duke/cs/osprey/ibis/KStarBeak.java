/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.confspace.TupleEnumerator;
import edu.duke.cs.osprey.control.ConfigFileParser;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.epic.EPICMatrix;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tools.ObjectIO;
import edu.duke.cs.osprey.tupexp.TupExpChooser;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

/**
 *
 * Special beak for K*, relies on tuple expander machinery
 * rather than dealing with repeated matrix entries
 * thus easier to do higher-order terms wrt sequence
 *
 *
 * @author mhall44
 */
public class KStarBeak {

    EnergyMatrix luteMatrix;
    PruningMatrix lutePruningMatrix;//Pruning matrix to go with the luteMatrix
    int numPos;

    ArrayList<ArrayList<String> > AATypes;//what amino acid types are allowed at each design position; defines indexing in seqEMat, etc
    ArrayList<ArrayList<ArrayList<Integer> > > RCsForAATypes;//list of RCs for each AA type at each position

    //pruning and energy matrix defining a "LUTE"/unpruned cluster expansion in seq space
    EnergyMatrix seqEMat = null;
    PruningMatrix seqPruneMat = null;


    SearchProblem sp = null;//DEBUG!!  used in checking EPIC


    public KStarBeak(EnergyMatrix luteMatrix, PruningMatrix lutePruningMatrix, ConfSpace confSpace) {
        this.luteMatrix = luteMatrix;
        this.lutePruningMatrix = lutePruningMatrix;
        numPos = luteMatrix.getNumPos();

        //let's generate the AA type info from a ConfSpace
        AATypes = new ArrayList<>();
        RCsForAATypes = new ArrayList<>();
        for(int pos=0; pos<numPos; pos++){
            HashMap<String,ArrayList<Integer> > rcLists = new HashMap<>();
            for(RC rc : confSpace.posFlex.get(pos).RCs){
                if(!rcLists.containsKey(rc.AAType))
                    rcLists.put(rc.AAType, new ArrayList<>());
                rcLists.get(rc.AAType).add(rc.RCIndex);
            }
            ArrayList<String> posAATypes = new ArrayList<>();
            ArrayList<ArrayList<Integer> > posRCsForAA = new ArrayList<>();
            for(String aa : rcLists.keySet()){
                posAATypes.add(aa);
                posRCsForAA.add(rcLists.get(aa));
            }
            AATypes.add(posAATypes);
            RCsForAATypes.add(posRCsForAA);
        }
    }



    public void makeUnprunedClusterExpansion(){
        //convert the LUTE matrix to a special energy matrix (indexed by AA's not RCs!)
        //that represents, for unpruned sequences

        //expansion at current stage of integration, defined as energy and pruning matrices
        PruningMatrix curPruneMat = lutePruningMatrix;
        EnergyMatrix curEMat = luteMatrix;

        double totResidsq = 0;

        for(int pos=0; pos<numPos; pos++){
            //Figure out what RCs are available for each amino-acid type

            //Figure out what AAs, RCs, and tuples (of each or mixed) are pruned
            PruningMatrix integPruneMat = integratedPruningMatrix(curPruneMat, pos);
            //Estimate integrated interactions for purposes of tuple selection
            EnergyMatrix estIntegEMat = estimateIntegEMat(curEMat, pos);

            IBISKStarTupleExpander expander = new IBISKStarTupleExpander(pos, curPruneMat, integPruneMat, curEMat, RCsForAATypes.get(pos));
            TupleEnumerator tupEnum = new TupleEnumerator(integPruneMat, estIntegEMat, numPos);
            TupExpChooser chooser = new TupExpChooser(expander, tupEnum);

            double resid = chooser.calcPairwiseExpansion();
            System.out.println("IBIS K* Residual for position "+pos+": "+resid);
            //DEBUG!! If high resid will probably want to iterate to higher order
            //may focus on triples of integrated AAs?
            //may want to use the expansion at each stage as est integ e mat for the next...
            //YEAH may not use the est integ e mat from here...null for pairawise, pairwise lute exp is good est mat for triples

            //NOTE even if cannot get good expansion for general K*,
            //may be able to get a (1-mutation-wise x LUTE for ligand??) matrix describing point mutations well
            //for use in priarm
            //priarm is hard in having multiple layers but easy in having a small protein mut space.
            //Exploit that.

            totResidsq += resid;

            //OK we have a good expansion now so let's update it to go to the next integration step
            curPruneMat = integPruneMat;
            curEMat = expander.getEnergyMatrix();
        }

        System.out.println("IBIS tot residsq: "+totResidsq);

        //store final integrated expansion
        seqEMat = curEMat;
        seqPruneMat = curPruneMat;
    }


    PruningMatrix integratedPruningMatrix(PruningMatrix p, int integPos){
        //convert a (partially integrated?) pruning matrix p to one at which position pos
        //(RC-based in p) is converted to AA-based
        //"integrating" means for tuples including pos,
        //if tuple pruned for all RCs at AA, it is pruned for that AA
        //this scheme obviously prunes no AA's that have valid confs
        //it might miss some AAs that need to be pruned,
        //but dfs during tuple expansion will handle that
        int numAllowedAtPos[] = p.getNumConfAtPos().clone();
        int integPosNumAA = RCsForAATypes.get(integPos).size();
        numAllowedAtPos[integPos] = integPosNumAA;
        PruningMatrix ans = new PruningMatrix(p.getNumPos(), numAllowedAtPos, p.getPruningInterval());

        //first deal with pruning tuples involving integration pos
        for(int aa=0; aa<integPosNumAA; aa++){
            boolean singlePruned = true;
            for(int rc : RCsForAATypes.get(integPos).get(aa)){
                if(!p.getOneBody(integPos, rc))
                    singlePruned = false;
            }
            ans.setOneBody(integPos, aa, singlePruned);
            for(int pos2=0; pos2<p.getNumPos(); pos2++){
                if(pos2!=integPos){
                    for(int rc2=0; rc2<p.getNumConfAtPos(pos2); rc2++){//may be rc's or aa's depending on pos2, doesn't matter
                        boolean pairPruned = true;
                        for(int rc : RCsForAATypes.get(integPos).get(aa)){
                            if(!p.getPairwise(integPos, rc, pos2, rc2))
                                pairPruned = false;
                        }
                        ans.setPairwise(integPos, aa, pos2, rc2, pairPruned);
                    }
                    //DEBUG!!!  Unlikely that an entire seq triple prunable by RC terms
                    //if there is one leave it to DFS to find
                }
            }
        }

        //then deal with all other tuples
        for(int pos1=0; pos1<p.getNumPos(); pos1++){
            if(pos1!=integPos){
                for(int rc=0; rc<p.getNumConfAtPos(pos1); rc++){
                    ans.setOneBody(pos1, rc, p.getOneBody(pos1,rc));

                    for(int pos2=0; pos2<pos1; pos2++){
                        if(pos2!=integPos){
                            for(int rc2=0; rc2<p.getNumConfAtPos(pos2); rc2++){
                                ans.setPairwise(pos1, rc, pos2, rc2, p.getPairwise(pos1,rc,pos2,rc2));
                                HigherTupleFinder<Boolean> htf = p.getHigherOrderTerms(pos1, rc, pos2, rc2);
                                if(htf!=null){
                                    for(RCTuple prunedTup : htf.listInteractionsWithValue(true) ){
                                        if(!prunedTup.pos.contains(integPos)){
                                            ans.setTupleValue(prunedTup.addRC(pos1,rc).addRC(pos2,rc2), true);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        return ans;
    }




    EnergyMatrix estimateIntegEMat(EnergyMatrix emat, int integPos){
        //Estimates emat integrated at integPos for purposes of choosing tuples
        //hence based on maximum interactions over RCs at AA
        int numAllowedAtPos[] = emat.getNumConfAtPos().clone();
        int integPosNumAA = RCsForAATypes.get(integPos).size();
        numAllowedAtPos[integPos] = integPosNumAA;
        EnergyMatrix ans = new EnergyMatrix(emat.getNumPos(), numAllowedAtPos, emat.getPruningInterval());

        //first deal with pruning tuples involving integration pos
        for(int aa=0; aa<integPosNumAA; aa++){
            double singleE = 0;
            for(int rc : RCsForAATypes.get(integPos).get(aa)){
                if(Math.abs(emat.getOneBody(integPos, rc)) > Math.abs(singleE))
                    singleE = emat.getOneBody(integPos, rc);
            }
            ans.setOneBody(integPos, aa, singleE);
            for(int pos2=0; pos2<emat.getNumPos(); pos2++){
                if(pos2!=integPos){
                    for(int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++){//may be rc's or aa's depending on pos2, doesn't matter
                        double pairE = 0;
                        for(int rc : RCsForAATypes.get(integPos).get(aa)){
                            if(Math.abs(emat.getPairwise(integPos, rc, pos2, rc2)) > Math.abs(pairE))
                                pairE = emat.getPairwise(integPos, rc, pos2, rc2);
                        }
                        ans.setPairwise(integPos, aa, pos2, rc2, pairE);
                    }
                    //DEBUG!!!  For estimates higher-order terms currently don't matter
                }
            }
        }

        //then deal with all other tuples
        for(int pos1=0; pos1<emat.getNumPos(); pos1++){
            if(pos1!=integPos){
                for(int rc=0; rc<emat.getNumConfAtPos(pos1); rc++){
                    ans.setOneBody(pos1, rc, emat.getOneBody(pos1,rc));

                    for(int pos2=0; pos2<pos1; pos2++){
                        if(pos2!=integPos){
                            for(int rc2=0; rc2<emat.getNumConfAtPos(pos2); rc2++){
                                ans.setPairwise(pos1, rc, pos2, rc2, emat.getPairwise(pos1,rc,pos2,rc2));
                                //DEBUG!!  higher-order not used for estimates for now
                            }
                        }
                    }
                }
            }
        }

        return ans;
    }



    public void testVsExhaustiveLUTE(boolean testEPICToo){
        //compare the partition function compute by IBIS
        //to that by exhaustive enumeration of conformations with LUTE
        //If testEPICToo also test vs EPIC directly (much more expensive probably...)
        int numTests = 10;
        Random rand = new Random();

        System.out.println("COMPARING IBIS TO EXHAUSTIVE SEARCH");

        for(int t=0; t<numTests; t++){
            //let's draw a random unpruned sequence
            int seq[] = new int[numPos];

            do {
                for(int pos=0; pos<numPos; pos++)
                    seq[pos] = rand.nextInt(RCsForAATypes.get(pos).size());
            } while(seqPruneMat.isPruned(new RCTuple(seq)));


            System.out.print("SEQUENCE: ");
            for(int pos=0; pos<numPos; pos++)
                System.out.print(AATypes.get(pos).get(seq[pos])+" ");
            System.out.println();

            System.out.println("IBIS FREE ENERGY: "+seqEMat.confE(seq));

            int conf[] = new int[numPos];
            double luteQ = computeExhaustivePartitionFunction(seq, conf, 0, true);
            System.out.println("LUTE FREE ENERGY: "+(-IBISKStarTupleExpander.RT*Math.log(luteQ)));

            if(testEPICToo){
                conf = new int[numPos];
                double epicQ = computeExhaustivePartitionFunction(seq, conf, 0, false);
                System.out.println("EPIC FREE ENERGY: "+(-IBISKStarTupleExpander.RT*Math.log(epicQ)));
            }
        }
    }


    public double computeExhaustivePartitionFunction(int seq[], int conf[], int posOffset,
                                                     boolean useLUTE){
        //compute the partition function using LUTE or EPIC for the specified sequence
        //if posOffset==0, compute only over confs with the specified RCs for pos < posOffset
        if(posOffset==numPos){
            //all RCs defined
            if(lutePruningMatrix.isPruned(new RCTuple(conf)))
                return 0;
            else {
                double E = useLUTE ? luteMatrix.confE(conf) : sp.EPICMinimizedEnergy(conf);
                return Math.exp(-E/IBISKStarTupleExpander.RT);
            }
        }
        else {
            double ans = 0;
            for(int rc : RCsForAATypes.get(posOffset).get(seq[posOffset])){
                if(!lutePruningMatrix.getOneBody(posOffset, rc)){
                    conf[posOffset] = rc;
                    ans += computeExhaustivePartitionFunction(seq,conf,posOffset+1,useLUTE);
                }
            }
            return ans;
        }
    }



    public static void main(String args[]){
        //for testing

        ConfigFileParser cfp = ConfigFileParser.makeFromFilePaths(args);
        cfp.loadData();
        SearchProblem searchSpace = cfp.getSearchProblem();
        searchSpace.pruneMat = new PruningMatrix(searchSpace.confSpace, Double.POSITIVE_INFINITY);
        searchSpace.plugMat = new PolytopeMatrix(searchSpace,true);//Use PLUG to prune

        ObjectIO.writeObject(searchSpace.pruneMat, searchSpace.name+".PRUNEMAT.dat");

        //make LUTE matrix
        searchSpace.loadEnergyMatrix();//The tuple enumerator wants this
        searchSpace.loadEPICMatrix();
        searchSpace.loadTupExpEMatrix();


        ConfSpace confSpace = searchSpace.confSpace;
        KStarBeak curvedBeak = new KStarBeak(searchSpace.tupExpEMat, searchSpace.pruneMat, confSpace);
        curvedBeak.makeUnprunedClusterExpansion();



        //DEBUG!!!  testing vs exhaustive WITH EPIC
        curvedBeak.sp = searchSpace;//for EPIC
        curvedBeak.testVsExhaustiveLUTE(true);

        //then of course the seq-based energy and pruning matrices can be used for design.
        System.out.println("The ibis has completed its work.");
    }
}
