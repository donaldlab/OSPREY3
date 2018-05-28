/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ibis;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.plug.PolytopeMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.tupexp.BasicPruningTupleExpander;
import edu.duke.cs.osprey.tupexp.TupleExpander;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * A tuple expander wrt sequences for some residues and RC's for others
 *
 *
 * @author mhall44
 */
public class IBISKStarTupleExpander extends BasicPruningTupleExpander {

    EnergyMatrix baseMatrix;//energy matrix with one less res integrated

    PruningMatrix basePruneMat;
    //the super pruning matrix describes pruning in the actual ("integrated") space of this expansion


    int curIntegPos;
    ArrayList<ArrayList<Integer> > RCsforAA;//which RCs are available for each amino acid type for curIntegPos


    public static final double RT = 1.9891/1000.0 * 298.15;//RT in kcal/mol..see PoissonBoltzmannEnergy



    public IBISKStarTupleExpander(int curIntegPos, PruningMatrix basePruneMat, PruningMatrix integPruneMat,
                                  EnergyMatrix baseMatrix, ArrayList<ArrayList<Integer>> RCsforAA) {
        super(integPruneMat, null);//will not be needing PLUG, that was handled at the LUTE fitting stage presumably
        this.baseMatrix = baseMatrix;
        this.basePruneMat = basePruneMat;
        this.curIntegPos = curIntegPos;
        this.RCsforAA = RCsforAA;
    }





    @Override
    public double scoreAssignmentList(int[] assignmentList) {
        double boltz = 0;
        int mann[] = assignmentList.clone();
        for(int rc : RCsforAA.get(assignmentList[curIntegPos])){
            mann[curIntegPos] = rc;
            if(!basePruneMat.isPruned(new RCTuple(mann)))//DEBUG!!  could be slow
                boltz += Math.exp(-baseMatrix.confE(mann)/RT);
        }
        return -RT*Math.log(boltz);
    }

    //Although the conf space is atypical for a BasicPruningTupleExpander
    //(having AA's instead of RC's at some positions),
    //we can still use the super pruning-related methods because we provide
    //a suitable "integrated" pruning matrix

}
