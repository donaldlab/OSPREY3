/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.ConfSpaceSuper;
import edu.duke.cs.osprey.confspace.PositionConfSpaceSuper;
import edu.duke.cs.osprey.confspace.RC;
import edu.duke.cs.osprey.confspace.SuperRC;
import edu.duke.cs.osprey.confspace.SuperRCTuple;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.dof.DegreeOfFreedom;
import edu.duke.cs.osprey.ematrix.epic.EPICFitter;
import edu.duke.cs.osprey.ematrix.epic.EPICSettings;
import edu.duke.cs.osprey.ematrix.epic.EPoly;
import edu.duke.cs.osprey.ematrix.epic.FitParams;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.energy.EnergyFunctionGenerator;
import edu.duke.cs.osprey.energy.MultiTermEnergyFunction;
import edu.duke.cs.osprey.handlempi.MPISlaveTask;
import edu.duke.cs.osprey.minimization.CCDMinimizer;
import edu.duke.cs.osprey.minimization.MolecEObjFunction;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.Residue;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author hmn5
 */
//This class is meant to precalculate the minimized energy of an energy term
//(e.g. intra+shell+pairwise between and within a position)
//it should be suitable for calculating at a slave node
///It allows for super-RCs
public class TermECalculatorSuper implements MPISlaveTask {

    boolean doingEPIC;//doing EPIC fit instead of just minimum computation
    boolean doingIntra;//doing just intra energy (only relevant for one-body energies)
    //if false, one-body energies will be intra+shell

    PruningMatrix pruneMat;//pruning matrix so EPIC doesn't handle pruned RCs or pairs; this can be null 
    //for scalar (non-EPIC) energy calculations

    EPICSettings epicSettings = null;//needed if doing EPIC

    //HMN: For Super-RCs allow for confSpaceSuper and replace int[] res with int[] pos
    ConfSpaceSuper confSpaceSuper;
    boolean useSuperRCs = false;
    int pos[];
    MultiTermEnergyFunction termE; //replaces termE

    //HMN: End addition here 
    ///It may be better for this to be its own class
    //We'll be calculating either one-body or pairwise energies
    //as either scalar or EPIC energies
    //So we're trying to fill in one of these four lists
    double constTerm;
    ArrayList<Double> oneBodyE = new ArrayList<>();
    ArrayList<ArrayList<Double>> twoBodyE = new ArrayList<>();
    ArrayList<EPoly> oneBodyPoly = new ArrayList<>();
    ArrayList<ArrayList<EPoly>> twoBodyPoly = new ArrayList<>();

    //HMN: TermECalculator for input ConfSpaceSuper
    public TermECalculatorSuper(ConfSpaceSuper s, ArrayList<Residue> shellResidues, boolean doEPIC, boolean doIntra, PruningMatrix prm, EPICSettings es, int... posToCalc) {

        confSpaceSuper = s;
        useSuperRCs = true;
        doingEPIC = doEPIC;
        doingIntra = doIntra;
        pos = posToCalc;
        pruneMat = prm;
        epicSettings = es;

        termE = new MultiTermEnergyFunction();

        if (pos.length == 0) {//shell-shell energy
            for (int shell1 = 0; shell1 < shellResidues.size(); shell1++) {
                Residue shellRes1 = shellResidues.get(shell1);
                EnergyFunction intraE = EnvironmentVars.curEFcnGenerator.singleResEnergy(shellRes1);
                termE.addTerm(intraE);
                for (int shell2 = 0; shell2 < shell1; shell2++) {
                    Residue shellRes2 = shellResidues.get(shell2);
                    EnergyFunction pairE = EnvironmentVars.curEFcnGenerator.resPairEnergy(shellRes1, shellRes2);
                    termE.addTerm(pairE);
                }
            }
        } else {
            //Since we do not simply have single and pair residues, we use a multi-term energy function
            PositionConfSpaceSuper firstPosition = confSpaceSuper.posFlexSuper.get(posToCalc[0]);
            int numResAtFirstPos = firstPosition.resList.size();
            if (doingIntra) { //And intra and pairwise terms WITHIN the position
                for (int resNum1 = 0; resNum1 < numResAtFirstPos; resNum1++) {
                    Residue res1 = firstPosition.resList.get(resNum1);
                    EnergyFunction eTerm1 = EnvironmentVars.curEFcnGenerator.singleResEnergy(res1);
                    termE.addTerm(eTerm1);

                    for (int resNum2 = 0; resNum2 < resNum1; resNum2++) {
                        Residue res2 = firstPosition.resList.get(resNum2);
                        EnergyFunction eTerm2 = EnvironmentVars.curEFcnGenerator.resPairEnergy(res1, res2);
                        termE.addTerm(eTerm2);
                    }
                }
            } else {
                if (pos.length == 1) {//intra+shell+pairwise
                    for (int resNum1 = 0; resNum1 < numResAtFirstPos; resNum1++) {
                        Residue res1 = firstPosition.resList.get(resNum1);
                        EnergyFunction eTerm1 = EnvironmentVars.curEFcnGenerator.intraAndShellEnergy(res1, shellResidues);
                        termE.addTerm(eTerm1);

                        for (int resNum2 = 0; resNum2 < resNum1; resNum2++) {
                            Residue res2 = firstPosition.resList.get(resNum2);
                            EnergyFunction eTerm2 = EnvironmentVars.curEFcnGenerator.resPairEnergy(res1, res2);
                            termE.addTerm(eTerm2);
                        }
                    }
                } else if (pos.length == 2) {//pairwise between resi in pos1 and resj in pos2
                    for (Residue res1 : firstPosition.resList) {
                        PositionConfSpaceSuper secondPosition = confSpaceSuper.posFlexSuper.get(posToCalc[1]);
                        for (Residue res2 : secondPosition.resList) {
                            EnergyFunction eTerm = EnvironmentVars.curEFcnGenerator.resPairEnergy(res1, res2);
                            termE.addTerm(eTerm);
                        }
                    }
                } else {
                    throw new RuntimeException("ERROR: Can only precompute energy for 1- and 2-body terms");
                }
            }
        }
    }

    @Override
    public Object doCalculation() {
        if (pos.length == 0) {//shell-shell calculation
            shellCalc();
            if (doingEPIC) {
                //return oneBodyPoly
                throw new RuntimeException("ERROR: Currently EPIC is not supported with super-RCs");
            } else {
                return constTerm;
            }
        } else if (pos.length == 1) {//1-body calculation

            oneBodyCalc();

            if (doingEPIC) {
                //return oneBodyPoly
                throw new RuntimeException("ERROR: Currently EPIC is not supported with super-RCs");
            } else {
                return oneBodyE;
            }
        } else if (pos.length == 2) {//pairwise calculation

            twoBodyCalc();

            if (doingEPIC) {
                //return pairwisePoly;
                throw new RuntimeException("ERROR: Currently EPIC is not supported with super-RCs");
            } else {
                return twoBodyE;
            }
        } else {
            throw new RuntimeException("ERROR: Trying to precompute term for " + pos.length + " bodies");
        }
    }

    public void shellCalc() {
        //no super-RCs, input empty tuple
        calcSuperTupleEnergy(new SuperRCTuple());
    }

    public void oneBodyCalc() {
        ArrayList<SuperRC> superRCList = confSpaceSuper.posFlexSuper.get(pos[0]).superRCs;

        for (int superRCNum = 0; superRCNum < superRCList.size(); superRCNum++) {
            SuperRCTuple superRCTup = new SuperRCTuple(pos[0], superRCNum);
            calcSuperTupleEnergy(superRCTup);
        }
    }

    public void twoBodyCalc() {
        //list minimized one-body energies for all the RCs in res

        ArrayList<SuperRC> RCList1 = confSpaceSuper.posFlexSuper.get(pos[0]).superRCs;
        ArrayList<SuperRC> RCList2 = confSpaceSuper.posFlexSuper.get(pos[1]).superRCs;

        for (int firstRCNum = 0; firstRCNum < RCList1.size(); firstRCNum++) {

            for (int secondRCNum = 0; secondRCNum < RCList2.size(); secondRCNum++) {
                SuperRCTuple superRCTup = new SuperRCTuple(pos[0], firstRCNum, pos[1], secondRCNum);
                calcSuperTupleEnergy(superRCTup);
            }
        }
    }

    public void calcSuperTupleEnergy(SuperRCTuple superRCs) {

        boolean skipTuple = false;
        double minEnergy = Double.POSITIVE_INFINITY;
        EPoly EPICFit = null;

        if (superRCs.pos.size() == 2) {//pair: need to check for parametric incompatibility
            //If there are DOFs spanning multiple residues, then parametric incompatibility
            //is whether the pair is mathematically possible (i.e. has a well-defined voxel)
            SuperRC rc1 = confSpaceSuper.posFlexSuper.get(superRCs.pos.get(0)).superRCs.get(superRCs.superRCs.get(0));
            SuperRC rc2 = confSpaceSuper.posFlexSuper.get(superRCs.pos.get(1)).superRCs.get(superRCs.superRCs.get(1));
            if (rc1.isParametricallyIncompatibleWith(rc2)) {
                skipTuple = true;
            }
        }

        if (pruneMat != null) {
            if (pruneMat.isPruned(superRCs)) {
                skipTuple = true;
            }
        }

        if (!skipTuple) {
            MolecEObjFunction mof = new MolecEObjFunction(termE, confSpaceSuper, superRCs);

            DoubleMatrix1D bestDOFVals;

            if (mof.getNumDOFs() > 0) {//there are continuously flexible DOFs to minimize
                CCDMinimizer ccdMin = new CCDMinimizer(mof, true);
                bestDOFVals = ccdMin.minimize();
            } else {
                bestDOFVals = DoubleFactory1D.dense.make(0);
            }
            minEnergy = mof.getValue(bestDOFVals);

            if (doingEPIC) {
                //EPICFit = compEPICFit(mof,minEnergy,bestDOFVals,RCs);
                throw new RuntimeException("ERROR: EPIC not supported with super-RCs");
            }
        }
        int numBodies = superRCs.pos.size();
        if (numBodies == 0) {//shell -shell 
            constTerm = minEnergy;
        } else if (numBodies == 1) {//one-body is one-position

            if (doingEPIC) {
                throw new RuntimeException("ERROR: EPIC not supported with super-RCs");
                //oneBodyPoly.add(EPICFit);
            } else {
                oneBodyE.add(minEnergy);
            }
        } else if (numBodies == 2) {//two-body term between two positoins
            int firstSuperRCNum = superRCs.superRCs.get(0);

            if (doingEPIC) {
                throw new RuntimeException("ERROR: EPIC not supported with super-RCs");
                /*
                 if(pairwisePoly.size()<=firstRCNum)//this is the first term for RC #firstRCNum
                 pairwisePoly.add(new ArrayList<EPoly>());
                
                 pairwisePoly.get(firstRCNum).add(EPICFit);   
                 */
            } else {
                if (twoBodyE.size() <= firstSuperRCNum) {
                    twoBodyE.add(new ArrayList<Double>());
                }
                twoBodyE.get(firstSuperRCNum).add(minEnergy);
            }
        } else {//only 1- and 2-body precomputations supported here
            throw new RuntimeException("ERROR: Trying to precompute term for " + numBodies + " bodies");
        }
    }
}
