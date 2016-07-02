/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupleMatrix;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

/**
 *
 * @author mhall44
 */
public class EnergyMatrix extends TupleMatrix<Double> {

    private double constTerm = 0;

    //we may want to have reference energies associated with this matrix
    ReferenceEnergies eRefMat = null;

    public EnergyMatrix(ConfSpace cSpace, double pruningInterval) {
        //For a normal, precomputed energy matrix we expect infinite pruning interval
        //(matrix valid for all RCs),
        //but EnergyMatrix objects from tup-exp may only be valid for a finite pruning interval
        super(cSpace, pruningInterval, 0.);//default higher interaction for energies is 0
    }

    public EnergyMatrix(int numPos, int[] numRCsAtPos, double pruningInterval) {
        super(numPos, numRCsAtPos, pruningInterval, 0.);
    }

    public EnergyMatrix(TupleMatrix<Double> tupMat) {
        super(tupMat);
    }
    
    public EnergyMatrix(){
        super(0.);
    }
    
    public double confE(int conf[]) {
        //value of energy represented in energy matrix, for specified conformation
        //expressed as residue-specific RC indices (as in the storage matrices)
        return getInternalEnergy(new RCTuple(conf)) + constTerm;
    }

    public double getInternalEnergy(RCTuple tup) {
        //internal energy of a tuple of residues when they're in the specified RCs

        int numPosInTuple = tup.pos.size();
        double E = 0;

        for (int indexInTuple = 0; indexInTuple < numPosInTuple; indexInTuple++) {
            int posNum = tup.pos.get(indexInTuple);
            int RCNum = tup.RCs.get(indexInTuple);

            double intraE = getOneBody(posNum, RCNum);
            E += intraE;

            for (int index2 = 0; index2 < indexInTuple; index2++) {
                int pos2 = tup.pos.get(index2);
                int rc2 = tup.RCs.get(index2);

                double pairwiseE = getPairwise(posNum, RCNum, pos2, rc2);
                E += pairwiseE;

                HigherTupleFinder<Double> htf = getHigherOrderTerms(posNum, RCNum, pos2, rc2);
                if (htf != null) {
                    double higherOrderE = internalEHigherOrder(tup, index2, htf);
                    E += higherOrderE;
                }
            }
        }

        return E;
    }

    double internalEHigherOrder(RCTuple tup, int curIndex, HigherTupleFinder<Double> htf) {
        //Computes the portion of the internal energy for tuple tup
        //that consists of interactions in htf (corresponds to some sub-tuple of tup)
        //with RCs whose indices in tup are < curIndex
        double E = 0;
        ArrayList<Integer> interactingPos = htf.getInteractingPos();

        for (int ipos : interactingPos) {

            //see if ipos is in tup with index < curIndex
            int iposIndex = -1;
            for (int ind = 0; ind < curIndex; ind++) {
                if (tup.pos.get(ind) == ipos) {
                    iposIndex = ind;
                    break;
                }
            }

            if (iposIndex > -1) {//ipos interactions need to be counted
                int iposRC = tup.RCs.get(iposIndex);
                E += htf.getInteraction(ipos, iposRC);

                //see if need to go up to highers order again...
                HigherTupleFinder htf2 = htf.getHigherInteractions(ipos, iposRC);
                if (htf2 != null) {
                    E += internalEHigherOrder(tup, iposIndex, htf2);
                }
            }
        }

        return E;
    }

    /*
     public double getPairwiseE(int res1, int AA1, int rot1, int res2, int AA2, int rot2){
     //lookup by residue number, amino-acid index (for the residue's rotamer library),
     //and RC index for that residue and AA type (rotamer index for AA type if rigid backbone,
     //otherwise defined in the ConfSearchSpace)
     //RETURN ERROR IF RES1 AND RES2 ARE NOT SINGLE RESIDUES
     return getPairwise(res1, index1, res2, index2);
     }
     */
    //intra+shell similar...
    public double getConstTerm() {
        return constTerm;
    }

    public void setConstTerm(double constTerm) {
        this.constTerm = constTerm;
    }

    public double[][] topPairwiseInteractions() {
        //return the top absolute values of the pairwise interactions
        //between all pairs of positions
        int numPos = oneBody.size();

        double strongestPairE[][] = new double[numPos][numPos];

        for (int pos = 0; pos < numPos; pos++) {
            for (int pos2 = 0; pos2 < pos; pos2++) {
                for (int rc = 0; rc < numRCsAtPos(pos); rc++) {
                    for (int rc2 = 0; rc2 < numRCsAtPos(pos2); rc2++) {
                        strongestPairE[pos][pos2] = Math.max(strongestPairE[pos][pos2], Math.abs(getPairwise(pos, rc, pos2, rc2)));
                        strongestPairE[pos2][pos] = strongestPairE[pos][pos2];
                    }
                }
            }
        }

        return strongestPairE;
    }

    public ReferenceEnergies geteRefMat() {
        return eRefMat;
    }

    //TODO: CAN SPEED UP Using PruneMat
    //Add higher-order term functionality
    public void negatePairwiseEnergies(boolean[][] interactionGraph) {
        for (int i = 0; i < numPos(); i++) {
            for (int j = 0; j < i; j++) {
                if (interactionGraph[i][j]) {
                    for (int rotI = 0; rotI < numRCsAtPos(i); rotI++) {
                        for (int rotJ = 0; rotJ < numRCsAtPos(j); rotJ++) {
                            double currentE = this.getPairwise(i, rotI, j, rotJ);
                            this.setPairwise(i, rotI, j, rotJ, -currentE);
                        }
                    }
                }
            }
        }
    }

    /**
     * HMN: This method is used for partial search spaces. Given an interaction
     * graph, this method updates the energy matrix to only have non-zero
     * entries for the pairwise energies defined by the interactionGraph
     *
     * @param interactionGraph interaction graph that maps posI and posJ to true
     * if we should keep the pairwise energies between posI and posJ
     */
    public void updateMatrixCrossTerms(boolean[][] interactionGraph) {
        int numPos = this.oneBody.size();

        int[] numRCsPerPos = new int[numPos];
        for (int pos = 0; pos < numPos; pos++) {
            numRCsPerPos[pos] = this.oneBody.get(pos).size();
        }

        for (int i = 0; i < numPos; i++) {
            //Set all one-body terms to "zero"
            for (int rot = 0; rot < numRCsPerPos[i]; rot++) {
                this.setOneBody(i, rot, 0.0);
            }

            for (int j = 0; j < i; j++) {
                //If they don't interact, set all pairwise to "zero"
                //This is the defaultHigherInteraction (e.g. false, null, 0)
                if (!interactionGraph[i][j]) {
                    for (int rcI = 0; rcI < numRCsPerPos[i]; rcI++) {
                        for (int rcJ = 0; rcJ < numRCsPerPos[j]; rcJ++) {
                            //Set pairwise to "zero"
                            this.setPairwise(i, rcI, j, rcJ, 0.0);
                            //now we need to check if we should include higher-order terms
                            HigherTupleFinder<Double> htf = getHigherOrderTerms(i, rcI, j, rcJ);
                            if (htf != null) {
                                //Get all positions that interact with i or j
                                ArrayList<Integer> interactingPositions = new ArrayList<>();
                                for (int pos = 0; pos < numPos; pos++) {
                                    if (interactionGraph[i][pos] || interactionGraph[j][pos]) {
                                        interactingPositions.add(pos);
                                    }
                                }
                                HigherTupleFinder<Double> crossHtf = htf.getHigherOrderCrossTerms(interactingPositions, numRCsPerPos);
                                this.higherTerms.get(i).get(j).get(rcI).set(rcJ, crossHtf);
                            }
                        }
                    }
                }
            }
        }
    }

    public void addCrossTermInternalEnergies(EnergyMatrix boundEmat, EnergyMatrix unboundEmat, ArrayList<Integer> posNumsToSubtractFrom, HashMap<Integer, Integer> boundPosNumToUnboundPosNum) {
        for (int pos : posNumsToSubtractFrom) {
            int unboundPos = boundPosNumToUnboundPosNum.get(pos);
            for (int rc = 0; rc < numRCsAtPos(pos); rc++) {
                double boundE = boundEmat.getOneBody(pos, rc);
                double unboundE = unboundEmat.getOneBody(unboundPos, rc);
                this.oneBody.get(pos).set(rc, boundE - unboundE);
            }
        }
    }

    public void addInternalEnergies(EnergyMatrix originalEmat, ArrayList<Integer> posToAddE) {
        for (int pos : posToAddE) {
            for (int rc = 0; rc < this.numRCsAtPos(pos); rc++) {
                double currentE = this.getOneBody(pos, rc);
                double toAddE = originalEmat.getOneBody(pos, rc);
                this.setOneBody(pos, rc, currentE + toAddE);
            }
        }
    }

    /**
     * HMN: Negates the energies - useful for finding Global Maximum Energy
     * Conformation
     */
    public void negateEnergies() {
        for (int i = 0; i < this.numPos(); i++) {
            for (int rot = 0; rot < this.numRCsAtPos(i); rot++) {
                double currentE = this.getOneBody(i, rot);
                this.setOneBody(i, rot, -currentE);
            }

            for (int j = 0; j < i; j++) {
                for (int rotI = 0; rotI < this.numRCsAtPos(i); rotI++) {
                    for (int rotJ = 0; rotJ < this.numRCsAtPos(j); rotJ++) {
                        double pairE = this.getPairwise(i, rotI, j, rotJ);
                        this.setPairwise(i, rotI, j, rotJ, -pairE);
                        //Negate higher-order tuple energies
                        HigherTupleFinder<Double> htf = this.getHigherOrderTerms(i, rotI, j, rotJ);
                        if (htf != null) {
                            negateHTFEnergies(htf);
                        }
                    }
                }
            }
        }
    }

    /**
     * HMN: Negates the values of all entries in the higher order tuple,
     * including the higher tuples that it holds
     *
     * @param htf the higher order tuple
     */
    private void negateHTFEnergies(HigherTupleFinder<Double> htf) {
        ArrayList<Integer> interactingPos = htf.getInteractingPos();
        ArrayList<TreeMap<Integer, Double>> interactions = htf.getInteractions();
        ArrayList<TreeMap<Integer, HigherTupleFinder<Double>>> higher = htf.getHigher();
        for (int i = 0; i < interactingPos.size(); i++) {
            TreeMap<Integer, Double> interactionsAtPos = interactions.get(i);
            //Iterate over entries in the map to negate the energies
            Iterator iterI = interactionsAtPos.entrySet().iterator();
            while (iterI.hasNext()) {
                Map.Entry<Integer, Double> pair = (Map.Entry<Integer, Double>) iterI.next();
                pair.setValue(-pair.getValue());
            }

            //Now we recurse on the higher order tuples
            //NB: Since we only include four position tuples, this should only go one
            /// level deep. If this changes, it could cause a stack overflow
            TreeMap<Integer, HigherTupleFinder<Double>> higherAtPos = higher.get(i);
            Iterator iterJ = higherAtPos.entrySet().iterator();
            while (iterJ.hasNext()) {
                Map.Entry<Integer, HigherTupleFinder<Double>> higherPair = (Map.Entry<Integer, HigherTupleFinder<Double>>) iterJ.next();
                negateHTFEnergies(higherPair.getValue());
            }
        }
    }

}
