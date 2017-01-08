package edu.duke.cs.osprey.astar.kadee;

import java.util.List;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
/**
*
* @author PGC
*/
public class EmatBoundMAP2 extends EnergyMatrix {
    EnergyMatrix parent;
    double onebody_subtraction_matrix_ [][];
    double pairwise_subtraction_matrix_ [][][][];

    public EmatBoundMAP2(EnergyMatrix parent, double onebody_subtraction_matrix [][], double pairwise_subtraction_matrix [][][][]) {
        this.parent = parent;
        onebody_subtraction_matrix_ = onebody_subtraction_matrix;
        pairwise_subtraction_matrix_ = pairwise_subtraction_matrix;
    }

    @Override
    public double confE(int conf[]) {
        //value of energy represented in energy matrix, for specified conformation
        //expressed as residue-specific RC indices (as in the storage matrices)
        return getInternalEnergy(new RCTuple(conf)) + parent.getConstTerm();
    }

    @Override
    public Double getOneBody(int res, int index) {

         return parent.getOneBody(res, index) - this.onebody_subtraction_matrix_[res][index];

    }

    @Override
    public Double getPairwise(int res1, int index1, int res2, int index2) {

        return parent.getPairwise(res1, index1, res2, index2) - this.pairwise_subtraction_matrix_[res1][index1][res2][index2];
    }

    @Override
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

                //TODO ADD SUPPORT FOR HOT
            }
        }

        return E;
    }

    @Override
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

    @Override
    public int numPos() {
        return parent.numPos();
    }

    @Override
    public HigherTupleFinder<Double> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
        //working with residue-specific RC indices directly.  
        return null;
    }
}
