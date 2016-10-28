package edu.duke.cs.osprey.astar.kadee;

import java.util.List;
import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.HigherTupleFinder;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

// Energy matrix with virtual residues representing different amino acids.
public class EmatUnboundMAP2 extends EnergyMatrix {

    EnergyMatrix parent;
    int num_virtual_residues_;
    int rcs_per_virtual_residue_[];
    int map_from_virtual_res_to_parent_res_ [];
    int map_from_virtual_rc_to_parent_rc_ [][];

    public EmatUnboundMAP2(EnergyMatrix parent, int num_virtual_residues, int rcs_per_virtual_residue[], 
    		int map_from_virtual_res_to_parent_res [], int map_from_virtual_rc_to_parent_rc[][]) {

        this.parent = parent;
        num_virtual_residues_ = num_virtual_residues;
        rcs_per_virtual_residue_ = rcs_per_virtual_residue;
        map_from_virtual_res_to_parent_res_ = map_from_virtual_res_to_parent_res;
        map_from_virtual_rc_to_parent_rc_ = map_from_virtual_rc_to_parent_rc;
    }

    @Override
    public double confE(int conf[]) {
        //value of energy represented in energy matrix, for specified conformation
        //expressed as residue-specific RC indices (as in the storage matrices)
        return getInternalEnergy(new RCTuple(conf)) + parent.getConstTerm();
    }

    @Override
    public Double getOneBody(int virtual_res, int virtual_rc) {
    	// We need to map the virtual residue to the actual values in the parent energy matrix.
    	int parent_res = map_from_virtual_res_to_parent_res_[virtual_res];
    	int parent_rc = map_from_virtual_rc_to_parent_rc_[virtual_res][virtual_rc]; 			
        return parent.getOneBody(parent_res, parent_rc);
        
    }

    @Override
    public Double getPairwise(int virtual_res1, int virtual_rc1, int virtual_res2, int virtual_rc2) {
    	int parent_res1 = map_from_virtual_res_to_parent_res_[virtual_res1];
    	int parent_rc1 = map_from_virtual_rc_to_parent_rc_[virtual_res1][virtual_rc1]; 
    	int parent_res2 = map_from_virtual_res_to_parent_res_[virtual_res2];
    	int parent_rc2 = map_from_virtual_rc_to_parent_rc_[virtual_res2][virtual_rc2]; 
    	if(parent_res1 == parent_res2){
    		return 0.0;
    	}
    	else{
    		return parent.getPairwise(parent_res1, parent_rc1, parent_res2, parent_rc2);
    	}
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
        return this.num_virtual_residues_;
    }

    @Override
    public HigherTupleFinder<Double> getHigherOrderTerms(int res1, int index1, int res2, int index2) {
        //working with residue-specific RC indices directly.  
        return null;
    }
}
