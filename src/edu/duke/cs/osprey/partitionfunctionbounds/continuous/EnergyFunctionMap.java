package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.ArrayList;
import java.util.LinkedHashMap;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SearchProblem;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.pruning.PruningMatrix;

/**
 * Represents a mapping between points in the CRMF domain to rotamers and energies
 * in the protein voxel space
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */

public class EnergyFunctionMap {

	SearchProblem search;
	LinkedHashMap<double[], int[]> point2Conf;
	LinkedHashMap<RCTuple, Double> oneBody2Energy;//unary RCTuple -> energy
	LinkedHashMap<RCTuple, Double> pairWise2Energy;//pairwise RCTuple -> energy

	ArrayList<ArrayList<String>> AATypeOptions;//AA type options by position
	RCDatum[][] rcUnary;//2d array of unary RCTuples. each row corresponds to rotamers of a residue
	ArrayList<RCData> rcPairWise;

	public EnergyFunctionMap(SearchProblem search, 
			LinkedHashMap<double[], int[]> point2Conf) {
		this.search = search;
		this.point2Conf = point2Conf;
		oneBody2Energy = new LinkedHashMap<>();
		pairWise2Energy = new LinkedHashMap<>();

		AATypeOptions = new ArrayList<>(); 
		for(int i=0; i < search.confSpace.numPos; ++i) AATypeOptions.add(new ArrayList<>());
		rcUnary = new RCDatum[AATypeOptions.size()][];
		rcPairWise = new ArrayList<>();
	}

	public int getNumResidues() {
		return rcUnary.length;
	}

	/**
	 * Populate unary data
	 */
	public void populateOneBodyRCData() {
		PruningMatrix pmat = search.pruneMat;
		EnergyMatrix emat = search.emat;
		double energy;

		for(int pos = 0; pos < search.confSpace.numPos; ++pos) {
			ArrayList<Integer> posarr = new ArrayList<>(); posarr.add(pos);
			ArrayList<RCTuple> tuples = pmat.unprunedRCTuplesAtPos(posarr);
			rcUnary[pos] = new RCDatum[tuples.size()];

			for(int index = 0; index < tuples.size(); ++index) {
				RCTuple tuple = tuples.get(index);

				energy = emat.getOneBody(tuple.pos.get(0), tuple.RCs.get(0));
				RCDatum rcd = new RCDatum(search.confSpace, tuple, energy);
				rcUnary[pos][index] = rcd;

				String AAType = search.confSpace.posFlex.get(tuple.pos.get(0)).RCs.get(tuple.RCs.get(0)).AAType;
				if(!AATypeOptions.get(pos).contains(AAType)) AATypeOptions.get(pos).add(AAType);
			}
		}
	}
	
	public RCData addPairWise(RCDatum one, RCDatum two) {
		RCData pw = new RCData(one, two);
		rcPairWise.add(pw);
		return pw;
	}

	public int getNumNodes() {
		return rcUnary.length;
	}

	public RCDatum getNode(int pos) {
		return rcUnary[pos][0];
	}

	double getEnergy(RCTuple tup, double[] point) {
		return search.getEnergy(tup, point);
	}

	/**
	 * Get all intra energies
	 */
	public void populateOneBody2Energy() {
		PruningMatrix pmat = search.pruneMat;
		EnergyMatrix emat = search.emat;
		double energy = 0;

		for(int pos = 0; pos < search.confSpace.numPos; ++pos) {

			ArrayList<Integer> posarr = new ArrayList<>(); posarr.add(pos);
			ArrayList<RCTuple> tuples = pmat.unprunedRCTuplesAtPos(posarr);

			for(RCTuple tuple : tuples) {
				energy = emat.getOneBody(tuple.pos.get(0), tuple.RCs.get(0));
				oneBody2Energy.put(tuple, energy);
			}
		}
	}

	/**
	 * Get all pairwise energies
	 */
	public void populatePairWise2Energy() {

		PruningMatrix pmat = search.pruneMat;
		EnergyMatrix emat = search.emat;
		double energy = 0;

		for(int pos1 = 0; pos1 < search.confSpace.numPos; ++pos1) {
			ArrayList<Integer> unpruned1 = pmat.unprunedRCsAtPos(pos1);
			for(int rc1 : unpruned1) {
				for(int pos2 = pos1+1; pos2 < search.confSpace.numPos; ++pos2) {
					ArrayList<Integer> unpruned2 = pmat.unprunedRCsAtPos(pos2);
					for(int rc2 : unpruned2) {
						energy = emat.getPairwise(pos1, rc1, pos2, rc2);
						RCTuple tuple = new RCTuple();
						tuple.pos.add(pos1); tuple.pos.add(pos2); tuple.pos.trimToSize();
						tuple.RCs.add(rc1); tuple.RCs.add(rc2); tuple.RCs.trimToSize();
						pairWise2Energy.put(tuple, energy);
					}
				}
			}
		}
	}

	public double getPairWiseEnergy(RCDatum one, RCDatum two) {
		EnergyMatrix emat = search.emat;
		double energy = emat.getPairwise(one.rcTuple.pos.get(0), one.rcTuple.RCs.get(0), two.rcTuple.pos.get(0), two.rcTuple.RCs.get(0));
		return energy;
	}

	public double[][] getKernelDomainBounds(int pos) {
		if(rcUnary[pos] == null) throw new RuntimeException("ERROR: RCData was not populated at pos "+pos);

		//for now, assuming just one aa per pos...
		//returns n x 2 matrix, where n is the number of dihedrals per AA at pos, and the 
		//2 dimensions are the lower and upper DOF bounds. [][0] and [][1] are lower and 
		//upper bounds, respectively
		double[][] ans = new double[rcUnary[pos][0].getNumDOFs()][2];
		for(int i=0; i<ans.length; ++i) {
			ans[i][0] = Double.POSITIVE_INFINITY;
			ans[i][1] = Double.NEGATIVE_INFINITY;
		}

		for(RCDatum datum : rcUnary[pos]) {
			if(datum == null) throw new RuntimeException("ERROR: RCDatum was not populated");
			double[] dofMin = datum.getDOFMin();
			double[] dofMax = datum.getDOFMax();

			for(int i=0; i<ans.length; ++i) {
				if(dofMin[i]<ans[i][0]) ans[i][0] = dofMin[i];
				if(dofMax[i]>ans[i][1]) ans[i][1] = dofMax[i];
			}
		}

		return ans;
	}
}
