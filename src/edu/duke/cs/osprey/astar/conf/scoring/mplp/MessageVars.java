package edu.duke.cs.osprey.astar.conf.scoring.mplp;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class MessageVars {
	
	private RCs rcs;
	private ConfIndex confIndex;
	
	private double[][] sums;
	
	// decision variables, msg_ij(xj)
	// for every node,node,rotamer tuple, where nodes are ordered
	private double[][] vars;
	
	public MessageVars(RCs rcs, ConfIndex confIndex, EnergyMatrix emat) {
		
		this.rcs = rcs;
		this.confIndex = confIndex;
		
		int n = confIndex.getNumUndefined();
		
		// allocate space for the sums
		sums = new double[n][];
		for (int posi1=0; posi1<n; posi1++) {
			int pos1 = confIndex.getUndefinedPos()[posi1];
			sums[posi1] = new double[rcs.getNum(pos1)];
		}
		
		// allocate space for the messages
		vars = new double[n*n][];
		for (int posi1=0; posi1<n; posi1++) {
			for (int posi2=0; posi2<n; posi2++) {
				
				if (posi1 == posi2) {
					continue;
				}
				
				int pos2 = confIndex.getUndefinedPos()[posi2];
				vars[getNodePairIndex(posi1, posi2)] = new double[rcs.getNum(pos2)];
			}
		}
		
		// initialize the sums
		for (int posi1=0; posi1<n; posi1++) {
			int pos1 = confIndex.getUndefinedPos()[posi1];
			
			for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
				int rc1 = rcs.get(pos1, rci1);
		
				// start with single energy
				double sum = emat.getOneBody(pos1, rc1);
				
				// add undefined-defined energies
				for (int posi2=0; posi2<confIndex.getNumDefined(); posi2++) {
					int pos2 = confIndex.getDefinedPos()[posi2];
					int rc2 = confIndex.getDefinedRCs()[posi2];
					sum += emat.getPairwise(pos1, rc1, pos2, rc2);
				}

				sums[posi1][rci1] = sum;
			}
		}
	}
	
	public RCs getRCs() {
		return rcs;
	}
	
	public ConfIndex getConfIndex() {
		return confIndex;
	}
	
	public double get(int posi1, int posi2, int rci2) {
		return vars[getNodePairIndex(posi1, posi2)][rci2];
	}
	
	public void set(int posi1, int posi2, int rci2, double val) {
		int index = getNodePairIndex(posi1, posi2);
		sums[posi2][rci2] -= vars[index][rci2];
		vars[index][rci2] = val;
		sums[posi2][rci2] += val;
	}
	
	public double getEnergy(int posi1, int rci1) {
		
		// b_i(xi) = sum_k lambda_ki(xi)
		// energy_posi1(rci1) = sum_{posi2 != posi1} lambda_{posi2,posi1}(rci1)
		
		return sums[posi1][rci1];
	}
	
	public double getEnergyWithout(int posi1, int rci1, int posi2) {
		return getEnergy(posi1, rci1) - get(posi2, posi1, rci1);
	}
	
	public double getTotalEnergy() {
		double energy = 0;
		for (int posi1=0; posi1<confIndex.getNumUndefined(); posi1++) {
			int pos1 = confIndex.getUndefinedPos()[posi1];
			
			double minPosEnergy = Double.POSITIVE_INFINITY;
			
			for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
				minPosEnergy = Math.min(minPosEnergy, getEnergy(posi1, rci1));
			}
			
			energy += minPosEnergy;
		}
		return energy;
	}

	private int getNodePairIndex(int posi1, int posi2) {
		assert (posi1 != posi2);
		return posi1*confIndex.getNumUndefined() + posi2;
	}
}
