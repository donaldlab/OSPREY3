package edu.duke.cs.osprey.astar.conf.scoring.mplp;

import java.util.Arrays;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class MessageVars {
	
	private RCs rcs;
	private ConfIndex confIndex;
	
	private double[][] sums;
	
	// decision variables, msg_ij(xj)
	// for every (node,node,rotamer) tuple, where nodes are ordered
	private double[][] vars;
	
	public MessageVars(RCs rcs, ConfIndex confIndex) {
		this(rcs, confIndex, true);
	}
	
	public MessageVars(RCs rcs, ConfIndex confIndex, boolean usePrecomputedSums) {
		
		this.rcs = rcs;
		this.confIndex = confIndex;
		
		int n = confIndex.numUndefined;
		
		// allocate space for the sums
		sums = new double[n][];
		for (int posi1=0; posi1<n; posi1++) {
			int pos1 = confIndex.undefinedPos[posi1];
			sums[posi1] = new double[rcs.getNum(pos1)];
			if (usePrecomputedSums) {
				Arrays.fill(sums[posi1], 0);
			} else {
				Arrays.fill(sums[posi1], Double.POSITIVE_INFINITY);
			}
		}
		
		// allocate space for the messages
		vars = new double[n*n][];
		for (int posi1=0; posi1<n; posi1++) {
			for (int posi2=0; posi2<n; posi2++) {
				int pos2 = confIndex.undefinedPos[posi2];
				int index = getNodePairIndex(posi1, posi2);
				vars[index] = new double[rcs.getNum(pos2)];
				Arrays.fill(vars[index], 0);
			}
		}
	}
	
	public void initTraditionalAStar(EnergyMatrix emat) {
		
		for (int posi1=0; posi1<confIndex.numUndefined; posi1++) {
			int pos1 = confIndex.undefinedPos[posi1];
			
			for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
				int rc1 = rcs.get(pos1, rci1);
				
				// init i,i messages with single and defined-undefined energies
				double sum = emat.getOneBody(pos1, rc1);
				for (int posi2=0; posi2<confIndex.numDefined; posi2++) {
					int pos2 = confIndex.definedPos[posi2];
					int rc2 = confIndex.definedRCs[posi2];
					sum += emat.getPairwise(pos1, rc1, pos2, rc2);
				}
				set(posi1, posi1, rci1, sum);
				
				// init i,j messages with undefined-undefined energies
				for (int posi2=0; posi2<confIndex.numUndefined; posi2++) {
					int pos2 = confIndex.undefinedPos[posi2];
					
					if (pos2 < pos1) {
					
						// min over the other RC
						double minEnergy = Double.POSITIVE_INFINITY;
						for (int rc2 : rcs.get(pos2)) {
							minEnergy = Math.min(minEnergy, emat.getPairwise(pos1, rc1, pos2, rc2));
						}
						set(posi2, posi1, rci1, minEnergy);
					
					} else if (pos2 > pos1) {
						
						set(posi2, posi1, rci1, 0);
					}
				}
				
				if (canUsePrecomputedSums(posi1, rci1)) {
					
					sum = 0;
					for (int posi2=0; posi2<confIndex.numUndefined; posi2++) {
						sum += get(posi2, posi1, rci1);
					}
					sums[posi1][rci1] = sum;
				}
			}
		}
	}
	
	public RCs getRCs() {
		return rcs;
	}
	
	public ConfIndex getConfIndex() {
		return confIndex;
	}
	
	private boolean canUsePrecomputedSums(int posi, int rci) {
		return Double.isFinite(sums[posi][rci]);
	}
	
	private boolean canUsePrecomputedSums(int posi, int rci, double val) {
		return canUsePrecomputedSums(posi, rci) && Double.isFinite(val);
	}
	
	public double get(int posi1, int posi2, int rci2) {
		return vars[getNodePairIndex(posi1, posi2)][rci2];
	}
	
	public void set(int posi1, int posi2, int rci2, double val) {
		
		if (Double.isNaN(val)) {
			throw new IllegalArgumentException("val shouldn't be NaN");
		}
		if (val == Double.NEGATIVE_INFINITY) {
			throw new IllegalArgumentException("val shouldn't be -inf");
		}
		
		int index = getNodePairIndex(posi1, posi2);
		if (canUsePrecomputedSums(posi2, rci2, val)) {
			sums[posi2][rci2] -= vars[index][rci2];
			sums[posi2][rci2] += val;
		} else {
			sums[posi2][rci2] = Double.POSITIVE_INFINITY;
		}
		vars[index][rci2] = val;
	}
	
	public double getEnergy(int posi1, int rci1) {
		
		// b_i(xi) = sum_k lambda_ki(xi)
		// energy_posi1(rci1) = sum_{posi2 != posi1} lambda_{posi2,posi1}(rci1)
		
		if (canUsePrecomputedSums(posi1, rci1)) {
			
			return sums[posi1][rci1];
			
		} else {
			
			double sum = 0;
			for (int posi2=0; posi2<confIndex.numUndefined; posi2++) {
				sum += vars[getNodePairIndex(posi2, posi1)][rci1];
			}
			return sum;
		}
	}
	
	public double getEnergyWithout(int posi1, int rci1, int posi2Out) {
		if (canUsePrecomputedSums(posi1, rci1)) {
			
			return getEnergy(posi1, rci1) - get(posi2Out, posi1, rci1);
			
		} else {
			
			double sum = 0;
			for (int posi2=0; posi2<confIndex.numUndefined; posi2++) {
				if (posi2 != posi2Out) {
					sum += vars[getNodePairIndex(posi2, posi1)][rci1];
				}
			}
			return sum;
		}
	}
	
	public double getTotalEnergy() {
		double energy = 0;
		for (int posi1=0; posi1<confIndex.numUndefined; posi1++) {
			int pos1 = confIndex.undefinedPos[posi1];
			
			double minPosEnergy = Double.POSITIVE_INFINITY;
			
			for (int rci1=0; rci1<rcs.getNum(pos1); rci1++) {
				minPosEnergy = Math.min(minPosEnergy, getEnergy(posi1, rci1));
			}
			
			energy += minPosEnergy;
		}
		return energy;
	}

	private int getNodePairIndex(int posi1, int posi2) {
		return posi1*confIndex.numUndefined + posi2;
	}
}
