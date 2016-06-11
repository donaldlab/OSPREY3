package edu.duke.cs.osprey.astar.conf.scoring.mplp;

import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

public class MessageVars {
	
	private RCs rcs;
	private ConfIndex confIndex;
	private EnergyMatrix emat;
	
	// decision variables, msg_ij(xj)
	// for every node,node,rotamer tuple, where nodes are ordered
	private double[][] vars;
	
	public MessageVars(RCs rcs, ConfIndex confIndex, EnergyMatrix emat) {
		
		this.rcs = rcs;
		this.confIndex = confIndex;
		this.emat = emat;
		
		int numPos = rcs.getNumPos();
		
		// allocate space for the messages
		vars = new double[numPos*numPos][];
		for (int i=0; i<numPos; i++) {
			for (int j=0; j<numPos; j++) {
				
				if (i == j) {
					continue;
				}
				
				vars[getNodePairIndex(i, j)] = new double[rcs.get(j).length];
			}
		}
	}
	
	public RCs getRCs() {
		return rcs;
	}
	
	public ConfIndex getConfIndex() {
		return confIndex;
	}
	
	public double get(int pos1, int pos2, int rci2) {
		return vars[getNodePairIndex(pos1, pos2)][rci2];
	}
	
	public void set(int pos1, int pos2, int rci2, double val) {
		vars[getNodePairIndex(pos1, pos2)][rci2] = val;
	}
	
	public double getEnergy(int pos1, int rci1) {
		
		// b_i(xi) = sum_k lambda_ki(xi)
		// energy_pos1(rci1) = sum_{pos2 != pos1} lambda_{pos2,pos1}(rci1)
		
		// TODO: speed this up by precomputing sums
		
		// start with single energy
		int rc1 = rcs.get(pos1)[rci1];
		double sum = emat.getOneBody(pos1, rc1);
		
		// add undefined-defined energies
		for (int i=0; i<confIndex.getNumDefined(); i++) {
			int pos2 = confIndex.getDefinedPos()[i];
			int rc2 = confIndex.getDefinedRCs()[i];
			sum += emat.getPairwise(pos1, rc1, pos2, rc2);
		}
		
		// undefined-undefined energies
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos2 = confIndex.getUndefinedPos()[i];
			if (pos2 != pos1) {
				sum += get(pos2, pos1, rci1);
			}
		}
		
		return sum;
	}
	
	public double getEnergyWithout(int pos1, int rci1, int pos2) {
		return getEnergy(pos1, rci1) - get(pos2, pos1, rci1);
	}
	
	public double getTotalEnergy() {
		double energy = 0;
		for (int i=0; i<confIndex.getNumUndefined(); i++) {
			int pos1 = confIndex.getUndefinedPos()[i];
			
			double minPosEnergy = Double.POSITIVE_INFINITY;
			
			for (int rci1=0; rci1<rcs.get(pos1).length; rci1++) {
				minPosEnergy = Math.min(minPosEnergy, getEnergy(pos1, rci1));
			}
			
			energy += minPosEnergy;
		}
		return energy;
	}

	private int getNodePairIndex(int pos1, int pos2) {
		assert (pos1 != pos2);
		return pos1*rcs.getNumPos() + pos2;
	}
}
