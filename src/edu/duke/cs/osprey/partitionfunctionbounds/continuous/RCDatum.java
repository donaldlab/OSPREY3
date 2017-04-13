package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.RCTuple;

public class RCDatum {

	private ConfSpace confSpace;
	public RCTuple rcTuple;
	private double energy;
	public static EnergyFunctionMap EFM;

	public RCDatum(ConfSpace confSpace, RCTuple rcTuple, double energy) {
		this.confSpace = confSpace;
		this.rcTuple = rcTuple;
		this.energy = energy;
	}
	
	public int getNumDOFs() {
		return confSpace.posFlex.get(rcTuple.pos.get(0)).RCs.get(rcTuple.RCs.get(0)).DOFmin.size();
	}

	public double[] getDOFMin() {
		ArrayList<Double> DOFMinArr = confSpace.posFlex.get(rcTuple.pos.get(0)).RCs.get(rcTuple.RCs.get(0)).DOFmin;
		double[] ans = new double[DOFMinArr.size()]; 
		for(int i=0; i<ans.length; ++i) ans[i] = DOFMinArr.get(i);
		return ans;
	}

	public double[] getDOFMax() {
		ArrayList<Double> DOFMaxArr = confSpace.posFlex.get(rcTuple.pos.get(0)).RCs.get(rcTuple.RCs.get(0)).DOFmax;
		double[] ans = new double[DOFMaxArr.size()]; 
		for(int i=0; i<ans.length; ++i) ans[i] = DOFMaxArr.get(i);
		return ans;
	}

	public double getEnergy(double[] point) {
		//check that dihedrals are within dofmin and dofmax
		//sample around dihedrals; return energy average
		double[] dofMin = getDOFMin();
		double[] dofMax = getDOFMax();
		double[] x = new double[dofMin.length];
		double sumE = 0;
		int numSamples = 25;

		for(int i=0; i<numSamples; ++i) {
			for(int dof = 0; dof < dofMin.length; ++dof) {
				double rand = Math.random();
				if(rand > 0.5) {
					x[dof] = point[dof] + Math.random()*(dofMax[dof]-point[dof]);
				}
				else {
					x[dof] = point[dof] - Math.random()*(point[dof]-dofMin[dof]);
				}

				x[dof] = x[dof] > dofMax[dof] ? dofMax[dof] : x[dof];
				x[dof] = x[dof] < dofMin[dof] ? dofMin[dof] : x[dof];

				sumE += EFM.getEnergy(rcTuple, x);
			}
		}
		
		return sumE/(double)numSamples;
	}

	public RCTuple getRCTuple() {
		return rcTuple;
	}

	public double getEnergy() {
		return energy;
	}
}
