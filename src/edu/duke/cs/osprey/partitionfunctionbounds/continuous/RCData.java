package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

import edu.duke.cs.osprey.confspace.RCTuple;

public class RCData {
	
	public RCTuple rcTuple;
	public static EnergyFunctionMap EFM;
	private RCDatum one;
	private RCDatum two;
	double energy;
	
	public RCData(RCDatum one, RCDatum two) {
		this.rcTuple = new RCTuple();
		this.rcTuple.pos.addAll(one.getRCTuple().pos);
		this.rcTuple.RCs.addAll(one.getRCTuple().RCs);
		this.rcTuple.pos.addAll(two.getRCTuple().pos);
		this.rcTuple.RCs.addAll(two.getRCTuple().RCs);
		this.one = one;
		this.two = two;
		this.energy = EFM.getPairWiseEnergy(one, two);
	}
	
	public double getEnergy(double[] point) {
		double[] dofMinOne = one.getDOFMin(); double[] dofMaxOne = one.getDOFMax();
		double[] dofMinTwo = two.getDOFMin(); double[] dofMaxTwo = two.getDOFMax();
		double[] x = new double[dofMinOne.length+dofMinTwo.length];
		double sumE = 0;
		int numSamples = EnergyFunctionMap.numSamples;
		
		for(int i=0; i<numSamples; ++i) {
			for(int dof = 0; dof < dofMinOne.length; ++dof) {
				double rand = Math.random();
				if(rand > 0.5) {
					x[dof] = point[dof] + Math.random()*(dofMaxOne[dof]-point[dof]);
				}
				else {
					x[dof] = point[dof] - Math.random()*(point[dof]-dofMinOne[dof]);
				}

				x[dof] = x[dof] > dofMaxOne[dof] ? dofMaxOne[dof] : x[dof];
				x[dof] = x[dof] < dofMinOne[dof] ? dofMinOne[dof] : x[dof];
			}
		}
		
		int offset = dofMinOne.length;
		for(int i=0; i<numSamples; ++i) {
			for(int dof = 0; dof < dofMinTwo.length; ++dof) {
				double rand = Math.random();
				if(rand > 0.5) {
					x[offset+dof] = point[offset+dof] + Math.random()*(dofMaxTwo[dof]-point[offset+dof]);
				}
				else {
					x[offset+dof] = point[offset+dof] - Math.random()*(point[offset+dof]-dofMinTwo[dof]);
				}

				x[offset+dof] = x[offset+dof] > dofMaxTwo[dof] ? dofMaxTwo[dof] : x[offset+dof];
				x[offset+dof] = x[offset+dof] < dofMinTwo[dof] ? dofMinTwo[dof] : x[offset+dof];
				
				sumE += EFM.getEnergy(rcTuple, x);
			}
		}
		
		return (sumE+energy)/(double)(numSamples+1);
	}
	
	public RCTuple getRCTuple() {
		return rcTuple;
	}
	
	public double getEnergy() {
		return energy;
	}
}
