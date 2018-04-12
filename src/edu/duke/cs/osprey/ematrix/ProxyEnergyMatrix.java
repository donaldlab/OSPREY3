package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public class ProxyEnergyMatrix extends EnergyMatrix {

	public final EnergyMatrix target;

	public ProxyEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
		super(confSpace);
		this.target = target;
	}

	@Override
	protected void allocate(int numOneBody, int numPairwise) {
		// don't allocate anything
		// we're just going to proxy all methods to the target emat
	}

	@Override
	public Double getOneBody(int pos, int rc) {
		return target.getOneBody(pos, rc);
	}

	@Override
	public void setOneBody(int pos, int rc, Double val) {
		target.setOneBody(pos, rc, val);
	}

	@Override
	public Double getPairwise(int pos1, int rc1, int pos2, int rc2) {
		return target.getPairwise(pos1, rc1, pos2, rc2);
	}

	@Override
	public void setPairwise(int pos1, int rc1, int pos2, int rc2, Double val) {
		target.setPairwise(pos1, rc1, pos2, rc2, val);
	}

	// TODO: need to proxy anything else?
}
