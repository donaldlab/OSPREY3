package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;

public class NegatedEnergyMatrix extends ProxyEnergyMatrix {

	public NegatedEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix target) {
		super(confSpace, target);
	}

	@Override
	public double getConstTerm() {
		return -super.getConstTerm();
	}

	@Override
	public void setConstTerm(double val) {
		super.setConstTerm(-val);
	}

	@Override
	public Double getOneBody(int pos, int rc) {
		return -super.getOneBody(pos, rc);
	}

	@Override
	public void setOneBody(int pos, int rc, Double val) {
		super.setOneBody(pos, rc, -val);
	}

	@Override
	public Double getPairwise(int pos1, int rc1, int pos2, int rc2) {
		return -super.getPairwise(pos1, rc1, pos2, rc2);
	}

	@Override
	public void setPairwise(int pos1, int rc1, int pos2, int rc2, Double val) {
		super.setPairwise(pos1, rc1, pos2, rc2, -val);
	}
}
