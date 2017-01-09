package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.ConfSpace;

public class LazyEnergyMatrix extends EnergyMatrix {

	private static final long serialVersionUID = -6363593602037378135L;
	
	private static final double DefaultNullVal = Double.NaN;
	
	private SimpleEnergyCalculator ecalc;
	private double nullVal;
	
	public LazyEnergyMatrix(LazyEnergyMatrix other) {
		super(other);
		init(other.ecalc, other.nullVal);
	}
	
	public LazyEnergyMatrix(EnergyMatrix emat, SimpleEnergyCalculator ecalc) {
		this(emat, ecalc, DefaultNullVal);
	}
	
	public LazyEnergyMatrix(EnergyMatrix emat, SimpleEnergyCalculator ecalc, double nullVal) {
		super(emat);
		init(ecalc, nullVal);
	}

	public LazyEnergyMatrix(ConfSpace confSpace, double pruningInterval, SimpleEnergyCalculator ecalc) {
		this(confSpace, pruningInterval, ecalc, DefaultNullVal);
	}
	
	public LazyEnergyMatrix(ConfSpace confSpace, double pruningInterval, SimpleEnergyCalculator ecalc, double nullVal) {
		super(confSpace, pruningInterval);
		init(ecalc, nullVal);
		clear();
	}
	
	private void init(SimpleEnergyCalculator ecalc, double nullVal) {
		this.ecalc = ecalc;
		this.nullVal = nullVal;
	}
	
    @Override
    public Double getOneBody(int res, int conf) {
    	double val = super.getOneBody(res, conf);
    	if (hasVal(val)) {
    		return val;
    	}
    	val = ecalc.calcSingle(res, conf).energy;
    	super.setOneBody(res, conf, val);
    	return val;
    }

	@Override
    public Double getPairwise(int res1, int conf1, int res2, int conf2) {
    	double val = super.getPairwise(res1, conf1, res2, conf2);
    	if (hasVal(val)) {
    		return val;
    	}
    	val = ecalc.calcPair(res1, conf1, res2, conf2).energy;
    	super.setPairwise(res1, conf1, res2, conf2, val);
    	return val;
    }
	
	public boolean hasOneBody(int res, int conf) {
		return hasVal(super.getOneBody(res, conf));
	}
	
	public boolean hasPairwise(int res1, int conf1, int res2, int conf2) {
		return hasVal(super.getPairwise(res1, conf1, res2, conf2));
	}
	
	public void clear() {
		for (int res1=0; res1<getNumPos(); res1++) {
			for (int conf1=0; conf1<getNumConfAtPos(res1); conf1++) {
				
				clearOneBody(res1, conf1);
				
				for (int res2=0; res2<res1; res2++) {
					for (int conf2=0; conf2<getNumConfAtPos(res2); conf2++) {
						
						clearPairwise(res1, conf1, res2, conf2);
					}
				}
			}
		}
	}
	
	public void clearOneBody(int res, int conf) {
		super.setOneBody(res, conf, nullVal);
	}
	
	public void clearPairwise(int res1, int conf1, int res2, int conf2) {
		super.setPairwise(res1, conf1, res2, conf2, nullVal);
	}
	
	private boolean hasVal(double val) {
		
		// is our "null" value NaN? treat is specially since ==,!= don't work
		if (Double.isNaN(nullVal)) {
			return !Double.isNaN(val);
		} else {
			return val != nullVal;
		}
	}
}
