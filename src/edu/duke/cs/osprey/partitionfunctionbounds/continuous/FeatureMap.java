package edu.duke.cs.osprey.partitionfunctionbounds.continuous;

public class FeatureMap {
	
	// kernel that this feature map is constructed from
	Kernel k;
	
	// value in domain from which feature map was generated
	double[] loc; 

	/** 
	 * A feature map is determined by a kernel function and a variable in the
	 * domain
	 * @param k
	 * @param x
	 */
	public FeatureMap(Kernel k, double[] x) {
		this.k = k;
		this.loc = x;
	}
	
	public double eval(double[] y) {
		return k.eval(this.loc, y);
	}
	
	public double[] getLoc() {
		return this.loc;
	}
	
	public Kernel getKernel() {
		return this.k;
	}
}
