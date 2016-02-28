package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.util.Comparator;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings({ "rawtypes", "serial" })
public class KSConf implements Comparable, Serializable {

	public static final double ERROR = 0.0001;

	private int conf[];
	private double minEnergyLB = Double.POSITIVE_INFINITY;
	private double minEnergy = Double.POSITIVE_INFINITY;
	private double minEnergyUB = Double.POSITIVE_INFINITY;
	
	
	public KSConf( int conf[], double energyLB, double energyUB ) {
		this.conf = conf;
		this.minEnergyLB = energyLB;
		this.minEnergyUB = energyUB;
	}


	public double getMinEnergyLB() {
		return minEnergyLB;
	}
	
	
	public double getMinEnergy() {
		return minEnergy;
	}
	
	
	public double getMinEnergyUB() {
		return minEnergyUB;
	}
	
	
	public void setMinEnergy( double e ) {
		minEnergy = e;
	}
	
	
	public void setMinEnergyUB( double e ) {
		minEnergyUB = e;
	}


	public int[] getConf() {
		return conf;
	}


	@Override
	public int compareTo(Object rhs) {
		if( Math.abs(this.minEnergyLB - ((KSConf)rhs).getMinEnergyLB()) <  ERROR ) return 0;

		return this.minEnergyLB > ((KSConf)rhs).getMinEnergyLB() ? 1 : -1;
	}
	
	public class KSConfMinEComparator implements Comparator<KSConf> {

		@Override
		public int compare(KSConf o1, KSConf o2) {
			
			double diff = o1.getMinEnergy() - o2.getMinEnergy();
			
			if(diff > ERROR) return -1;
			
			if(diff < ERROR) return 1;
			
			return 0;
		}
		
	}

}
