package edu.duke.cs.osprey.kstar;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Comparator;

/**
 * 
 * @author Adegoke Ojewole (ao68@duke.edu)
 *
 */
@SuppressWarnings({ "serial", "rawtypes" })
public class KSConf implements Comparable, Serializable {

	private ArrayList<Integer> conf = null;
	private double minEnergyLB = Double.POSITIVE_INFINITY;
	private double minEnergy = Double.POSITIVE_INFINITY;


	public KSConf( int[] conf, double energyLB ) {

		this.conf = array2List(conf);
		this.minEnergyLB = energyLB;
	}
	
	
	public KSConf( ArrayList<Integer> conf, double energyLB ) {
		this.conf = conf;
		this.minEnergyLB = energyLB;
	}


	public double getMinEnergyLB() {
		return minEnergyLB;
	}


	public double getMinEnergy() {
		return minEnergy;
	}


	public void setMinEnergy( double e ) {
		minEnergy = e;
	}
	
	
	public ArrayList<Integer> getConf() {
		return conf;
	}


	public static int[] list2Array(ArrayList<Integer> in) {
		int[] ans = new int[in.size()];
		for(int i = 0; i < in.size(); ++i) ans[i] = in.get(i);
		return ans;
	}
	
	
	public static ArrayList<Integer> array2List(int[] in) {
		ArrayList<Integer> ans = new ArrayList<>();
		for(int i : in) ans.add(i);
		ans.trimToSize();
		return ans;
	}
	
	
	public int[] getConfArray() {
		return list2Array(conf);
	}


	@Override
	public int compareTo(Object rhs) {
		if( this.minEnergyLB == ((KSConf)rhs).getMinEnergyLB() ) return 0;

		return this.minEnergyLB > ((KSConf)rhs).getMinEnergyLB() ? 1 : -1;
	}


	public class KSConfMinEComparator implements Comparator<KSConf>, Serializable {

		public int compare(KSConf o1, KSConf o2) {
			return o1.getMinEnergy() <= o2.getMinEnergy() ? 1 : -1;
		}	
	}
}
