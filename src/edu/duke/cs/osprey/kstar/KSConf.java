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
	private double energyBound = Double.POSITIVE_INFINITY;
	private double energy = Double.POSITIVE_INFINITY;


	public KSConf( int[] conf, double energyLB ) {

		this.conf = array2List(conf);
		this.energyBound = energyLB;
	}
	
	
	public KSConf( ArrayList<Integer> conf, double energyBound ) {
		this.conf = conf;
		this.energyBound = energyBound;
	}


	public double getEnergyBound() {
		return energyBound;
	}


	public void setEnergyBound( double e ) {
		energyBound = e;
	}
	
	
	public double getEnergy() {
		return energy;
	}


	public void setEnergy( double e ) {
		energy = e;
	}
	
	
	public ArrayList<Integer> getConf() {
		conf.trimToSize();
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
		if( this.energyBound == ((KSConf)rhs).getEnergyBound() ) return 0;

		return this.energyBound > ((KSConf)rhs).getEnergyBound() ? 1 : -1;
	}


	public class KSConfMinEComparator implements Comparator<KSConf>, Serializable {

		public int compare(KSConf o1, KSConf o2) {
			return o1.getEnergy() <= o2.getEnergy() ? 1 : -1;
		}	
	}
}
