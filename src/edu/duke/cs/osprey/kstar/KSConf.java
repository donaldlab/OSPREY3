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

	private ArrayList<Integer> conf; 
	private double energyBound;
	private double energy;


	public KSConf( int[] conf, double energyLB ) {
		this(array2List(conf), energyLB, Double.POSITIVE_INFINITY);
	}
	
	
	public KSConf(int[] conf, double energyBound, double energy) {
		this(array2List(conf), energyBound, energy);
	}
	
	
	public KSConf( ArrayList<Integer> conf, double energyBound ) {
		this(conf, energyBound, Double.POSITIVE_INFINITY);
	}
	
	
	public KSConf(ArrayList<Integer> conf, double energyBound, double energy) {
		this.conf = conf;
		this.energyBound = energyBound;
		this.energy = energy;
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
