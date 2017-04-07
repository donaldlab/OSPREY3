/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import java.util.HashMap;
import java.util.Map;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.Position;
import edu.duke.cs.osprey.confspace.SimpleConfSpace.ResidueConf;

/**
 * This class stores a reference energy for each residue position & AA type
 * (which is the minimum intra-RC energy at that position & AA type)
 */
public class SimpleReferenceEnergies {
	
	private Map<String,Double> energies;
	
	public SimpleReferenceEnergies() {
		energies = new HashMap<>();
	}
	
	public Double get(int pos, String resType) {
		return energies.get(makeKey(pos, resType));
	}
	
	public void set(int pos, String resType, double val) {
		energies.put(makeKey(pos, resType), val);
	}
	
	private String makeKey(int pos, String resType) {
		return "" + pos + "-" + resType;
	}
	
	private double getOffset(int pos, String resType) {
		
		double energy = get(pos, resType);
		
		if (Double.isFinite(energy)) {
			
			// NOTE: negate reference energies here, so they can be added later like normal energy offsets
			return -energy;
			
		} else {
			// if all RCs for a residue type have infinite one-body energy
			// (i.e., are impossible),
			// then they stay at infinity after eRef correction
			return Double.POSITIVE_INFINITY;
		}
	}
	
	public double getConfEnergy(SimpleConfSpace confSpace, int[] conf) {
		double energy = 0;
		for (Position pos : confSpace.positions) {
			ResidueConf rc = pos.resConfs.get(conf[pos.index]);
			energy += getOffset(pos.index, rc.template.name);
		}
		return energy;
	}
	
	public double getFragmentEnergy(SimpleConfSpace confSpace, RCTuple frag) {
		double energy = 0;
		for (int i=0; i<frag.size(); i++) {
			Position pos = confSpace.positions.get(frag.pos.get(i));
			ResidueConf rc = pos.resConfs.get(frag.RCs.get(i));
			energy += getOffset(pos.index, rc.template.name);
		}
		return energy;
	}
}
