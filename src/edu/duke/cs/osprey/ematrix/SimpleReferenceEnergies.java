/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ematrix;

import java.util.HashMap;
import java.util.Map;

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
	
	public void updateEnergyMatrix(SimpleConfSpace confSpace, EnergyMatrix emat) {
		for (Position pos : confSpace.positions) {
			for (ResidueConf rc : pos.resConfs) {
				String resType = rc.template.name;
				double eref = get(pos.index, resType);
			
				if (eref == Double.POSITIVE_INFINITY) {
					// if all RCs for a residue type have infinite one-body energy
					// (i.e., are impossible),
					// then they stay at infinity after eRef correction
					eref = 0;
				}
				
				// update the emat
				double e = emat.getOneBody(pos.index, rc.index);
				e -= eref;
				emat.setOneBody(pos.index, rc.index, e);
			}
		}
	}
	
	public double getConfEnergy(SimpleConfSpace confSpace, int[] conf) {
		double totERef = 0;
		for (Position pos : confSpace.positions) {
			String resType = pos.resConfs.get(conf[pos.index]).template.name;
			totERef += get(pos.index, resType);
		}
		return totERef;
	}
}
