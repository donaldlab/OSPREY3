/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import java.math.BigDecimal;
import java.math.BigInteger;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;

/**
 *
 * @author Anna Lowegard(anna.lowegard@duke.edu)
 * Adegoke Ojewole (ao68@duke.edu)
 */
public class PartitionFunction {
	
	private String sequence;
	private BigDecimal z;
	private BigInteger numConfs;
	private BoltzmannCalculator bc;
	
	public PartitionFunction() {
		sequence = "";
		z = BigDecimal.ZERO;
		numConfs = BigInteger.ZERO;
		bc = new BoltzmannCalculator();
	}
	
	public void addEnergy(double energy) {
		//only include good, negative energies in calculations
		if (energy < 0) {
			BigDecimal boltzmannWeight = bc.calc(energy);
			z = z.add(boltzmannWeight);
			//increment numConfs with each conformation
			numConfs = numConfs.add(BigInteger.ONE);
		}
	}
	
	public void setSequence(String sequence) {
		this.sequence = sequence;
	}
	
	public String getSequence() {
		return sequence;
	}
	
	public BigDecimal getZ() {
		return z;
	}
	
	public BigInteger getNumConfs() {
		return numConfs;
	}
	
	public String toString() {
		return "seq: " + getSequence() + ", " + String.format("z*: %12e", z) + ", " + String.format("confs: %d", numConfs);
	}
	
    
}
