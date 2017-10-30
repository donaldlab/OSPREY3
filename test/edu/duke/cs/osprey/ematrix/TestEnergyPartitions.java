package edu.duke.cs.osprey.ematrix;

import static edu.duke.cs.osprey.TestBase.*;
import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyPartition;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.structure.PDBIO;

public class TestEnergyPartitions {
	
	private static final double ExpectedEnergy = -89.083499;
	private static final double EnergyEpsilon = 1e-6;
	
	private static SimpleConfSpace confSpace;
	private static RCTuple conf;
	
	@BeforeClass
	public static void beforeClass() {
		
		// read a protein
		Strand strand = new Strand.Builder(PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb")).build();
		
		// make a conf space that allows only the wild type conf
		for (int i=2; i<=6; i++) {
			strand.flexibility.get("A" + i).addWildTypeRotamers();
		}
		confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();
		conf = new RCTuple(new int[] { 0, 0, 0, 0, 0 });
	}
	
	// regardless of energy partition, the total conf energy should be the same
	// and the sum of singles and pairs should equal the total conf energy
	
	@Test
	public void checkPartitions() {
		new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
			.use((ecalc) -> {
				
				double confEnergy = new ConfEnergyCalculator.Builder(confSpace, ecalc)
					.build()
					.calcEnergy(conf, EnergyPartition.makeFragment(confSpace, null, false, conf)).energy;
				assertThat(confEnergy, isAbsolutely(ExpectedEnergy, EnergyEpsilon));
				
				for (EnergyPartition epart : EnergyPartition.values()) {
					
					ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
						.setEnergyPartition(epart)
						.build();

					double partsEnergy = 0;
					for (int pos1=0; pos1<confSpace.positions.size(); pos1++) {
						partsEnergy += confEcalc.calcEnergy(
							new RCTuple(pos1, 0),
							epart.makeSingle(confSpace, null, false, pos1, 0)
						).energy;
						for (int pos2=0; pos2<pos1; pos2++) {
							partsEnergy += confEcalc.calcEnergy(
								new RCTuple(pos1, 0, pos2, 0),
								epart.makePair(confSpace, null, false, pos1, 0, pos2, 0)
							).energy;
						}
					}
					assertThat(partsEnergy, isAbsolutely(ExpectedEnergy, EnergyEpsilon));
				}
			});
	}
}
