package edu.duke.cs.osprey.ematrix.compiled;


import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.confspace.compiled.PosInterDist;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.Test;

public class TestEmatCalculator {

	@Test
	public void test() {

		ConfSpace confSpace = new ConfSpace(FileTools.readResourceBytes("/confSpaces/dipeptide.5hydrophobic.ccs.toml.xz"));
		CPUConfEnergyCalculator confEcalc = new CPUConfEnergyCalculator(confSpace);
		EnergyMatrix emat = new EmatCalculator.Builder(confEcalc)
			.setMinimize(false)
			.setPosInterDist(PosInterDist.DesmetEtAl1992)
			.build()
			.calc();

		// TEMP
		System.out.println(emat.toString());

		// TODO: NEXTTIME: not getting any pair energies, no atom pairs in the ccs for some reason...
	}
}
