package edu.duke.cs.osprey.sparse;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.TestKStar;
import org.junit.Test;

import java.io.File;

public class TestSparseGraphs {

    @Test
    public void testComputeGraph () {
        SimpleConfSpace confSpace = TestKStar.make2RL0().complex;
        EnergyCalculator ecalc = new EnergyCalculator.Builder(confSpace, new ForcefieldParams())
                .setIsMinimizing(false)
                .build();
        ConfEnergyCalculator confEcalc = new ConfEnergyCalculator.Builder(confSpace, ecalc)
                .setReferenceEnergies(
                        new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalc)
                                .setCacheFile(new File("complex.eref"))
                                .build()
                                .calcReferenceEnergies())
                .build();
        EnergyMatrix emat = new  SimplerEnergyMatrixCalculator.Builder(confEcalc)
                .setCacheFile(new File("complex.emat"))
                .build()
                .calcEnergyMatrix();

        ResidueInteractions graph  = SparseGraphGenerator.generateGraph(confSpace, emat, 0.5, 0.2);
        System.out.println("Done!");
    }
}
