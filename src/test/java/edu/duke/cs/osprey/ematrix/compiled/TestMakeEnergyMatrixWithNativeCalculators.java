package edu.duke.cs.osprey.ematrix.compiled;

import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.energy.compiled.CPUConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.CudaConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.NativeConfEnergyCalculator;
import edu.duke.cs.osprey.gpu.Precision;
import edu.duke.cs.osprey.tools.FileTools;
import org.junit.jupiter.api.Test;

import java.util.function.Supplier;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assumptions.assumeTrue;

public class TestMakeEnergyMatrixWithNativeCalculators {

    // This confspace has 0 mutations, 0 amino acid flexibility, but translation/rotation on one of the two molecules.
    private static final ConfSpace confSpace = ConfSpace.fromBytes(FileTools.readResourceBytes("/confSpaces/6og7.transrot.minimal.ccsx"));

    private void testCanMakeEnergyMatrixWithCalculator(Supplier<ConfEnergyCalculator> calculatorSupplier) {
        try (var calculator = calculatorSupplier.get()) {
            var emat = new EmatCalculator.Builder(calculator)
                    .build()
                    .calc();

            assertNotNull(emat);
        }
    }

    @Test
    public void testCanMakeEnergyMatrixWithCpuCalculator() {testCanMakeEnergyMatrixWithCalculator(() -> new CPUConfEnergyCalculator(confSpace));}

    @Test
    public void testCanMakeEnergyMatrixWithNative64Calculator() { testCanMakeEnergyMatrixWithCalculator(() -> new NativeConfEnergyCalculator(confSpace, Precision.Float64)); }

    @Test
    public void testCanMakeEnergyMatrixWithNative32Calculator() { testCanMakeEnergyMatrixWithCalculator(() -> new NativeConfEnergyCalculator(confSpace, Precision.Float32)); }

    @Test
    public void testCanMakeEnergyMatrixWithCuda64Calculator() {
        assumeTrue(CudaConfEnergyCalculator.isSupported());
        testCanMakeEnergyMatrixWithCalculator(() -> new CudaConfEnergyCalculator(confSpace, Precision.Float64));
    }

    @Test
    public void testCanMakeEnergyMatrixWithCuda32Calculator() {
        assumeTrue(CudaConfEnergyCalculator.isSupported());
        testCanMakeEnergyMatrixWithCalculator(() -> new CudaConfEnergyCalculator(confSpace, Precision.Float32));
    }
}
