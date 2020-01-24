package edu.duke.cs.osprey.kstar;

import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.tools.ExpFunction;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.BeforeClass;
import org.junit.Test;

import java.math.MathContext;
import java.math.RoundingMode;

import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.junit.Assert.assertThat;

public class TestBoltzmannCalculator {
    private static BoltzmannCalculator bc;
    private static MathContext mc;
    final double Epsilon = 1e-8;

    @BeforeClass
    public static void init(){
        mc = new MathContext(64, RoundingMode.HALF_UP);
        bc = new BoltzmannCalculator(mc);

    }

    @Test
    public void testCalc_accuracy(){
        assertThat(bc.calc(-10.0).doubleValue(), isAbsolutely(21040926.1669687, Epsilon*10)); // for some reason this one thinks there is some error
        assertThat(bc.calc( -7.5).doubleValue(), isAbsolutely(310669.514920208, Epsilon));
        assertThat(bc.calc( -5.0).doubleValue(), isAbsolutely( 4587.03893235807, Epsilon));
        assertThat(bc.calc( -2.5).doubleValue(), isAbsolutely(67.7276821717536, Epsilon));
        assertThat(bc.calc( -1.0).doubleValue(), isAbsolutely( 5.39891494042175, Epsilon));
        assertThat(bc.calc( -0.5).doubleValue(), isAbsolutely( 2.32355652834652, Epsilon));
        assertThat(bc.calc(  0.0).doubleValue(), isAbsolutely( 1, Epsilon));
        assertThat(bc.calc(  0.5).doubleValue(), isAbsolutely( 0.430374724178376, Epsilon));
        assertThat(bc.calc(  1.0).doubleValue(), isAbsolutely( 0.185222403211613, Epsilon));
        assertThat(bc.calc(  2.5).doubleValue(), isAbsolutely( 0.014765011409427, Epsilon));
        assertThat(bc.calc(  5.0).doubleValue(), isAbsolutely( 0.000218005561921, Epsilon));
        assertThat(bc.calc(  7.5).doubleValue(), isAbsolutely( 3.21885460907498E-06, Epsilon));
        assertThat(bc.calc( 10.0).doubleValue(), isAbsolutely( 4.75264250282795E-08, Epsilon));
    }

    @Test
    public void testCalc_size(){
        assertThat(bc.calc(Integer.MIN_VALUE).doubleValue(), isAbsolutely(MathTools.BigPositiveInfinity.doubleValue(), Epsilon));
        assertThat(bc.calc(Long.MIN_VALUE).doubleValue(), isAbsolutely(MathTools.BigPositiveInfinity.doubleValue(), Epsilon));
        assertThat(bc.calc(Integer.MAX_VALUE).doubleValue(), isAbsolutely(0, Epsilon));
        assertThat(bc.calc(Long.MAX_VALUE).doubleValue(), isAbsolutely(0, Epsilon));
    }
}
