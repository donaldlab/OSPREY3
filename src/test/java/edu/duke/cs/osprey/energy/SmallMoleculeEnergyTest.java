package edu.duke.cs.osprey.energy;

import EDU.oswego.cs.dl.util.concurrent.FJTask;
import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.ConfSpaceIteration;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.energy.forcefield.amber.ForcefieldFileParser;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.kstar.pfunc.GradientDescentPfunc;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.restypes.ResidueTemplateLibrary;
import edu.duke.cs.osprey.structure.PDBIO;
import org.junit.jupiter.api.Test;

import java.nio.file.Paths;
import java.util.List;

public class SmallMoleculeEnergyTest {

    static String dabrafenib = """
HETATM 1721  C1  P06 B 801       0.356  -5.287   5.535  1.00 46.32           C  
HETATM 1722  C2  P06 B 801      -1.074  -3.695   6.429  1.00 45.82           C  
HETATM 1723  N3  P06 B 801      -0.057  -4.584   6.615  1.00 47.01           N  
HETATM 1724  C4  P06 B 801      -1.663  -3.520   5.172  1.00 43.93           C  
HETATM 1725  N6  P06 B 801      -0.172  -5.197   4.290  1.00 45.92           N  
HETATM 1726  C7  P06 B 801      -1.176  -4.301   4.125  1.00 44.81           C  
HETATM 1727  N9  P06 B 801       1.433  -6.144   5.718  1.00 45.68           N  
HETATM 1728  C12 P06 B 801      -1.463  -2.876   7.620  1.00 44.51           C  
HETATM 1729  S13 P06 B 801      -0.416  -2.830   9.023  1.00 48.07           S  
HETATM 1730  C14 P06 B 801      -1.485  -1.707   9.793  1.00 48.17           C  
HETATM 1731  N15 P06 B 801      -2.528  -1.430   9.048  1.00 46.93           N  
HETATM 1732  C16 P06 B 801      -2.568  -2.048   7.879  1.00 44.09           C  
HETATM 1733  C17 P06 B 801      -1.283  -1.112  11.195  1.00 47.87           C  
HETATM 1734  C18 P06 B 801       0.077  -1.527  11.809  1.00 49.44           C  
HETATM 1735  C22 P06 B 801      -2.405  -1.619  12.127  1.00 47.86           C  
HETATM 1736  C26 P06 B 801      -1.332   0.433  11.125  1.00 44.68           C  
HETATM 1737  C30 P06 B 801      -3.770  -1.804   7.019  1.00 42.67           C  
HETATM 1738  C31 P06 B 801      -5.317  -0.234   5.935  1.00 42.58           C  
HETATM 1739  C32 P06 B 801      -4.120  -0.506   6.625  1.00 43.36           C  
HETATM 1740  C33 P06 B 801      -6.127  -1.331   5.600  1.00 42.24           C  
HETATM 1741  C35 P06 B 801      -4.592  -2.871   6.648  1.00 42.72           C  
HETATM 1742  C37 P06 B 801      -5.762  -2.636   5.932  1.00 41.82           C  
HETATM 1743  F39 P06 B 801      -3.280   0.509   6.941  1.00 44.60           F  
HETATM 1744  N40 P06 B 801      -5.622   1.082   5.663  1.00 40.82           N  
HETATM 1745  S42 P06 B 801      -7.030   1.673   5.287  1.00 46.34           S2+
HETATM 1746  C43 P06 B 801      -7.025   1.467   3.498  1.00 47.46           C  
HETATM 1747  C44 P06 B 801      -7.873   0.443   1.471  1.00 47.37           C  
HETATM 1748  C46 P06 B 801      -7.864   0.540   2.860  1.00 47.17           C  
HETATM 1749  C47 P06 B 801      -7.054   1.265   0.703  1.00 47.73           C  
HETATM 1750  C49 P06 B 801      -6.214   2.302   2.714  1.00 48.30           C  
HETATM 1751  C50 P06 B 801      -6.222   2.192   1.324  1.00 48.29           C  
HETATM 1752  F52 P06 B 801      -5.422   3.234   3.286  1.00 49.98           F  
HETATM 1753  F53 P06 B 801      -8.693  -0.267   3.558  1.00 48.53           F  
HETATM 1754  O54 P06 B 801      -8.174   0.918   5.807  1.00 44.79           O1-
HETATM 1755  O55 P06 B 801      -7.098   3.125   5.484  1.00 47.12           O1-
HETATM 1756  H4  P06 B 801      -2.443  -2.800   4.994  1.00  0.00           H  
HETATM 1757  H7  P06 B 801      -1.579  -4.198   3.131  1.00  0.00           H  
HETATM 1758  H91 P06 B 801       1.846  -6.613   4.922  1.00  0.00           H  
HETATM 1759  H92 P06 B 801       1.866  -6.194   6.628  1.00  0.00           H  
HETATM 1760 H181 P06 B 801       0.223  -1.091  12.799  1.00  0.00           H  
HETATM 1761 H182 P06 B 801       0.157  -2.609  11.919  1.00  0.00           H  
HETATM 1762 H183 P06 B 801       0.913  -1.202  11.188  1.00  0.00           H  
HETATM 1763 H221 P06 B 801      -2.278  -1.241  13.141  1.00  0.00           H  
HETATM 1764 H222 P06 B 801      -3.394  -1.314  11.781  1.00  0.00           H  
HETATM 1765 H223 P06 B 801      -2.410  -2.708  12.188  1.00  0.00           H  
HETATM 1766 H261 P06 B 801      -1.197   0.890  12.107  1.00  0.00           H  
HETATM 1767 H262 P06 B 801      -0.553   0.828  10.470  1.00  0.00           H  
HETATM 1768 H263 P06 B 801      -2.283   0.793  10.735  1.00  0.00           H  
HETATM 1769  H33 P06 B 801      -7.053  -1.194   5.068  1.00  0.00           H  
HETATM 1770  H35 P06 B 801      -4.336  -3.881   6.933  1.00  0.00           H  
HETATM 1771  H37 P06 B 801      -6.399  -3.462   5.652  1.00  0.00           H  
HETATM 1772  H44 P06 B 801      -8.525  -0.270   0.993  1.00  0.00           H  
HETATM 1773  H47 P06 B 801      -7.074   1.193  -0.373  1.00  0.00           H  
HETATM 1774  H50 P06 B 801      -5.595   2.838   0.728  1.00  0.00           H  
""";

    static String dabrafenibRotamers = """
1
P06 10 1
H91 N9  C1  N3
N3  C2  C12 C16
C12 C16 C30 C32
C32 C31 N40 S42
C31 N40 S42 C43
N40 S42 C43 C46
N15 C14 C17 C18
C14 C17 C18 H181
C14 C17 C22 H221
C14 C17 C26 H261
-173 -168 -124 -162 -88 111 -176 179 178 180
""";

    static String dabrafenibPrepi = """
    0    0    2

This is a remark line
molecule.res
P06   INT  0
CORRECT     OMIT DU   BEG
  0.0000
   1  DUMM  DU    M    0  -1  -2     0.000      .0        .0      .00000
   2  DUMM  DU    M    1   0  -1     1.449      .0        .0      .00000
   3  DUMM  DU    M    2   1   0     1.523   111.21       .0      .00000
   4  N3    NC    M    3   2   1     1.540   111.208  -180.000 -0.764000
   5  C1    CQ    M    4   3   2     1.353   107.770    56.939  0.835900
   6  N9    N2    B    5   4   3     1.388   116.885  -138.655 -0.911500
   7  H91   H     E    6   5   4     1.012   119.929  -173.148  0.403800
   8  H92   H     E    6   5   4     1.009   118.837    -1.069  0.403800
   9  N6    NC    M    5   4   3     1.355   125.432    43.104 -0.751000
  10  C7    CA    M    9   5   4     1.356   116.371     1.592  0.454200
  11  H7    H4    E   10   9   5     1.078   116.899   177.948  0.032100
  12  C4    CA    M   10   9   5     1.394   122.518    -0.845 -0.363600
  13  H4    HA    E   12  10   9     1.076   120.247   178.284  0.192000
  14  C2    CA    M   12  10   9     1.399   117.233    -0.194  0.563800
  15  C12   C*    M   14  12  10     1.497   122.480   177.202 -0.268900
  16  S13   S     M   15  14  12     1.751   119.766  -163.349 -0.076300
  17  C14   CA    M   16  15  14     1.731    90.241   177.535  0.386700
  18  C17   CT    3   17  16  15     1.536   125.165   179.485 -0.078300
  19  C18   CT    3   18  17  16     1.549   111.928     5.544 -0.089100
  20  H181  HC    E   19  18  17     1.092   111.714   179.195  0.041367
  21  H182  HC    E   19  18  17     1.091   111.732   -60.767  0.041367
  22  H183  HC    E   19  18  17     1.091   111.556    59.297  0.041367
  23  C22   CT    3   18  17  16     1.544   109.152  -114.089 -0.089100
  24  H221  HC    E   23  18  17     1.090   111.290   178.287  0.041367
  25  H222  HC    E   23  18  17     1.091   112.047   -61.081  0.041367
  26  H223  HC    E   23  18  17     1.091   111.401    58.775  0.041367
  27  C26   CT    3   18  17  16     1.547   109.952   125.904 -0.089100
  28  H261  HC    E   27  18  17     1.092   111.927   179.693  0.041367
  29  H262  HC    E   27  18  17     1.092   111.457   -59.927  0.041367
  30  H263  HC    E   27  18  17     1.089   111.956    59.600  0.041367
  31  N15   NB    M   17  16  15     1.311   112.054     0.856 -0.614000
  32  C16   CC    M   31  17  16     1.323   115.302     0.169  0.465800
  33  C30   CA    M   32  31  17     1.498   117.096   177.574 -0.143500
  34  C35   CA    B   33  32  31     1.397   120.011  -119.655 -0.192000
  35  C37   CA    B   34  33  32     1.392   120.153   175.293 -0.100000
  36  C33   CA    S   35  34  33     1.395   120.021    -0.950 -0.238000
  37  H33   HA    E   36  35  34     1.077   117.490  -178.439  0.127000
  38  H37   HA    E   35  34  33     1.080   120.006  -179.601  0.109000
  39  H35   HA    E   34  33  32     1.080   120.304    -3.358  0.120000
  40  C32   CA    M   33  32  31     1.401   120.855    57.524  0.024900
  41  F39   F     E   40  33  32     1.355   118.270     5.759 -0.146900
  42  C31   CA    M   40  33  32     1.408   121.949  -173.319  0.475400
  43  N40   N2    M   42  40  33     1.378   117.991   175.872 -0.934900
  44  S42   SO    M   43  42  40     1.573   127.176  -161.499  1.514000
  45  O54   O     E   44  43  42     1.466   114.855    27.186 -0.674700
  46  O55   O     E   44  43  42     1.467   112.420   161.529 -0.718300
  47  C43   CA    M   44  43  42     1.801   101.072   -87.815 -0.302500
  48  C46   CA    M   47  44  43     1.404   121.697   110.669  0.228400
  49  F53   F     E   48  47  44     1.351   121.760    -1.667 -0.123900
  50  C44   CA    M   48  47  44     1.392   120.217   177.468 -0.227500
  51  H44   HA    E   50  48  47     1.078   119.494  -179.737  0.133500
  52  C47   CA    M   50  48  47     1.391   120.372    -0.154 -0.074000
  53  H47   HA    E   52  50  48     1.079   120.017  -179.286  0.119000
  54  C50   CA    M   52  50  48     1.392   119.935    -0.122 -0.227500
  55  H50   HA    E   54  52  50     1.080   119.958  -179.278  0.133500
  56  C49   CA    M   54  52  50     1.394   120.050    -0.333  0.228400
  57  F52   F     M   56  54  52     1.350   118.694  -178.310 -0.123900


LOOP
   C2   N3
  C16  C12
  C31  C33
  C49  C43

IMPROPER
   N9   N3   C1   N6
   C1  H91   N9  H92
   C4   H7   C7   N6
   C2   C7   C4   H4
  C12   C4   C2   N3
   C2  C16  C12  S13
  C17  N15  C14  S13
  C12  C30  C16  N15
  C35  C32  C30  C16
  C30  C37  C35  H35
  C35  C33  C37  H37
  C37  C31  C33  H33
  C30  C31  C32  F39
  C33  C32  C31  N40
  C46  C49  C43  S42
  C43  C44  C46  F53
  C46  C47  C44  H44
  C44  C50  C47  H47
  C47  C49  C50  H50
  C43  C50  C49  F52

DONE
STOP
""";

    static String dabrafenibPrepiCoordinates = """
P06 54
C1  0.356  -5.287   5.535  
C2  -1.074  -3.695   6.429  
N3  -0.057  -4.584   6.615  
C4  -1.663  -3.520   5.172  
N6  -0.172  -5.197   4.290  
C7  -1.176  -4.301   4.125  
N9  1.433  -6.144   5.718  
C12 -1.463  -2.876   7.620  
S13 -0.416  -2.830   9.023  
C14 -1.485  -1.707   9.793  
N15 -2.528  -1.430   9.048  
C16 -2.568  -2.048   7.879  
C17 -1.283  -1.112  11.195  
C18 0.077  -1.527  11.809  
C22 -2.405  -1.619  12.127  
C26 -1.332   0.433  11.125  
C30 -3.770  -1.804   7.019  
C31 -5.317  -0.234   5.935  
C32 -4.120  -0.506   6.625  
C33 -6.127  -1.331   5.600  
C35 -4.592  -2.871   6.648  
C37 -5.762  -2.636   5.932  
F39 -3.280   0.509   6.941  
N40 -5.622   1.082   5.663  
S42 -7.030   1.673   5.287  
C43 -7.025   1.467   3.498  
C44 -7.873   0.443   1.471  
C46 -7.864   0.540   2.860  
C47 -7.054   1.265   0.703  
C49 -6.214   2.302   2.714  
C50 -6.222   2.192   1.324  
F52 -5.422   3.234   3.286  
F53 -8.693  -0.267   3.558  
O54 -8.174   0.918   5.807  
O55 -7.098   3.125   5.484  
H4  -2.443  -2.800   4.994  
H7  -1.579  -4.198   3.131  
H91 1.846  -6.613   4.922  
H92 1.866  -6.194   6.628  
H181 0.223  -1.091  12.799  
H182 0.157  -2.609  11.919  
H183 0.913  -1.202  11.188  
H221 -2.278  -1.241  13.141  
H222 -3.394  -1.314  11.781  
H223 -2.410  -2.708  12.188  
H261 -1.197   0.890  12.107  
H262 -0.553   0.828  10.470  
H263 -2.283   0.793  10.735  
H33 -7.053  -1.194   5.068  
H35 -4.336  -3.881   6.933  
H37 -6.399  -3.462   5.652  
H44 -8.525  -0.270   0.993  
H47 -7.074   1.193  -0.373  
H50 -5.595   2.838   0.728  
ENDRES
""";

    String dabrafenibFrcmodPath = "/dabrafenib.frcmod";
    String parm96DatPath = "/config/parm96.dat";

    @Test
    public void testComputedEnergyOfDabrafenib() {
        var parser = new ForcefieldFileParser(getClass().getResourceAsStream(parm96DatPath), Paths.get(getClass().getResource(dabrafenibFrcmodPath).getFile()));
        var ffParams = new ForcefieldParams(ForcefieldParams.Forcefield.AMBER, parser);
        var templateLib = new ResidueTemplateLibrary.Builder(ffParams)
                .addTemplateCoords(dabrafenibPrepiCoordinates)
                .addTemplates(dabrafenibPrepi)
                .addRotamers(dabrafenibRotamers)
                .build();

        var strand = PDBIO.read(dabrafenib);
        var dabrafenib = new Strand.Builder(strand)
                .setResidueMutability("B801", List.of(), true, true)
                .setTemplateLibrary(templateLib)
                .build();

        var confSpace = new SimpleConfSpace.Builder()
                .addStrand(dabrafenib)
                .build();
        var ecalc = new EnergyCalculator.Builder(confSpace, ffParams)
                .setParallelism(Parallelism.SingleThreaded)
                .build();
        var referenceEnergies = new SimpleReferenceEnergies.Builder(confSpace, ecalc)
                .build();
        var confEnergyCalculator = new ConfEnergyCalculator.Builder(confSpace, ecalc)
                .setReferenceEnergies(referenceEnergies)
                .build();
        var minimizedEnergyMatrix = new SimplerEnergyMatrixCalculator.Builder(confEnergyCalculator)
                .build()
                .calcEnergyMatrix();

        var discreteEnergyCalc = new EnergyCalculator.Builder(confSpace, ffParams)
                .setParallelism(Parallelism.SingleThreaded)
                .setIsMinimizing(false)
                .build();
        var discreteConfEnergyCalc = new ConfEnergyCalculator.Builder(confSpace, discreteEnergyCalc)
                .setReferenceEnergies(referenceEnergies)
                .build();
        var discreteEnergyMatrix = new SimplerEnergyMatrixCalculator.Builder(discreteConfEnergyCalc)
                .build()
                .calcEnergyMatrix();

        var diff = discreteEnergyMatrix.diff(minimizedEnergyMatrix);

        var rcs = new RCs((ConfSpaceIteration) confSpace);
        var pfn = new GradientDescentPfunc(
                confEnergyCalculator,
                new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                new ConfAStarTree.Builder(minimizedEnergyMatrix, rcs).setTraditional().build(),
                rcs.getNumConformations()
        );

        try (var ctxGroup = ecalc.tasks.contextGroup()) {
            pfn.setInstanceId(0);
            pfn.init(0.00001);
            pfn.putTaskContexts(ctxGroup);
            pfn.compute();
        }
        System.out.println(pfn.makeResult());

        var analyzer = new ConfAnalyzer(confEnergyCalculator);
        var analysis = analyzer.analyze(new int[]{0});

        System.out.println(analysis.toString());
        System.out.println(analysis.breakdownEnergyByPosition(ResidueForcefieldBreakdown.Type.All));
        System.out.println(analysis.toString());
    }
}
