package edu.duke.cs.osprey.sharkstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarNode;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.markstar.framework.MARKStarNode;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.BeforeClass;
import org.junit.Test;
import org.ojalgo.matrix.transformation.Rotation;

import java.io.FileNotFoundException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static edu.duke.cs.osprey.sharkstar.TestSHARKStar.loadFromCFS;


public class TestSHARKStarNodeScorer {

    private static String design = "test-resources/3bua_B_10res_4.363E+11.cfs";
    private static MultiSequenceSHARKStarNode.Node rootNode;
    private static RCs rcs;
    private static int[] order = { 4, 5, 6, 0, 2, 1, 3};
    private static List<MultiSequenceSHARKStarNode.Node> levelOne;
    private static BoltzmannCalculator bc;
    private static ConfIndex confIndex;

    private static SHARKStarNodeScorer lbScorer;
    private static SHARKStarNodeScorer ubScorer;
    private static AStarScorer lbGScorer;
    private static AStarScorer ubGScorer;

    @BeforeClass
    public static void beforeClass() throws FileNotFoundException {
        SimpleConfSpace confSpace = loadFromCFS(design).protein.makeFlexibleCopy();

        ForcefieldParams ffparams = new ForcefieldParams();


        // Define the minimizing energy calculator
        EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(Parallelism.makeCpu(2))
                .build();
        // Define the rigid energy calculator
        EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(Parallelism.makeCpu(2))
                .setIsMinimizing(false)
                .build();
        // Define the method to make confEnergyCalculators (with reference energies)
        BBSHARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
            return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                    .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
                            .build()
                            .calcReferenceEnergies()
                    )
                    .build();
        };
        // Make the confEnergyCalculators
        ConfEnergyCalculator rigidConfEcalc = confEcalcFactory.make(confSpace, rigidEcalc);
        ConfEnergyCalculator minimizingConfEcalc = confEcalcFactory.make(confSpace, minimizingEcalc);

        // Make the energy matrices
        EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc).build().calcEnergyMatrix();
        EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc).build().calcEnergyMatrix();
        EnergyMatrix correctionEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat);

        Sequence wildType = confSpace.makeWildTypeSequence();
        rcs = wildType.makeRCs(confSpace);
        rootNode = new MultiSequenceSHARKStarNode.Node(confSpace.positions.size(), 0, new MathTools.DoubleBounds());
        levelOne = getLevelOne(rootNode);

        bc = new BoltzmannCalculator(new MathContext(BigDecimal.ROUND_HALF_UP));

        lbScorer = new SHARKStarNodeScorer(minimizingEmat, false);
        ubScorer = new SHARKStarNodeScorer(rigidEmat, true);
        lbGScorer = new PairwiseGScorer(minimizingEmat);
        ubGScorer = new PairwiseGScorer(rigidEmat);

        confIndex = new ConfIndex(rcs.getNumPos());
        System.out.println("RCs:"+rcs.toString());
    }

    private static List<MultiSequenceSHARKStarNode.Node> getLevelOne(MultiSequenceSHARKStarNode.Node rootNode){
        List<MultiSequenceSHARKStarNode.Node> children = new ArrayList<>();
        int nextPos = order[0];
        for(int nextRC : rcs.get(nextPos)){
            MultiSequenceSHARKStarNode.Node child = rootNode.assign(nextPos, nextRC);
            children.add(child);
        }
        return children;
    }

    @Test
    public void testRootBound(){
        rootNode.index(confIndex);
        //int[] debugConf = {-1, -1, -1, -1, -1, -1, -1};
        //lbScorer.setDebugConf(debugConf);
        double lb = lbGScorer.calc(confIndex, rcs) + lbScorer.calc(confIndex,rcs);
        double ub = ubGScorer.calc(confIndex, rcs) + ubScorer.calc(confIndex,rcs);

        System.out.println(String.format("%s bounds: [%.3e, %.3e] --> Energy [%.20f, %.20f]", rootNode.confToString(), bc.calc(ub), bc.calc(lb), lb, ub));
    }

    @Test
    public void testRootBound_allatonce(){
        rootNode.index(confIndex);
        //int[] debugConf = {-1, -1, -1, -1, -1, -1, -1};
        //lbScorer.setDebugConf(debugConf);
        BigDecimal pfuncUpper = lbScorer.calcCombinedScore(confIndex,rcs);
        BigDecimal pfuncLower =  ubScorer.calcCombinedScore(confIndex,rcs);

        System.out.println(String.format("%s bounds: [%.3e, %.3e] --> Energy [%.20f, %.20f]", rootNode.confToString(), pfuncLower.doubleValue(), pfuncUpper.doubleValue(), bc.freeEnergy(pfuncUpper),bc.freeEnergy(pfuncLower) ));
    }

    @Test
    public void testBoltzmannAndFreeEnergy(){
        rootNode.index(confIndex);
        BigDecimal pfuncUpper = lbScorer.calcCombinedScore(confIndex,rcs);
        double energy = bc.freeEnergy(pfuncUpper);
        BigDecimal weightedEnergy = bc.calc(energy);
        double energy2 = bc.freeEnergy(weightedEnergy);
        BigDecimal weightedEnergy2 = bc.calc(energy2);

        System.out.println(String.format("Directly calculated pfunc: %.3e --> \n" +
                "energy via bc.freeEnergy: %.3f --> \n" +
                        "pfunc via bc.calc: %.3e --> \n" +
                        "energy via bc.freeEnergy: %.3f --> \n" +
                        "pfunc via bc.calc: %.3e",
                pfuncUpper.doubleValue(), energy, weightedEnergy, energy2, weightedEnergy2));

        BigDecimal testValue = BigDecimal.valueOf(7e75).setScale(64);
        double testenergy = bc.freeEnergy(testValue);
        BigDecimal testValue2 = bc.calc(testenergy);
        System.out.println(String.format("test value: %.3e --> \n" +
                        "energy via bc.freeEnergy: %.3f --> \n" +
                        "pfunc via bc.calc: %.3e --> \n",
                testValue, testenergy, testValue2));
    }

    @Test
    public void testLevel1Bound(){
        BigDecimal pfuncSum_upper = BigDecimal.ZERO;
        BigDecimal pfuncSum_lower = BigDecimal.ZERO;
        for (MultiSequenceSHARKStarNode.Node node : levelOne){
            node.index(confIndex);
            double g_lb = lbGScorer.calc(confIndex, rcs);
            double g_ub = ubGScorer.calc(confIndex, rcs);
            double h_lb = lbScorer.calc(confIndex,rcs);
            double h_ub = ubScorer.calc(confIndex,rcs);
            double lb = g_lb + h_lb;
            double ub = g_ub + h_ub;
            System.out.println(String.format("%s bounds: [%.3e, %.3e] --> Energy [%.3f, %.3f] Gscore: [%.3f, %.3f] Hscore: [%.3f, %.3f]",
                    node.confToString(), bc.calc(ub), bc.calc(lb), lb, ub, g_lb, g_ub, h_lb, h_ub));
            pfuncSum_upper = pfuncSum_upper.add(bc.calc(lb));
            pfuncSum_lower = pfuncSum_lower.add(bc.calc(ub));
        }
        System.out.println(String.format("Level 1 sum bounds: [%.3e, %.3e]", pfuncSum_lower.doubleValue(), pfuncSum_upper.doubleValue()));
    }

    @Test
    public void testLevel1Bound_allatonce(){
        BigDecimal pfuncSum_upper = BigDecimal.ZERO;
        BigDecimal pfuncSum_lower = BigDecimal.ZERO;
        for (MultiSequenceSHARKStarNode.Node node : levelOne){
            node.index(confIndex);
            BigDecimal pfuncUpper = lbScorer.calcCombinedScore(confIndex, rcs);
            BigDecimal pfuncLower = ubScorer.calcCombinedScore(confIndex, rcs);
            System.out.println(String.format("%s bounds: [%.3e, %.3e] --> Energy [%.3f, %.3f]",
                    node.confToString(), pfuncLower, pfuncUpper, bc.freeEnergy(pfuncUpper), bc.freeEnergy(pfuncLower)));
            pfuncSum_upper = pfuncSum_upper.add(pfuncUpper);
            pfuncSum_lower = pfuncSum_lower.add(pfuncLower);
        }
        System.out.println(String.format("Level 1 sum bounds: [%.3e, %.3e]", pfuncSum_lower.doubleValue(), pfuncSum_upper.doubleValue()));
    }

    private class TestContext{
        private String design;
        private String state;
        private int[] order;

        private MultiSequenceSHARKStarNode.Node rootNode;
        private List<MultiSequenceSHARKStarNode.Node> list;

        private SimpleConfSpace confSpace;
        private RCs rcs;
        private BoltzmannCalculator bc;
        private ConfIndex confIndex;

        private SHARKStarNodeScorer lbScorer;
        private SHARKStarNodeScorer ubScorer;
        private AStarScorer lbGScorer;
        private AStarScorer ubGScorer;

        private TestContext( String design, String state, int[] order ){
            this.design = design;
            this.state = state;
            this.order = order;
        }

        private void init() throws FileNotFoundException{
            // Make confspace
            if(state.equals("complex"))
                this.confSpace = loadFromCFS(design).complex;
            else if(state.equals("protein"))
                this.confSpace = loadFromCFS(design).protein;
            else
                this.confSpace = loadFromCFS(design).ligand;
            // Define forcefield parameters
            ForcefieldParams ffparams = new ForcefieldParams();
            // Define the minimizing energy calculator
            EnergyCalculator minimizingEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                    .setParallelism(Parallelism.makeCpu(2))
                    .build();
            // Define the rigid energy calculator
            EnergyCalculator rigidEcalc = new EnergyCalculator.Builder(confSpace, ffparams)
                    .setParallelism(Parallelism.makeCpu(2))
                    .setIsMinimizing(false)
                    .build();
            // Define the method to make confEnergyCalculators (with reference energies)
            BBSHARKStar.ConfEnergyCalculatorFactory confEcalcFactory = (confSpaceArg, ecalcArg) -> {
                return new ConfEnergyCalculator.Builder(confSpaceArg, ecalcArg)
                        .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpaceArg, minimizingEcalc)
                                .build()
                                .calcReferenceEnergies()
                        )
                        .build();
            };
            // Make the confEnergyCalculators
            ConfEnergyCalculator rigidConfEcalc = confEcalcFactory.make(confSpace, rigidEcalc);
            ConfEnergyCalculator minimizingConfEcalc = confEcalcFactory.make(confSpace, minimizingEcalc);
            // Make the energy matrices
            EnergyMatrix rigidEmat = new SimplerEnergyMatrixCalculator.Builder(rigidConfEcalc).build().calcEnergyMatrix();
            EnergyMatrix minimizingEmat = new SimplerEnergyMatrixCalculator.Builder(minimizingConfEcalc).build().calcEnergyMatrix();
            EnergyMatrix correctionEmat = new UpdatingEnergyMatrix(confSpace, minimizingEmat);

            Sequence wildType = confSpace.makeWildTypeSequence();
            this.rcs = wildType.makeRCs(confSpace);
            rootNode = new MultiSequenceSHARKStarNode.Node(confSpace.positions.size(), 0, new MathTools.DoubleBounds());
            levelOne = getLevelOne(rootNode);

            this.bc = new BoltzmannCalculator(new MathContext(BigDecimal.ROUND_HALF_UP));
            this.lbScorer = new SHARKStarNodeScorer(minimizingEmat, false);
            this.ubScorer = new SHARKStarNodeScorer(rigidEmat, true);
            this.lbGScorer = new PairwiseGScorer(minimizingEmat);
            this.ubGScorer = new PairwiseGScorer(rigidEmat);
            this.confIndex = new ConfIndex(rcs.getNumPos());
        }
    }
}
