package edu.duke.cs.osprey.sharkstar_refactor;

import cern.colt.matrix.DoubleMatrix1D;
import edu.duke.cs.osprey.astar.conf.ConfIndex;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.conf.scoring.AStarScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.PairwiseRigidGScorer;
import edu.duke.cs.osprey.astar.conf.scoring.ScorerFactory;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.DependentScoreCorrector;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.gmec.ConfAnalyzer;
import edu.duke.cs.osprey.minimization.ObjectiveFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.sharkstar.SHARKStarNodeScorer;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.tools.MathTools;
import org.junit.BeforeClass;
import org.junit.Test;
import static edu.duke.cs.osprey.TestBase.isAbsolutely;
import static org.junit.Assert.assertThat;

import java.util.*;

public class TestUpperBoundCorrections {
    private static Molecule metallochaperone;
    private static Molecule protein_1a0r;

    private static final List<TestConfInfo> flexibleConfList = new ArrayList<>();
    private static final List<TestConfInfo> fullConfList = new ArrayList<>();

    private static class TestConfInfo{
        public int[] assignments;
        public double ELB;
        public double EUB;
        public double minE;

        public TestConfInfo(int[] assignments, double ELB, double EUB, double minE){
            this.assignments = assignments;
            this.ELB = ELB;
            this.EUB = EUB;
            this.minE = minE;
        }
    }

    @BeforeClass
    public static void beforeClass() {
        ArrayList<int[]> flexAssList = new ArrayList<>();
        ArrayList<int[]> fullAssList = new ArrayList<>();

        flexAssList.add(new int[] {4,2,4,3,2,3,1});
        fullAssList.add(new int[] {-1,-1,4,-1,2,4,3,2,3,1});
        flexAssList.add(new int[] {4,4,0,3,3,4,1});
        fullAssList.add(new int[] {-1,-1,4,-1,4,0,3,3,4,1});
        flexAssList.add(new int[] {4,4,0,3,3,4,1});
        fullAssList.add(new int[] {-1,-1,4,-1,4,0,3,3,4,1});
        flexAssList.add(new int[] {4,4,0,3,3,3,1});
        fullAssList.add(new int[] {-1,-1,4,-1,4,0,3,3,3,1});
        flexAssList.add(new int[] {4,4,1,3,3,3,1});
        fullAssList.add(new int[] {-1,-1,4,-1,4,1,3,3,3,1});
        flexAssList.add(new int[] {4,2,0,3,3,4,1});
        fullAssList.add(new int[] {-1,-1,4,-1,2,0,3,3,4,1});

        fullAssList.add(new int[] {0,-1,4,-1,2,0,3,3,4,1});
        fullAssList.add(new int[] {0,3,4,-1,2,0,3,3,4,1});
        fullAssList.add(new int[] {0,3,4,5,2,0,3,3,4,1});
        fullAssList.add(new int[] {0,3,4,5,2,5,3,3,4,1});
        fullAssList.add(new int[] {2,-1,4,-1,2,4,3,2,3,1});
        fullAssList.add(new int[] {2,0,4,-1,2,4,3,2,3,1});
        fullAssList.add(new int[] {2,0,4,0,2,4,3,2,3,1});

                metallochaperone = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
        //protein_1a0r = PDBIO.readFile("examples//1a0r_prepped.pdb");

        SimpleConfSpace confSpace = make1CC8Mutable();
        SimpleConfSpace flexSpace = confSpace.makeFlexibleCopy();

        Sequence wildType = confSpace.makeWildTypeSequence();
        RCs fullRCs = wildType.makeRCs(confSpace);
        RCs flexRCs = wildType.makeRCs(flexSpace);

        ConfAnalyzer flexMinimizer = makeMinimizerForConfSpace(flexSpace);
        ConfAnalyzer fullMinimizer = makeMinimizerForConfSpace(confSpace);

        ScorerFactory LBgScorerFactory = (emat) -> new PairwiseGScorer(emat);
        ScorerFactory LBhScorerFactory = (emat) -> new SHARKStarNodeScorer(emat, false);
        ScorerFactory UBgScorerFactory = (emat) -> new PairwiseRigidGScorer(emat);
        ScorerFactory UBhScorerFactory = (emat) -> new SHARKStarNodeScorer(emat, true);

        EnergyMatrix flexMinEmat = makeMinimizingEmatForConfSpace(flexSpace);
        EnergyMatrix flexRigidEmat = makeRigidEmatForConfSpace(flexSpace);
        EnergyMatrix fullMinEmat = makeMinimizingEmatForConfSpace(confSpace);
        EnergyMatrix fullRigidEmat = makeRigidEmatForConfSpace(confSpace);

        AStarScorer flexGScorerLB = LBgScorerFactory.make(flexMinEmat);
        AStarScorer flexHScorerLB = LBhScorerFactory.make(flexMinEmat);
        AStarScorer flexGScorerUB = UBgScorerFactory.make(flexRigidEmat);
        AStarScorer flexHScorerUB = UBhScorerFactory.make(flexRigidEmat);

        AStarScorer fullGScorerLB = LBgScorerFactory.make(fullMinEmat);
        AStarScorer fullHScorerLB = LBhScorerFactory.make(fullMinEmat);
        AStarScorer fullGScorerUB = UBgScorerFactory.make(fullRigidEmat);
        AStarScorer fullHScorerUB = UBhScorerFactory.make(fullRigidEmat);

        // Get info for flex confs
        ConfIndex flexIndex = new ConfIndex(flexSpace.getNumPos());
        for (int[] assignment : flexAssList){
            RCTuple tup = new RCTuple(assignment);
            tup.pasteToIndex(flexIndex);
            flexibleConfList.add(new TestConfInfo(
                    assignment,
                    flexGScorerLB.calc(flexIndex, flexRCs) + flexHScorerLB.calc(flexIndex,flexRCs),
                    flexGScorerUB.calc(flexIndex, flexRCs) + flexHScorerUB.calc(flexIndex,flexRCs),
                    flexMinimizer.analyze(new ConfSearch.ScoredConf(assignment, 0.0)).epmol.energy
            ));
        }

        ConfIndex fullIndex = new ConfIndex(confSpace.getNumPos());
        for (int[] assignment : fullAssList){
            RCTuple tup = new RCTuple(assignment);
            tup.pasteToIndex(fullIndex);
            fullConfList.add(new TestConfInfo(
                    assignment,
                    fullGScorerLB.calc(fullIndex, fullRCs) + fullHScorerLB.calc(fullIndex,fullRCs),
                    fullGScorerUB.calc(fullIndex, fullRCs) + fullHScorerUB.calc(fullIndex,fullRCs),
                    fullMinimizer.analyze(new ConfSearch.ScoredConf(assignment, 0.0)).epmol.energy
            ));
        }
    }

    /**
     * Creates a mutable confspace with one chain using 10 residues from 1CC8
     */
    private static SimpleConfSpace make1CC8Mutable(){
        Strand strand1 = new Strand.Builder(metallochaperone).setResidues("A2", "A20").build();
        strand1.flexibility.get("A2").setLibraryRotamers(Strand.WildType, "ARG").setContinuous();
        strand1.flexibility.get("A3").setLibraryRotamers(Strand.WildType, "ILE").setContinuous();
        strand1.flexibility.get("A4").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A5").setLibraryRotamers(Strand.WildType, "MET").setContinuous();
        strand1.flexibility.get("A6").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A7").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A8").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A9").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A10").setLibraryRotamers(Strand.WildType).setContinuous();
        strand1.flexibility.get("A11").setLibraryRotamers(Strand.WildType).setContinuous();

        return new SimpleConfSpace.Builder()
                .addStrands(strand1)
                .build();

    }

    /**
     * Create a minimizing ConfAnalyzer for a particular confspace
     * @param confSpace     The confSpace over which to compute energies
     * @return              A ConfAnalyzer object that minimizes
     */
    private static ConfAnalyzer makeMinimizerForConfSpace(SimpleConfSpace confSpace){
        Parallelism parallelism = Parallelism.makeCpu(4);
        ForcefieldParams ffparams = new ForcefieldParams();

        // Set up the minimizing ecalc
        EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();
        // Set up the minimizing confEcalc
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();

        // Set up the confAnalyzer
        return new ConfAnalyzer(confEcalcMinimized);
    }

    /**
     * Create a rigid ConfAnalyzer for a particular confspace
     * @param confSpace     The confSpace over which to compute energies
     * @return              A ConfAnalyzer object that minimizes
     */
    private static ConfAnalyzer makeRigidAnalyzerForConfSpace(SimpleConfSpace confSpace){
        // Set up requirements for energy calculation
        Parallelism parallelism = Parallelism.makeCpu(4);
        ForcefieldParams ffparams = new ForcefieldParams();

        // Set up the minimizing ecalc
        EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();

        // Set up the minimizing confEcalc
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();
        // Set up the rigid ecalc and confEcalc
        EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                .setIsMinimizing(false)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, ecalcRigid);

        // Set up the confAnalyzer
        return new ConfAnalyzer(confEcalcRigid);
    }

    private static EnergyMatrix makeRigidEmatForConfSpace(SimpleConfSpace confSpace){
        // Set up requirements for energy calculation
        Parallelism parallelism = Parallelism.makeCpu(4);
        ForcefieldParams ffparams = new ForcefieldParams();

        // Set up the minimizing ecalc
        EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();

        // Set up the minimizing confEcalc
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();
        // Set up the rigid ecalc and confEcalc
        EnergyCalculator ecalcRigid = new EnergyCalculator.SharedBuilder(ecalcMinimized)
                .setIsMinimizing(false)
                .build();
        ConfEnergyCalculator confEcalcRigid = new ConfEnergyCalculator(confEcalcMinimized, ecalcRigid);

        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalcRigid)
                .build()
                .calcEnergyMatrix();
        return emat;
    }

    private static EnergyMatrix makeMinimizingEmatForConfSpace(SimpleConfSpace confSpace) {
        // Set up requirements for energy calculation
        Parallelism parallelism = Parallelism.makeCpu(4);
        ForcefieldParams ffparams = new ForcefieldParams();

        // Set up the minimizing ecalc
        EnergyCalculator ecalcMinimized = new EnergyCalculator.Builder(confSpace, ffparams)
                .setParallelism(parallelism)
                .build();

        // Set up the minimizing confEcalc
        ConfEnergyCalculator confEcalcMinimized = new ConfEnergyCalculator.Builder(confSpace, ecalcMinimized)
                .setReferenceEnergies(new SimplerEnergyMatrixCalculator.Builder(confSpace, ecalcMinimized)
                        .build()
                        .calcReferenceEnergies()
                )
                .build();
        EnergyMatrix emat = new SimplerEnergyMatrixCalculator.Builder(confEcalcMinimized)
                .build()
                .calcEnergyMatrix();
        return emat;
    }

    @Test
    public void testMinimizeFlexible(){
        // Make the confspace (flexible only)
        SimpleConfSpace confSpace = make1CC8Mutable().makeFlexibleCopy();
        // Get the analyzer
        ConfAnalyzer analyzer = makeMinimizerForConfSpace(confSpace);

        for (int i=0; i < flexibleConfList.size(); i++){
            ConfAnalyzer.ConfAnalysis analysis = analyzer.analyze(new ConfSearch.ScoredConf(flexibleConfList.get(i).assignments, flexibleConfList.get(i).ELB));
            assertThat(analysis.epmol.energy, isAbsolutely(flexibleConfList.get(i).minE, 0.001));
        }

    }

    @Test
    public void testGetDOFs(){
        // Make the confspace (flexible only)
        SimpleConfSpace confSpace = make1CC8Mutable().makeFlexibleCopy();
        // Get the analyzer
        ConfAnalyzer analyzer = makeMinimizerForConfSpace(confSpace);

        flexibleConfList.stream().map( c -> analyzer.analyze(new ConfSearch.ScoredConf(c.assignments, c.ELB))).forEach( analysis -> {
            System.out.println(String.format("Original DOF bounds: %s", analysis.epmol.pmol.dofBounds));
            System.out.println(String.format("Minimized DOF values: %s", analysis.epmol.params.toString()));
        });
    }

    @Test
    public void testRCDOFs(){
        // Make the confspace (flexible only)
        SimpleConfSpace confSpace = make1CC8Mutable().makeFlexibleCopy();
        confSpace.positions.get(2).resConfs.stream().forEach( rc -> {
            System.out.println(String.format("Index %d, Rotamer index %d",rc.index, rc.rotamerIndex));
            rc.dofBounds.forEach( (k,v) -> System.out.println(String.format("DOF %s bounds: [%.3f, %.3f]",k,v[0], v[1])));
        });
    }
    @Test
    public void testAssignDOFsFromEPMOL(){
        // Make the confspace (flexible only)
        SimpleConfSpace confSpace = make1CC8Mutable().makeFlexibleCopy();
        // Get the analyzer
        ConfAnalyzer analyzer = makeMinimizerForConfSpace(confSpace);

        //Record minimized DOF values
        List<DoubleMatrix1D> minDOFs = new ArrayList<>();
        List<ObjectiveFunction.DofBounds> rigidDOFs = new ArrayList<>();

        //Record new rotamerAssignments
        List<RCTuple> mapTupList = new ArrayList<>();

        // Add new RCs from confspace
        for (TestConfInfo conf: flexibleConfList){
            ConfAnalyzer.ConfAnalysis analysis = analyzer.analyze((new ConfSearch.ScoredConf(conf.assignments, conf.ELB)));
            mapTupList.add(confSpace.addResidueConfFromMinimization(new RCTuple(analysis.assignments), analysis.epmol.pmol.dofs, analysis.epmol.params));
            minDOFs.add(analysis.epmol.params);
        }
        minDOFs.forEach(System.out::println);
        mapTupList.forEach(System.out::println);


        // For testing, make molecules from "rotamers" for minimized conformations
        for (RCTuple tup : mapTupList){
            ConfIndex index = new ConfIndex(confSpace.getNumPos());
            tup.pasteToIndex(index);
            rigidDOFs.add(confSpace.makeMolecule(index.makeConf()).dofBounds);
        }

        // Check to make sure that the DOFS for these new "rotamers" match the minmized DOF values
        for(int i=0; i<minDOFs.size(); i++){
            for (int j=0; j<minDOFs.get(i).size(); j++){
                assertThat(minDOFs.get(i).get(j), isAbsolutely(rigidDOFs.get(i).getBounds()[0].get(j)));
            }
        }

        // Check to make sure the energies look the same
        ConfAnalyzer rigidAnalyzer = makeRigidAnalyzerForConfSpace(confSpace);
        for (int i=0; i < flexibleConfList.size(); i++){
            ConfIndex index = new ConfIndex(confSpace.getNumPos());
            mapTupList.get(i).pasteToIndex(index);
            ConfAnalyzer.ConfAnalysis rigidAnalysis = rigidAnalyzer.analyze(index.makeConf());
            assertThat(rigidAnalysis.epmol.energy, isAbsolutely(flexibleConfList.get(i).minE, 0.001));
        }
    }

    @Test
    public void testAssignDOFsToFullConfSpace(){

        // Make the confspace
        SimpleConfSpace confSpace = make1CC8Mutable();
        SimpleConfSpace confSpaceUntouched = make1CC8Mutable();
        SimpleConfSpace flexCopy = confSpace.makeFlexibleCopy();
        // Get the analyzer
        ConfAnalyzer analyzer = makeMinimizerForConfSpace(flexCopy);

        // make the original, usual rigidEmat
        EnergyMatrix rigidEmat = makeRigidEmatForConfSpace(confSpaceUntouched);

        //Record minimized DOF values
        List<int[]> flexAssignments = new ArrayList<>();
        List<int[]> fullAssignments = new ArrayList<>();
        List<DoubleMatrix1D> minDOFs = new ArrayList<>();
        List<ObjectiveFunction.DofBounds> rigidDOFs = new ArrayList<>();

        //Record new rotamerAssignments
        List<RCTuple> mapTupList = new ArrayList<>();

        //Generate the array to map from the flex confspace to the full confspace
        int[] permArray = flexCopy.positions.stream()
                .mapToInt(confSpace.positions::indexOf)
                .toArray();

        // Score the flexible confspace leaves
        for (TestConfInfo conf : flexibleConfList){
            ConfAnalyzer.ConfAnalysis analysis = analyzer.analyze(new ConfSearch.ScoredConf(conf.assignments, conf.ELB));
            flexAssignments.add(analysis.assignments);
            minDOFs.add(analysis.epmol.params);

            // permute the flex leaf into the full confspace internal node
            int[] permutedAss = new int[confSpace.getNumPos()];
            Arrays.fill(permutedAss, -1);
            for (int j =0; j < analysis.assignments.length; j++){
                permutedAss[permArray[j]] = analysis.assignments[j];
            }
            fullAssignments.add(permutedAss);

            //Add the new RCs from the partial minimization into the full confspace
            mapTupList.add(confSpace.addResidueConfFromMinimization(new RCTuple(permutedAss), analysis.epmol.pmol.dofs,analysis.epmol.params));
        }


        // For testing, make molecules from "rotamers" for minimized conformations
        for (RCTuple tup : mapTupList){
            ConfIndex index = new ConfIndex(confSpace.getNumPos());
            tup.pasteToIndex(index);
            rigidDOFs.add(confSpace.makeMolecule(index.makeConf()).dofBounds);
        }

        // Check to make sure that the DOFS for these new "rotamers" match the minmized DOF values
        for(int i=0; i<minDOFs.size(); i++){
            for (int j=0; j<minDOFs.get(i).size(); j++){
                assertThat(minDOFs.get(i).get(j), isAbsolutely(rigidDOFs.get(i).getBounds()[0].get(j)));
            }
        }

        // Check to make sure the energies look the same
        ConfAnalyzer rigidAnalyzer = makeRigidAnalyzerForConfSpace(confSpace);
        for (int i=0; i < flexibleConfList.size(); i++){
            ConfIndex index = new ConfIndex(confSpace.getNumPos());
            mapTupList.get(i).pasteToIndex(index);
            ConfAnalyzer.ConfAnalysis rigidAnalysis = rigidAnalyzer.analyze(index.makeConf());
            assertThat(rigidAnalysis.epmol.energy, isAbsolutely(flexibleConfList.get(i).minE, 0.001));
        }
        // Check to make sure the energies look the same
        Sequence wildType = confSpace.makeWildTypeSequence();
        RCs seqRCs = wildType.makeRCs(confSpace);
        EnergyMatrix fancyRigidEmat = makeRigidEmatForConfSpace(confSpace);
        PairwiseGScorer betterGscorer = new PairwiseGScorer(fancyRigidEmat);
        for (int i=0; i < flexibleConfList.size(); i++){
            //update the index
            ConfIndex index = new ConfIndex(confSpace.getNumPos());
            mapTupList.get(i).pasteToIndex(index);
            assertThat(betterGscorer.calc(index,seqRCs), isAbsolutely(flexibleConfList.get(i).minE, 0.001));
        }

        // Test hscorers
        PairwiseGScorer worseGscorer = new PairwiseGScorer(rigidEmat);
        SHARKStarNodeScorer worseHScorer = new SHARKStarNodeScorer(rigidEmat, true);
        SHARKStarNodeScorer betterHScorer = new SHARKStarNodeScorer(fancyRigidEmat, true);
        for (int i=0; i < flexibleConfList.size(); i++){
            // First, score the original style
            //update the index
            ConfIndex index = new ConfIndex(seqRCs.getNumPos());
            for (int j =0; j< fullAssignments.get(i).length; j++){
                if(fullAssignments.get(i)[j] > 0)
                    index.assignInPlace(j,fullAssignments.get(i)[j]);
            }
            double oldgscore = worseGscorer.calc(index,seqRCs);
            double oldhscore = worseHScorer.calc(index,seqRCs);

            // Now, score the new style
            //update the index
            ConfIndex index_better = new ConfIndex(seqRCs.getNumPos());
            mapTupList.get(i).pasteToIndex(index_better);
            double newGScore = betterGscorer.calc(index_better, seqRCs);
            double newHScore = betterHScorer.calc(index_better, seqRCs);
            System.out.println(String.format("Old way: %.3f + %.3f = %.3f", oldgscore, oldhscore, oldgscore+oldhscore));
            System.out.println(String.format("New way: %.3f + %.3f = %.3f", newGScore, newHScore,newGScore+newHScore));
        }

        // Try one fullyAssigned
        System.out.println("Fully assigned:");
        for (int i=0; i < flexibleConfList.size(); i++){
            // First, score the original style
            //update the index
            ConfIndex index = new ConfIndex(seqRCs.getNumPos());
            for (int j =0; j< fullAssignments.get(i).length; j++){
                if(fullAssignments.get(i)[j] > 0)
                    index.assignInPlace(j,fullAssignments.get(i)[j]);
                else
                    index.assignInPlace(j,0);
            }
            double oldgscore = worseGscorer.calc(index,seqRCs);
            double oldhscore = worseHScorer.calc(index,seqRCs);

            // Now, score the new style
            //update the index
            ConfIndex index_better = new ConfIndex(seqRCs.getNumPos());
            for (int j =0; j< fullAssignments.get(i).length; j++){
                if(fullAssignments.get(i)[j] > 0)
                    index_better.assignInPlace(j,fullAssignments.get(i)[j]);
                else
                    index_better.assignInPlace(j,0);
            }
            mapTupList.get(i).pasteToIndex(index_better);
            double newGScore = betterGscorer.calc(index_better, seqRCs);
            double newHScore = betterHScorer.calc(index_better, seqRCs);
            System.out.println(String.format("Old way: %.3f + %.3f = %.3f", oldgscore, oldhscore, oldgscore+oldhscore));
            System.out.println(String.format("New way: %.3f + %.3f = %.3f", newGScore, newHScore,newGScore+newHScore));
            //assertThat(gscorer.calc(index,seqRCs), isAbsolutely(minConfEs.get(i), 0.001));
        }
    }

    private String RCtoString(SimpleConfSpace.ResidueConf rc){
        StringBuilder out = new StringBuilder(String.format("%s%d, %d: ", rc.type.letter, rc.index, rc.rotamerIndex));
        for (Map.Entry<String, double[]> entry : rc.dofBounds.entrySet()) {
            String k = entry.getKey();
            double[] v = entry.getValue();
            out.append((String.format("%s: [%.3f, %.3f]", k, v[0], v[1])));
        }
        return out.toString();
    }

    @Test
    public void testScoreCorrector(){
        // Make the confspace
        SimpleConfSpace confSpace = make1CC8Mutable();
        SimpleConfSpace confSpaceUntouched = make1CC8Mutable();
        SimpleConfSpace flexCopy = confSpace.makeFlexibleCopy();
        // Get the analyzer
        ConfAnalyzer analyzer = makeMinimizerForConfSpace(flexCopy);

        // make the original, usual rigidEmat
        EnergyMatrix rigidEmat = makeRigidEmatForConfSpace(confSpaceUntouched);

        List<int[]> fullAssignments = new ArrayList<>();
        //Corrections list
        ArrayList<TupEMapping> corrList = new ArrayList<>();

        //Generate the array to map from the flex confspace to the full confspace
        int[] permArray = flexCopy.positions.stream()
                .mapToInt(confSpace.positions::indexOf)
                .toArray();

        // Score the flexible confspace leaves
        for (TestConfInfo conf : flexibleConfList){
            ConfAnalyzer.ConfAnalysis analysis = analyzer.analyze(new ConfSearch.ScoredConf(conf.assignments, conf.ELB));

            // permute the flex leaf into the full confspace internal node
            int[] permutedAss = new int[confSpace.getNumPos()];
            Arrays.fill(permutedAss, -1);
            for (int j =0; j < analysis.assignments.length; j++){
                permutedAss[permArray[j]] = analysis.assignments[j];
            }
            fullAssignments.add(permutedAss);
            RCTuple originalTup = new RCTuple(permutedAss);

            //Add the new RCs from the partial minimization into the full confspace
            RCTuple mappedTup = confSpace.addResidueConfFromMinimization(originalTup, analysis.epmol.pmol.dofs,analysis.epmol.params);
            corrList.add(new TupEMapping(originalTup, mappedTup, analysis.epmol.energy - analysis.score));
        }

        // Make the new emat
        Sequence wildType = confSpace.makeWildTypeSequence();
        RCs seqRCs = wildType.makeRCs(confSpace);
        EnergyMatrix fancyRigidEmat = makeRigidEmatForConfSpace(confSpace);
        ScorerFactory gScorerFactory = (emat) -> new PairwiseRigidGScorer(emat);
        ScorerFactory hScorerFactory = (emat) -> new SHARKStarNodeScorer(emat, true);

        DependentScoreCorrector upperBoundCorrector = new DependentScoreCorrector(confSpace, MathTools.Optimizer.Minimize, fancyRigidEmat, gScorerFactory, hScorerFactory);
        upperBoundCorrector.insertAllCorrections(corrList);

        // Do testing
        for (TestConfInfo fullConf: fullConfList){
            RCTuple tupQuery = new RCTuple(fullConf.assignments);
            ConfIndex index = new ConfIndex(seqRCs.getNumPos());
            tupQuery.pasteToIndex(index);

            double correction = upperBoundCorrector.getCorrection(tupQuery,fullConf.EUB,index, seqRCs);
            assert(correction <= 0);
            System.out.println("------");
            System.out.println(String.format("Correction for %s -> %.3f", tupQuery.toString(), correction));
            System.out.println(String.format("Results in %.3f in [%.3f, %.3f] corrected -> [%.3f, %.3f]",
                    fullConf.minE,
                    fullConf.ELB, fullConf.EUB,
                    fullConf.ELB, fullConf.EUB+correction
            ));
        }

    }

    @Test
    public void testScoreCorrectorEarlyMake(){
        // Make the confspace
        SimpleConfSpace confSpace = make1CC8Mutable();
        SimpleConfSpace confSpaceUntouched = make1CC8Mutable();
        SimpleConfSpace flexCopy = confSpace.makeFlexibleCopy();
        // Get the analyzer
        ConfAnalyzer analyzer = makeMinimizerForConfSpace(flexCopy);

        // make the upperBoundCorrector early
        DependentScoreCorrector upperBoundCorrector = new DependentScoreCorrector(confSpace, MathTools.Optimizer.Minimize);

        List<int[]> fullAssignments = new ArrayList<>();
        //Corrections list
        ArrayList<TupEMapping> corrList = new ArrayList<>();

        //Generate the array to map from the flex confspace to the full confspace
        int[] permArray = flexCopy.positions.stream()
                .mapToInt(confSpace.positions::indexOf)
                .toArray();

        // Score the flexible confspace leaves
        for (TestConfInfo conf : flexibleConfList){
            ConfAnalyzer.ConfAnalysis analysis = analyzer.analyze(new ConfSearch.ScoredConf(conf.assignments, conf.ELB));

            // permute the flex leaf into the full confspace internal node
            int[] permutedAss = new int[confSpace.getNumPos()];
            Arrays.fill(permutedAss, -1);
            for (int j =0; j < analysis.assignments.length; j++){
                permutedAss[permArray[j]] = analysis.assignments[j];
            }
            fullAssignments.add(permutedAss);
            RCTuple originalTup = new RCTuple(permutedAss);

            //Add the new RCs from the partial minimization into the full confspace
            RCTuple mappedTup = confSpace.addResidueConfFromMinimization(originalTup, analysis.epmol.pmol.dofs,analysis.epmol.params);
            upperBoundCorrector.insertCorrection(new TupEMapping(originalTup, mappedTup, analysis.epmol.energy - analysis.score));
        }

        // Make the new emat
        Sequence wildType = confSpace.makeWildTypeSequence();
        RCs seqRCs = wildType.makeRCs(confSpace);
        EnergyMatrix fancyRigidEmat = makeRigidEmatForConfSpace(confSpace);
        ScorerFactory gScorerFactory = (emat) -> new PairwiseRigidGScorer(emat);
        ScorerFactory hScorerFactory = (emat) -> new SHARKStarNodeScorer(emat, true);

        // Set the stuff
        upperBoundCorrector.setMappedEmat(fancyRigidEmat);
        upperBoundCorrector.setGScorerFactory(gScorerFactory);
        upperBoundCorrector.setHScorerFactory(hScorerFactory);

        // Do testing
        for (TestConfInfo fullConf: fullConfList){
            RCTuple tupQuery = new RCTuple(fullConf.assignments);
            ConfIndex index = new ConfIndex(seqRCs.getNumPos());
            tupQuery.pasteToIndex(index);

            double correction = upperBoundCorrector.getCorrection(tupQuery,fullConf.EUB,index, seqRCs);
            assert(correction <= 0);
            System.out.println("------");
            System.out.println(String.format("Correction for %s -> %.3f", tupQuery.toString(), correction));
            System.out.println(String.format("Results in %.3f in [%.3f, %.3f] corrected -> [%.3f, %.3f]",
                    fullConf.minE,
                    fullConf.ELB, fullConf.EUB,
                    fullConf.ELB, fullConf.EUB+correction
            ));
        }

    }
}
