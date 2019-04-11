/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.ewakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.astar.seq.SeqAStarTree;
import edu.duke.cs.osprey.astar.seq.order.SequentialSeqAStarOrder;
import edu.duke.cs.osprey.astar.seq.scoring.NOPSeqAStarScorer;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.SimpleReferenceEnergies;
import edu.duke.cs.osprey.ematrix.SimplerEnergyMatrixCalculator;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.EnergyCalculator;
import edu.duke.cs.osprey.energy.forcefield.ForcefieldParams;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.pruning.PruningMatrix;
import edu.duke.cs.osprey.pruning.SimpleDEE;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.util.function.Function;
import java.util.stream.Collectors;

import static edu.duke.cs.osprey.tools.Log.formatBig;


/**
 * Implementation of the COMETS multi-state algorithm to predict protein sequence mutations that improve binding affinity.
 * {@cite Hallen2016 Mark A. Hallen, Bruce R. Donald, 2016.
 * Comets (Constrained Optimization of Multistate Energies by Tree Search): A provable and efficient protein design
 * algorithm to optimize binding affinity and specificity with respect to sequence.
 * Journal of Computational Biology, 23(5), 311-321.}.
 */
public class EwakstarDoer {

    public static class Results {
        private EWAKStarBBKStar bbkstar;
        public List<Sequence> sequences;
    }

    /**
     * A state for a multi-state design
     * (e.g., an unbound state, or a complex state)
     *
     * For COMETS, all states must share the same sequence space.
     */
    public static class State {

        public static class InitException extends RuntimeException {

            public InitException(State state, String name) {
                super(String.format("set %s for state %s before running", name, state.name));
            }
        }


        /** A name for the state */
        public final String name;

        /** The conformation space for the state that describes the flexibility */
        public final SimpleConfSpace confSpace;

        public FragmentEnergies fragmentEnergies;
        public ConfEnergyCalculator confEcalc;
        public ConfEnergyCalculator confRigidEcalc;
        public Function<RCs,ConfAStarTree> confTreeFactoryMin;
        public Function<RCs,ConfAStarTree> confTreeFactoryRigid;
        public EnergyMatrix ematRigid;
        public EnergyMatrix emat;

        public State(String name, SimpleConfSpace confSpace) {
            this.name = name;
            this.confSpace = confSpace;
        }

        /**
         * make sure the state is fully configured
         */
        public void checkConfig() {
            if (fragmentEnergies == null) {
                throw new InitException(this, "fragmentEnergies");
            }
            if (confEcalc == null) {
                throw new InitException(this, "confEcalc");
            }
            if (confTreeFactoryMin == null) {
                throw new InitException(this, "confTreeFactoryMin");
            }
            if (confTreeFactoryRigid == null) {
                throw new InitException(this, "confTreeFactoryRigid");
            }
        }

        public double getSingleEnergy(int pos, int rc) {
            return fragmentEnergies.getEnergy(pos, rc);
        }

        public double getPairEnergy(int pos1, int rc1, int pos2, int rc2) {
            return fragmentEnergies.getEnergy(pos1, rc1, pos2, rc2);
        }
    }

    /**
     * implements A* heuristic for partially-defined sequences
     * as described in COMETS paper, SI section B.2
     */
    private class SeqHScorer implements SeqAStarScorer {

        MathTools.Optimizer opt = MathTools.Optimizer.Minimize;

        // collect conf space positions from all states
        List<SimpleConfSpace.Position> allPositions = new ArrayList<>();

        // the states associated with each position
        State statesByAllPosition;

        SeqHScorer() {
                for (SimpleConfSpace.Position pos : state.confSpace.positions) {
                    allPositions.add(pos);
                    statesByAllPosition = state;
                }
        }

        @Override
        public double calc(SeqAStarNode.Assignments assignments) {

            // TODO: if inner loops are independent of assignments, pre-compute them somehow?

            double score = 0;
            // sum over all positions
            for (int i1=0; i1<allPositions.size(); i1++) {
                SimpleConfSpace.Position pos1 = allPositions.get(i1);

                // optimize over res types at pos1
                double bestPos1Energy = opt.initDouble();
                for (SeqSpace.ResType rt1 : getRTs(pos1, assignments)) {

                    // min over RCs at (pos1,rt1,state)
                    double bestRC1Energy = opt.initDouble();
                    for (SimpleConfSpace.ResidueConf rc1 : getRCs(pos1, rt1, state)) {

                        double rc1Energy = 0.0;

                        // singles
                        rc1Energy += state.getSingleEnergy(pos1.index, rc1.index);

                        // pairs
                        for (int i2=0; i2<pos1.index; i2++) {
                            SimpleConfSpace.Position pos2 = state.confSpace.positions.get(i2);

                            // min over RTs at pos2
                            double bestRT2Energy = opt.initDouble();
                            for (SeqSpace.ResType rt2 : getRTs(pos2, assignments)) {

                                // min over RCs at (pos2,rt2,state)
                                double bestRC2Energy = opt.initDouble();
                                for (SimpleConfSpace.ResidueConf rc2 : getRCs(pos2, rt2, state)) {

                                    double rc2Energy = state.getPairEnergy(pos1.index, rc1.index, pos2.index, rc2.index);

                                    bestRC2Energy = opt.opt(bestRC2Energy, rc2Energy);
                                }

                                bestRT2Energy = opt.opt(bestRT2Energy, bestRC2Energy);
                            }

                            rc1Energy += bestRT2Energy;
                        }

                        bestRC1Energy = opt.opt(bestRC1Energy, rc1Energy);
                    }

                    bestPos1Energy = opt.opt(bestPos1Energy, bestRC1Energy);
                }

                score += bestPos1Energy;
            }

            return score;
        }

        List<SeqSpace.ResType> getRTs(SimpleConfSpace.Position confPos, SeqAStarNode.Assignments assignments) {

            // TODO: pre-compute this somehow?

            // map the conf pos to a sequence pos
            SeqSpace.Position seqPos = seqSpace.getPosition(confPos.resNum);
            if (seqPos != null) {

                Integer assignedRT = assignments.getAssignment(seqPos.index);
                if (assignedRT != null) {
                    // use just the assigned res type
                    return Collections.singletonList(seqPos.resTypes.get(assignedRT));
                } else {
                    // use all the res types at the pos
                    return seqPos.resTypes;
                }

            } else {

                // immutable position, use all the res types (should just be one)
                assert (confPos.resTypes.size() == 1);

                // use the null value to signal there's no res type here
                return Collections.singletonList(null);
            }
        }

        List<SimpleConfSpace.ResidueConf> getRCs(SimpleConfSpace.Position pos, SeqSpace.ResType rt, State state) {
            // TODO: pre-compute this somehow?
            if (rt != null) {
                // mutable pos, grab the RCs that match the RT
                return pos.resConfs.stream()
                        .filter(rc -> rc.template.name.equals(rt.name))
                        .collect(Collectors.toList());
            } else {
                // immutable pos, use all the RCs
                return pos.resConfs;
            }
        }
    }

    /**
     * essentially, an iterative mini-GMEC-finder for a sequence and a state
     */
    private static class StateConfs {

        private static class Key {

            final Sequence sequence;
            final State state;

            Key(Sequence sequence, State state) {
                this.sequence = sequence;
                this.state = state;
            }

            @Override
            public int hashCode() {
                return HashCalculator.combineHashes(
                        sequence.hashCode(),
                        state.hashCode()
                );
            }

            @Override
            public boolean equals(Object other) {
                return other instanceof Key && equals((Key)other);
            }

            public boolean equals(Key other) {
                return this.sequence.equals(other.sequence)
                        && this.state == other.state;
            }
        }

        final State state;
        final Sequence sequence; // filtered to the state

        ConfAStarTree confTree = null;
        ConfSearch.ScoredConf lowestScoringConf = null;
        ConfSearch.EnergiedConf wtConf = null;
        boolean isWT = false;

        List<ConfSearch.ScoredConf> confs = new ArrayList<>();

        StateConfs(Sequence sequence, State state) {

            this.state = state;
            this.sequence = sequence;
            // make the conf tree
            RCs rcs = sequence.makeRCs(state.confSpace);
            confTree = state.confTreeFactoryMin.apply(rcs);
            setIsWildType(this.sequence.toString().equals(this.state.confSpace.makeWildTypeSequence().toString()));
        }

        void setIsWildType(boolean val){
            isWT = val;
        }

        void getLowestScoringConf() {

            // get the next batch of confs
            confs.clear();

            //if we have the wild-type sequence, we want it's minimized energy
            if(isWT){
                lowestScoringConf = confTree.nextConf();
                wtConf = state.confEcalc.calcEnergy(lowestScoringConf);
            }
            else {
                lowestScoringConf = confTree.nextConf();
            }

            confs.clear();
        }

        double getObjectiveLowerBound() {
            return lowestScoringConf.getScore();
        }

//        double getObjectiveUpperBound() {
//            if (gmec != null) {
//                return gmec.getEnergy();
//            }
//            return minEnergyConf.getEnergy();
//        }
    }


    /**
     * storage for conf trees at each sequence node
     * also tracks GMECs for each state
     */
    private class SeqConfs {

        final StateConfs stateConfs;

        SeqConfs(SeqAStarNode seqNode) {

                Sequence sequence = seqNode.makeSequence(state.confSpace.seqSpace)
                        .filter(state.confSpace.seqSpace);

                // get state confs from the global cache if possible, otherwise make a new one
                StateConfs.Key key = new StateConfs.Key(sequence, state);
                StateConfs collectedConfs = stateConfsCache.get(key);
                if (collectedConfs == null) {
                    collectedConfs = new StateConfs(sequence, state);
                    stateConfsCache.put(key, collectedConfs);
                }

            stateConfs = collectedConfs;
        }

        /**
         * implements A* heuristic for fully-defined sequences
         * as described in COMETS paper, SI section B.1
         *
         * returns the new score for the seqeunce node
         *
         * also flags that GMECs are found, when applicable
         */
        public double refineBounds() {

            // refine the GMEC bounds for the state
            this.stateConfs.getLowestScoringConf();

            double score = 0;

            StateConfs stateConfs = this.stateConfs;
            score += stateConfs.getObjectiveLowerBound(); //get the score for the current conformation for the current state

            return score;

        }
    }

    public class SequenceInfo {

        public final Sequence sequence;
        public final ConfSearch.ScoredConf lowestConf;
        public final double pfUB;
        public final boolean isWT;
        public boolean canPrune = false;

        public SequenceInfo(SeqAStarNode node, SeqConfs confs, double pfUB) {
            this.sequence = node.makeSequence(seqSpace);
            this.isWT = confs.stateConfs.isWT;
            this.lowestConf = confs.stateConfs.lowestScoringConf;
            this.pfUB = pfUB;

        }
        public void setCanPrune(boolean val){
            this.canPrune = val;
        }
    }


    public static class Builder {

        public static class InitException extends RuntimeException {

            public InitException(String mutType) {
                super(String.format("The value entered, %s,  is not a valid mutableType setting. Please set it to max, all, or exact.", mutType));
            }
        }

        private State state;

        /**
         * The energy window is actually necessary for COMETS to finish in a reasonable
         * amount of time in some edge cases. If the constraints restrict the sequence
         * space to fewer than the desired number of best sequences, without an energy window,
         * COMETS would enumerate every sequence while trying to find the desired number of
         * sequences. This would effectively make COMETS run in time that is linear in the
         * number of possible sequences, which could be very slow.
         *
         * The window size effectively adds an additional constraint that the difference between
         * the objective function value and the objective function value of the best sequence
         * must be below this value.
         */
        private double eW = 10.0;

        int numCPUs = 4;
        /**
         * this value can be max (mutations are <= the specified value), exact (exactly this many mutations),
         * or all (all positions set-up to be mutable can mutate).
         */
        private String mutableType = "all";

        /** The maximum number of simultaneous residue mutations to consider for each sequence mutant */
        private int numMutable = (int) Double.POSITIVE_INFINITY;

        private boolean printToConsole = true;

        private boolean seqFilterOnly = false;

        private Double numEWAKStarSeqs = Double.POSITIVE_INFINITY;
        private int orderOfMag = 5;
        private double pfEw = 1.0;
        private int numPfConfs = 500;
        private int numTopOverallSeqs = 10;
        private boolean useWtBenchmark = false;
        private double epsilon = 0.68;
        private int smaNodes;
        private boolean useSMA = false;
        private boolean printPDBs = false;

        /** File to which to log sequences as they are found */
        private File logFile = null;

        public Builder setEpsilon(double val){
            epsilon = val;
            return this;
        }

        public Builder setNumEWAKStarSeqs(Double val){
            numEWAKStarSeqs = val;
            return this;
        }

        public Builder setUseWtBenchmark(boolean val){
            useWtBenchmark = val;
            return this;
        }

        public Builder setOrderOfMag(int val){
            orderOfMag = val;
            return this;
        }

        public Builder setPfEw(double val){
            pfEw = val;
            return this;
        }

        public Builder setNumPfConfs(int val){
            numPfConfs = val;
            return this;
        }

        public Builder setNumTopOverallSeqs(int val){
            numTopOverallSeqs = val;
            return this;
        }

        public Builder setSeqFilterOnly(boolean val){
            seqFilterOnly = val;
            return this;
        }

        public Builder addState(State s){
            state = s;
            return this;
        }

        public Builder setNumMutable(int val) {
            numMutable = val;
            return this;
        }

        public Builder setEw(double val){
            eW = val;
            return this;
        }

        public Builder setMutableType(String type){
            if(!(type.equals("max")||type.equals("all")||type.equals("exact")))
                throw new Builder.InitException(type);
            mutableType = type;
            return this;
        }

        public Builder setNumCpus(int val){
            numCPUs = val;
            return this;
        }

        public Builder setPrintToConsole(boolean val) {
            printToConsole = val;
            return this;
        }

        public Builder setLogFile(File val) {
            logFile = val;
            return this;
        }

        public Builder setupSMA(boolean val, int num){
            useSMA = val;
            smaNodes = num;
            return this;
        }

        public Builder setPrintPDBs(boolean val){
            printPDBs = val;
            return this;
        }

        public EwakstarDoer build() {
            return new EwakstarDoer(printPDBs,
                    useSMA,
                    smaNodes,
                    state,
                    eW,
                    mutableType,
                    numMutable,
                    numCPUs,
                    printToConsole,
                    seqFilterOnly,
                    logFile,
                    numEWAKStarSeqs,
                    useWtBenchmark,
                    orderOfMag,
                    pfEw,
                    numPfConfs,
                    epsilon,
                    numTopOverallSeqs
                    );
        }
    }


    public final Double numEWAKStarSeqs;
    public final int numCPUs;
    public final double eW;
    public final boolean printToConsole;
    public final boolean seqFilterOnly;
    public final File logFile;
    public final int orderOfMag;
    public final double pfEw;
    public final int numPfConfs;
    public final int numTopOverallSeqs;
    public final boolean useWtBenchmark;
    public final double epsilon;

    public final State state;
    public final SeqSpace seqSpace;
    public String mutableType;
    public int numMutable;
    public boolean useExact;
    public EwakstarDoer ewakstarDoerP;
    public EwakstarDoer ewakstarDoerL;
    public EwakstarDoer ewakstarDoerPL;
    public final boolean useSMA;
    public final int smaNodes;
    public final boolean printPDBs;
    public List<SequenceInfo> infos = new ArrayList<>();

    private final Map<StateConfs.Key,StateConfs> stateConfsCache = new HashMap<>();

    private EwakstarDoer(boolean printPDBs, boolean useSMA, int smaNodes, State state, double eW, String mutableType, int numMutable, int numCPUs, boolean printToConsole, boolean seqFilterOnly, File logFile, Double numEWAKStarSeqs, boolean useWtBenchmark, int orderOfMag, double pfEw, int numPfConfs, double epsilon, int numTopOverallSeqs) {

        if(mutableType.equals("exact") || mutableType.equals("max")){
            this.numMutable = numMutable;
        } else if(mutableType.equals("all")){
            this.numMutable = state.confSpace.mutablePositions.size();
        }

        if(mutableType.equals("max")){
            this.useExact = false;
        }

        this.printPDBs = printPDBs;
        this.useSMA = useSMA;
        this.smaNodes = smaNodes;
        this.epsilon = epsilon;
        this.numEWAKStarSeqs = numEWAKStarSeqs;
        this.orderOfMag = orderOfMag;
        this.pfEw = pfEw;
        this.numPfConfs = numPfConfs;
        this.numTopOverallSeqs = numTopOverallSeqs;
        this.useWtBenchmark = useWtBenchmark;
        this.seqFilterOnly = seqFilterOnly;
        this.numCPUs = numCPUs;
        this.mutableType = mutableType;
        this.eW = eW;
        this.numMutable = numMutable;
        this.printToConsole = printToConsole;
        this.logFile = logFile;

        this.state = state;

        seqSpace = state.confSpace.seqSpace;

        log("sequence space has %s sequences\n%s", formatBig(new RTs(seqSpace).getNumSequences()), seqSpace);
    }

    /**
     * find the best sequences as ranked by the objective function
     *
     * searches all sequences within the objective window
     */
    public Set<Sequence> run(State PL) {

        ewakstarDoerPL = this;
        // reset any previous state
        stateConfsCache.clear();

        //Start timer
        long startEWAKStarTime = System.currentTimeMillis();

        System.out.println("Performing EWAK*");

        //calculate the size of the system
        BigDecimal totalConfs = BigDecimal.ONE;
        for(int i=0; i<PL.confSpace.positions.size(); i++) {
            BigDecimal newConfs = new BigDecimal (PL.confSpace.positions.get(i).resConfs.size());
            totalConfs = totalConfs.multiply(newConfs);
        }


        state.checkConfig();

        List<SequenceInfo> bestPLSeqs = extractPLSeqsByLB(numEWAKStarSeqs, orderOfMag, PL);

        //stop timer for energy matrix calculations
        long intermediateStopTime = System.currentTimeMillis();
        long timeTakenSoFar = intermediateStopTime - startEWAKStarTime;
        //limit sequence spaces for unbound states to those found in the bound complex above
        setupUnboundStates();
        long intermediateStartTime = System.currentTimeMillis();

        Set<Sequence> filteredSeqsP = new HashSet<>();
        Set<Sequence> filteredSeqsL = new HashSet<>();
        Set<String> filteredSeqsPL = new HashSet<>();

        for(SequenceInfo si : bestPLSeqs){
            filteredSeqsP.add(si.sequence.filter(ewakstarDoerP.seqSpace));
            filteredSeqsL.add(si.sequence.filter(ewakstarDoerL.seqSpace));
            filteredSeqsPL.add(si.sequence.toString());
        }

        EwakstarLimitedSequenceTrie elstP = new EwakstarLimitedSequenceTrie(ewakstarDoerP.seqSpace);
        if(filteredSeqsP.size()!=1) {
            for (Sequence s : filteredSeqsP) {
                elstP.addSeq(s.toString());
            }
        }

        EwakstarLimitedSequenceTrie elstL = new EwakstarLimitedSequenceTrie(ewakstarDoerL.seqSpace);
        if(filteredSeqsL.size()!=1) {
            for (Sequence s : filteredSeqsL) {
                elstL.addSeq(s.toString());
            }
        }

        List<SequenceInfo> bestPSeqs = null;
        List<SequenceInfo> bestLSeqs = null;
        stateConfsCache.clear();
        if(ewakstarDoerP.state.confSpace.mutablePositions.size()!=0)
            bestPSeqs = ewakstarDoerP.extractUnboundSeqsByLB(orderOfMag, filteredSeqsP, elstP);
        stateConfsCache.clear();
        if(ewakstarDoerL.state.confSpace.mutablePositions.size()!=0)
            bestLSeqs = ewakstarDoerL.extractUnboundSeqsByLB(orderOfMag, filteredSeqsL, elstL);


        Set<Sequence> newPSeqs = new HashSet<>();
        Set<Sequence> newLSeqs = new HashSet<>();

        if (bestPSeqs != null){
            for(SequenceInfo si : bestPSeqs){
                newPSeqs.add(si.sequence);
            }
        }
        if (bestLSeqs != null){
            for(SequenceInfo si : bestLSeqs) {
                newLSeqs.add(si.sequence);
            }
        }

        Set<Sequence> fullSeqs = combineSeqs(newPSeqs, newLSeqs, filteredSeqsPL);

        long stopEWAKStarTime;

        log("The original conformation space size is: %6.3e", totalConfs);
        if (!seqFilterOnly) {
            EwakstarLimitedSequenceTrie elstPL = new EwakstarLimitedSequenceTrie(ewakstarDoerPL.seqSpace);
            if(fullSeqs.size()!=0) {
                for (Sequence s : fullSeqs) {
                    elstPL.addSeq(s.toString());
                }
            }
            runEWAKStarBBKStar(elstPL, ewakstarDoerPL.state, ewakstarDoerP.state, ewakstarDoerL.state);
        }
        if(seqFilterOnly) {
            writeSeqsToFile(fullSeqs);
            System.out.println("Number of sequences filtered down to "+ fullSeqs.size()+" from "+formatBig(new RTs(seqSpace).getNumSequences()));
            stopEWAKStarTime = System.currentTimeMillis()-intermediateStartTime+timeTakenSoFar;
            System.out.println("Total OSPREY/EWAK* time (not including energy matrix time): "+(String.format("%d sec",
                    TimeUnit.MILLISECONDS.toSeconds(stopEWAKStarTime))));
            return fullSeqs;
        } else
            System.out.println("Number of sequences filtered down to "+ fullSeqs.size()+" from "+formatBig(new RTs(seqSpace).getNumSequences()));
        stopEWAKStarTime = System.currentTimeMillis()-intermediateStartTime+timeTakenSoFar;
        System.out.println("Total OSPREY/EWAK* time (not including energy matrix time): "+(String.format("%d sec",
                TimeUnit.MILLISECONDS.toSeconds(stopEWAKStarTime))));
        return null;

    }

    private void runEWAKStarBBKStar(EwakstarLimitedSequenceTrie plTrie, State PL, State P, State L){

        // run BBK*
        EWAKStar.Settings ewakstarSettings = new EWAKStar.Settings.Builder()
                .setEw(pfEw)
                .setUseExact(useExact)
                .setNumMutations(numMutable)
                .setEpsilon(epsilon)
                .setMaxNumConfs(numPfConfs)
                .setNumTopOverallSeqs(numTopOverallSeqs)
                .addScoreConsoleWriter()
                .setWTBenchmark(useWtBenchmark)
                .setPrintPDBs(printPDBs)
                .addScoreFileWriter(new File("ewakStar.results.txt"))
                .build();

        EWAKStarBBKStar.Settings bbkstarSettings = new EWAKStarBBKStar.Settings.Builder()
                .setNumBestSequences(numTopOverallSeqs)
                .setAllowedSeqs(plTrie).build();

        Results results = new Results();
        results.bbkstar = new EWAKStarBBKStar( P, L, PL, ewakstarSettings, bbkstarSettings, null);

        for (EWAKStarBBKStar.ConfSpaceInfo info : results.bbkstar.confSpaceInfos()) {

            // how should we define energies of conformations?
            if (info.id.equals("protein")) {
                info.confEcalcMinimized = P.confEcalc;
                info.confTreeFactoryMinimized = P.confTreeFactoryMin;
                info.confTreeFactoryRigid = P.confTreeFactoryRigid;
            } else if (info.id.equals("ligand")) {
                info.confEcalcMinimized = L.confEcalc;
                info.confTreeFactoryMinimized = L.confTreeFactoryMin;
                info.confTreeFactoryRigid = L.confTreeFactoryRigid;
            } else {
                info.confEcalcMinimized = PL.confEcalc;
                info.confTreeFactoryMinimized = PL.confTreeFactoryMin;
                info.confTreeFactoryRigid = PL.confTreeFactoryRigid;
            }

        }

        results.sequences = results.bbkstar.run();

    }

    private void writeSeqsToFile(Set<Sequence> seqs){

        File newFile = new File("ewakStar.filteredSeqs.txt");

        boolean append;
        boolean started = false;

        for (Sequence s: seqs){
            String curSeq = s.toString();
            if (!started) {
                started = true;
                append = false;
            }
            else{
                append = true;
            }
            try (FileWriter out = new FileWriter(newFile, append)) {
                out.write(curSeq);
                out.write("\n");
            } catch (IOException ex) {
                System.err.println("writing to file failed: " + newFile);
                ex.printStackTrace(System.err);
                System.err.println(curSeq);
            }
        }
    }

    public Set<Sequence> combineSeqs(Set<Sequence> P, Set<Sequence> L, Set<String> PL){

        List<String> resTypes;
        Set<Sequence> newPL = new HashSet<>();
        BigInteger totalConfs = BigInteger.ZERO;

        if(P.size()!=0 && L.size()!=0) {
            for (Sequence sP : P) {
                for (Sequence sL : L) {
                    resTypes = new ArrayList<>();
                    for (String p : ewakstarDoerPL.seqSpace.getResNums()) {
                        if (ewakstarDoerP.seqSpace.getResNums().contains(p)) {
                            resTypes.add(sP.get(p).name);
                        } else
                            resTypes.add(sL.get(p).name);
                    }
                    Sequence newSeq = ewakstarDoerPL.seqSpace.makeSequence(resTypes);
                    if (PL.contains(newSeq.toString())) {
                        newPL.add(newSeq);
                        totalConfs = totalConfs.add(getNumConfsForSeq(newSeq, ewakstarDoerPL.state.confSpace));
                    }
                }
            }
        } else if(P.size()!=0) {
            for (Sequence sP : P) {
                resTypes = new ArrayList<>();
                for (String p : ewakstarDoerPL.seqSpace.getResNums()) {
                    if (ewakstarDoerP.seqSpace.getResNums().contains(p)) {
                        resTypes.add(sP.get(p).name);
                    }
                }
                Sequence newSeq = ewakstarDoerPL.seqSpace.makeSequence(resTypes);
                if (PL.contains(newSeq.toString())) {
                    newPL.add(newSeq);
                    totalConfs = totalConfs.add(getNumConfsForSeq(newSeq, ewakstarDoerPL.state.confSpace));
                }

            }
        } else if(L.size()!=0) {
            for (Sequence sL : L) {
                resTypes = new ArrayList<>();
                for (String p : ewakstarDoerPL.seqSpace.getResNums()) {
                    if (ewakstarDoerL.seqSpace.getResNums().contains(p)) {
                        resTypes.add(sL.get(p).name);
                    }
                }
                Sequence newSeq = ewakstarDoerPL.seqSpace.makeSequence(resTypes);
                if (PL.contains(newSeq.toString())) {
                    newPL.add(newSeq);
                    totalConfs = totalConfs.add(getNumConfsForSeq(newSeq, ewakstarDoerPL.state.confSpace));
                }

            }
        }

        log("The remaining conformation space size is: %6.3e", totalConfs.doubleValue());
        return newPL;
    }

    public void setupUnboundStates(){

        Strand protein = state.confSpace.strands.get(0);
        Strand ligand = state.confSpace.strands.get(1);

        EwakstarDoer.State P = new EwakstarDoer.State(
                "P",
                new SimpleConfSpace.Builder()
                        .addStrand(protein)
                        .build()
        );

        EwakstarDoer.State L = new EwakstarDoer.State(
                "L",
                new SimpleConfSpace.Builder()
                        .addStrand(ligand)
                        .build()
        );

        ewakstarDoerP = new EwakstarDoer.Builder()
                .addState(P)
                .setEw(eW)
                .setNumCpus(numCPUs)
                .setMutableType("all")
                .build();

        ewakstarDoerL = new EwakstarDoer.Builder()
                .addState(L)
                .setEw(eW)
                .setNumCpus(numCPUs)
                .setMutableType("all")
                .build();

        ForcefieldParams ffparams = new ForcefieldParams();

        EnergyCalculator ecalcP = new EnergyCalculator.Builder(P.confSpace, ffparams)
                .setParallelism(Parallelism.makeCpu(ewakstarDoerP.numCPUs))
                .build();
        EnergyCalculator rigidEcalcP = new EnergyCalculator.Builder(P.confSpace, ffparams)
                .setIsMinimizing(false)
                .setParallelism(Parallelism.makeCpu(ewakstarDoerP.numCPUs))
                .build();

        // what are conformation energies?
        SimpleReferenceEnergies erefP = new SimplerEnergyMatrixCalculator.Builder(P.confSpace, ecalcP)
                .build()
                .calcReferenceEnergies();
        SimpleReferenceEnergies rigidErefP = new SimplerEnergyMatrixCalculator.Builder(P.confSpace, rigidEcalcP)
                .build()
                .calcReferenceEnergies();

        P.confEcalc = new ConfEnergyCalculator.Builder(P.confSpace, ecalcP)
                .setReferenceEnergies(erefP)
                .build();
        P.confRigidEcalc = new ConfEnergyCalculator.Builder(P.confSpace, rigidEcalcP)
                .setReferenceEnergies(rigidErefP)
                .build();

        // calc the energy matrix
        P.emat = new SimplerEnergyMatrixCalculator.Builder(P.confEcalc)
                .setCacheFile(new File(String.format("ewakstar.%s.emat", P.name)))
                .build()
                .calcEnergyMatrix();
        P.fragmentEnergies = P.emat;
        P.ematRigid = new SimplerEnergyMatrixCalculator.Builder(P.confRigidEcalc)
                .setCacheFile(new File(String.format("ewakstar.%s.ematRigid", P.name)))
                .build()
                .calcEnergyMatrix();


        // make the conf tree factory

        if (useSMA) {
            P.confTreeFactoryMin = (rcs) -> new ConfAStarTree.Builder(P.emat, rcs)
                    .setMaxNumNodes(smaNodes)
                    .setTraditional()
                    .build();

            P.confTreeFactoryRigid = (rcs) -> new ConfAStarTree.Builder(P.ematRigid, rcs)
                    .setMaxNumNodes(smaNodes)
                    .setTraditional()
                    .build();
        } else {
            P.confTreeFactoryMin = (rcs) -> new ConfAStarTree.Builder(P.emat, rcs)
                    .setTraditional()
                    .build();

            P.confTreeFactoryRigid = (rcs) -> new ConfAStarTree.Builder(P.ematRigid, rcs)
                    .setTraditional()
                    .build();
        }



        //do all of this for ligand also
        EnergyCalculator ecalcL = new EnergyCalculator.Builder(L.confSpace, ffparams)
                .setParallelism(Parallelism.makeCpu(ewakstarDoerL.numCPUs))
                .build();
        EnergyCalculator rigidEcalcL = new EnergyCalculator.Builder(L.confSpace, ffparams)
                .setIsMinimizing(false)
                .setParallelism(Parallelism.makeCpu(ewakstarDoerL.numCPUs))
                .build();

        // what are conformation energies?
        SimpleReferenceEnergies erefL = new SimplerEnergyMatrixCalculator.Builder(L.confSpace, ecalcL)
                .build()
                .calcReferenceEnergies();
        SimpleReferenceEnergies rigidErefL = new SimplerEnergyMatrixCalculator.Builder(L.confSpace, rigidEcalcL)
                .build()
                .calcReferenceEnergies();

        L.confEcalc = new ConfEnergyCalculator.Builder(L.confSpace, ecalcL)
                .setReferenceEnergies(erefL)
                .build();
        L.confRigidEcalc = new ConfEnergyCalculator.Builder(L.confSpace, rigidEcalcL)
                .setReferenceEnergies(rigidErefL)
                .build();

        // calc the energy matrix
        L.emat = new SimplerEnergyMatrixCalculator.Builder(L.confEcalc)
                .setCacheFile(new File(String.format("ewakstar.%s.emat", L.name)))
                .build()
                .calcEnergyMatrix();
        L.fragmentEnergies = L.emat;
        L.ematRigid = new SimplerEnergyMatrixCalculator.Builder(L.confRigidEcalc)
                .setCacheFile(new File(String.format("ewakstar.%s.ematRigid", L.name)))
                .build()
                .calcEnergyMatrix();

        // make the conf tree factory
        if (useSMA) {
            L.confTreeFactoryMin = (rcs) -> new ConfAStarTree.Builder(L.emat, rcs)
                    .setMaxNumNodes(smaNodes)
                    .setTraditional()
                    .build();
            L.confTreeFactoryRigid = (rcs) -> new ConfAStarTree.Builder(L.ematRigid, rcs)
                    .setMaxNumNodes(smaNodes)
                    .setTraditional()
                    .build();
        } else {
            L.confTreeFactoryMin = (rcs) -> new ConfAStarTree.Builder(L.emat, rcs)
                    .setTraditional()
                    .build();
            L.confTreeFactoryRigid = (rcs) -> new ConfAStarTree.Builder(L.ematRigid, rcs)
                    .setTraditional()
                    .build();
        }




    }

    public List<SequenceInfo> extractPLSeqsByLB(double numSequences, int orderMag, State PL) {

        // start the A* search over sequences
        SeqAStarTree seqTree = new SeqAStarTree.Builder(new RTs(seqSpace))
                .setHeuristics(
                        new SequentialSeqAStarOrder(),
                        new NOPSeqAStarScorer(),
                        new SeqHScorer()
                )
                .setMutableType(mutableType)
                .setNumMutable(numMutable)
                .build();

        double pfUB;
        double wtPfLB = Double.POSITIVE_INFINITY;
        boolean didSeqMax = false;

        double seqLB = Double.NEGATIVE_INFINITY;
        Boolean wtFound = false;
        double wtEnergy = 0.0;
        Boolean didEW = false;
        int wtSpot = 0;

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        for(int seqNum=0; seqNum<numSequences; seqNum++) {

            // get the next sequence from the tree
            SeqAStarNode node = seqTree.nextLeafNode();
            if (node == null) {
                break;
            }

            SeqConfs confs = (SeqConfs) node.getData();
            if (confs == null) {

                seqLB = node.getScore();

                log("Estimated Sequence: %s   Sequence minimized lower bound: %12.6f",
                        node.makeSequence(seqSpace),
                        seqLB
                );

                // don't have them yet, make them
                confs = new SeqConfs(node);
                double curScore = confs.refineBounds();
                node.setData(confs);

            }

            if(confs.stateConfs.isWT){
                wtFound = true;
                wtSpot = seqNum;
                wtEnergy = confs.stateConfs.wtConf.getEnergy();
                wtPfLB = Math.log10(bc.calc(wtEnergy).doubleValue());
                System.out.println("Found wild-type sequence in complex after " + seqNum + " sequences.");
                System.out.println("Finding all sequences within " + eW + " kcal of WT minimized energy: "
                        + wtEnergy + " or until " + numSequences + " sequences are enumerated.");
            }

            Sequence newSequence = node.makeSequence(PL.confSpace.seqSpace);
            pfUB = Math.log10(bc.calc(confs.stateConfs.lowestScoringConf.getScore()).multiply(new BigDecimal (getNumConfsForSeq(newSequence, PL.confSpace))).doubleValue());
            SequenceInfo info = new SequenceInfo(node, confs, pfUB);
            if (!wtFound) {
                infos.add(info);
                reportSequence(infos.size() == 1, info);
            } else {
                if(seqLB > wtEnergy+eW) { //if we have surpassed the energy window we want to quit
                    if(infos.size() >= numSequences)
                        didSeqMax = true;
                    didEW = true;
                    break;
                } else if (infos.size() >= numSequences) { //if we have surpassed the maximum number of sequence we want to find, we want to quit
                    didSeqMax = true;
                    break;
                } else if(confs.stateConfs.lowestScoringConf.getScore() < wtEnergy+eW && pfUB > wtPfLB - orderMag) { //otherwise, add the sequence if it's conformation minimized lower bound is within the energy window and it's partition function upper bound is within orderMag of the wt pf lower bound
                    infos.add(info);
                    reportSequence(infos.size() == 1, info);
                }
            }
        }

        List<SequenceInfo> badPfs = new ArrayList<>();
        for(int i=0; i<wtSpot; i++) {
            if (infos.get(i).pfUB < wtPfLB - orderMag || infos.get(i).lowestConf.getScore() > wtEnergy+eW) {
                badPfs.add(infos.get(i));
            }
        }

        infos.removeAll(badPfs);

        log("");
        if (infos.isEmpty()) {
            log("EWAK* didn't find any sequences within the window that satisfy all the constraints.");
        } else {
            log("EWAK* found the best %d sequences within %d orders of magnitude of the wild-type sequence partition function.", infos.size(), orderMag);
        }

        if(didEW){
            log("Found all sequences within %1.0f kcal/mol of the wild-type sequence.", eW);
        } else {
            log("The energy window of %1.0f kcal/mol was not fully enumerated.", eW);
        }

        if(didSeqMax){
            log("The sequence number max of %d was enumerated.", numSequences);
        } else if(didEW){
            log("The energy window was completed before enumerating the requested number of sequences.");
        } else if(!didSeqMax && !didEW){
            log("Exhaustively enumerated all sequence possibilities!");
        }

        return infos;

    }

    public List<SequenceInfo> extractUnboundSeqsByLB(int orderMag, Set<Sequence> seqsToUse, EwakstarLimitedSequenceTrie elst) {

        // start the A* search over sequences
        SeqAStarTree seqTree = new SeqAStarTree.Builder(new RTs(state.confSpace.seqSpace))
                .setSeqTrie(elst)
                .setHeuristics(
                        new SequentialSeqAStarOrder(),
                        new NOPSeqAStarScorer(),
                        new SeqHScorer()
                )
                .setMutableType(mutableType)
                .setNumMutable(numMutable)
                .build();

        List<SequenceInfo> infos = new ArrayList<>();
        double pfUB;
        double wtPfLB = Double.POSITIVE_INFINITY;
        boolean didSeqMax = false;

        double seqLB = Double.NEGATIVE_INFINITY;
        Boolean wtFound = false;
        double wtEnergy = 0.0;
        Boolean didEW = false;
        int wtSpot = 0;

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        int seqNum = 0;
        while(infos.size()<seqsToUse.size()) {

            // get the next sequence from the tree
            SeqAStarNode node = seqTree.nextLeafNode();
            if (node == null) {
                break;
            }

            SeqConfs confs = (SeqConfs) node.getData();
            if (confs == null) {

                seqLB = node.getScore();
                log("Estimated Sequence: %s   Sequence minimized lower bound: %12.6f",
                        node.makeSequence(state.confSpace.seqSpace),
                        seqLB
                );

                // don't have them yet, make them
                confs = new SeqConfs(node);
                double curScore = confs.refineBounds();
                node.setData(confs);

            }

            if(confs.stateConfs.isWT){
                wtFound = true;
                wtSpot = seqNum;
                wtEnergy = confs.stateConfs.wtConf.getEnergy();
                wtPfLB = Math.log10(bc.calc(wtEnergy).doubleValue());
                System.out.println("Found wild-type sequence in complex after " + seqNum + " sequences.");
                System.out.println("Finding all sequences within " + eW + " kcal of WT minimized energy: "
                        + wtEnergy + " or until " + seqsToUse.size() + " sequences are enumerated.");
            }

            Sequence newSequence = node.makeSequence(state.confSpace.seqSpace);
            pfUB = Math.log10(bc.calc(confs.stateConfs.lowestScoringConf.getScore()).multiply(new BigDecimal(getNumConfsForSeq(newSequence, state.confSpace))).doubleValue());
            SequenceInfo info = new SequenceInfo(node, confs, pfUB);
            if (!wtFound) {
                infos.add(info);
                seqNum++;
                reportSequence(infos.size() == 1, info);
            } else {
                if (seqLB > wtEnergy + eW) { //if we have surpassed the energy window we want to quit
                    if (infos.size() >= seqsToUse.size())
                        didSeqMax = true;
                    didEW = true;
                    break;
                } else if (infos.size() >= seqsToUse.size()) { //if we have surpassed the maximum number of sequence we want to find, we want to quit
                    didSeqMax = true;
                    break;
                } else if (confs.stateConfs.lowestScoringConf.getScore() < wtEnergy+eW && pfUB > wtPfLB - orderMag) { //otherwise, add the sequence if it's partition function upper bound is within orderMag of the wt pf lower bound
                    infos.add(info);
                    seqNum++;
                    reportSequence(infos.size() == 1, info);
                }
            }
        }

        List<SequenceInfo> badPfs = new ArrayList<>();
        for(int i=0; i<wtSpot; i++) {
            if (infos.get(i).pfUB < wtPfLB - orderMag || infos.get(i).lowestConf.getScore() > wtEnergy+eW) {
                badPfs.add(infos.get(i));
            }
        }

        infos.removeAll(badPfs);

        log("");
        if (infos.isEmpty()) {
            log("EWAK* didn't find any sequences within the window that satisfy all the constraints.");
        } else {
            log("EWAK* found the best %d sequences within %d orders of magnitude of the wild-type sequence partition function.", infos.size(), orderMag);
        }

        if(didEW){
            log("Found all sequences within %1.0f kcal/mol of the wild-type sequence.", eW);
        } else {
            log("The energy window of %1.0f kcal/mol was not fully enumerated.", eW);
        }

        if(didSeqMax){
            log("The sequence number max of %d was enumerated.", seqsToUse.size());
        } else if(didEW){
            log("The energy window was completed before enumerating the requested number of sequences.");
        } else if(!didSeqMax && !didEW){
            log("Exhaustively enumerated all sequence possibilities!");
        }

        return infos;

    }

    private void log(String msg, Object ... args) {
        if (printToConsole) {
            edu.duke.cs.osprey.tools.Log.log(msg, args);
        }
    }

    public BigInteger getNumConfsForSeq(Sequence seq, SimpleConfSpace confSpace){
        BigInteger myNum = seq.makeRCs(confSpace).getNumConformations();
        return myNum;
    }

    private void reportSequence(boolean isFirstSequence, SequenceInfo info) {

        if (printToConsole) {
            int cellSize = info.sequence.calcCellSize();
            log("\nFully calculated sequence: %s      %s\n                           %s",
                    info.sequence.toString(Sequence.Renderer.ResNum, cellSize),
                    info.sequence.isWildType() ? "Wild-type" : "",
                    info.sequence.toString(Sequence.Renderer.ResTypeMutations, cellSize)
            );
            log("\tState: %-3s    Conformation minimized lower bound: %12.6f\n",
                    state.name,
                    info.lowestConf.getScore()
            );
        }

        if (logFile != null) {

            // build the log record in TSV format

            // make the header if needed
            StringBuilder header = null;
            if (isFirstSequence) {
                header = new StringBuilder();
                for (SeqSpace.Position pos : seqSpace.positions) {
                    if (pos.index > 0) {
                        header.append("\t");
                    }
                    header.append(pos.resNum);
                }

                header.append("\t");
                header.append(state.name);
                header.append(" Lowest minBound Energy\t");
                header.append(state.name);
                header.append(" Conf");

            }

            // make the sequence entry
            StringBuilder buf = new StringBuilder();
            for (SeqSpace.Position pos : seqSpace.positions) {
                if (pos.index > 0) {
                    buf.append("\t");
                }
                buf.append(info.sequence.get(pos).mutationName());
            }
            buf.append("\t");
            ConfSearch.ScoredConf lowconf = info.lowestConf;
            buf.append("\t");
            buf.append(String.format("%.6f", lowconf.getScore()));
            buf.append("\t");
            buf.append(Conf.toString(lowconf.getAssignments()));

            // write to the log file, but don't keep the file open
            try (FileWriter out = new FileWriter(logFile, !isFirstSequence)) {

                if (header != null) {
                    out.write(header.toString());
                    out.write("\n");
                }
                out.write(buf.toString());
                out.write("\n");

            } catch (IOException ex) {

                // don't take down the whole job if we can't log the sequence for some reason
                // just report the error, and dump info to the console
                ex.printStackTrace(System.err);
                if (header != null) {
                    System.err.println(header);
                }
                System.err.println(buf);
            }
        }
    }
}
