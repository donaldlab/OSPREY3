package edu.duke.cs.osprey.newEwakstar;

import edu.duke.cs.osprey.astar.conf.ConfAStarTree;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.astar.seq.RTs;
import edu.duke.cs.osprey.astar.seq.nodes.SeqAStarNode;
import edu.duke.cs.osprey.astar.seq.SeqAStarTree;
import edu.duke.cs.osprey.astar.seq.order.SequentialSeqAStarOrder;
import edu.duke.cs.osprey.astar.seq.scoring.NOPSeqAStarScorer;
import edu.duke.cs.osprey.astar.seq.scoring.SeqAStarScorer;
import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.energy.ConfEnergyCalculator;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.tools.HashCalculator;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.util.*;
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
public class Ewakstar {

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
        public Function<RCs,ConfAStarTree> confTreeFactory;

        /**
         * set this to a File if you want to use ConfDB with COMETS
         */
        public File confDBFile = null;

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
            if (confTreeFactory == null) {
                throw new InitException(this, "confTreeFactory");
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

    private class ConfDBs extends ConfDB.DBs {

        public ConfDB.ConfTable table;

        public ConfDBs() {

            // open the DBs
            add(state.confSpace, state.confDBFile);

            // make the table
            ConfDB confdb = get(state.confSpace);
            if (confdb != null) {
                table = confdb.new ConfTable("COMETS");
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
            confTree = state.confTreeFactory.apply(rcs);
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
                wtConf = state.confEcalc.calcEnergy(confTree.nextConf());
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

                Sequence sequence = seqNode.makeSequence(seqSpace)
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
            stateConfs.getLowestScoringConf();

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
        private double boundEw = 10.0;

        /** The maximum number of simultaneous residue mutations to consider for each sequence mutant */
        private int maxSimultaneousMutations = 1;

        private boolean printToConsole = true;

        /** File to which to log sequences as they are found */
        private File logFile = null;

        public Builder addState(State s){
            state = s;
            return this;
        }

        public Builder setMaxSimultaneousMutations(int val) {
            maxSimultaneousMutations = val;
            return this;
        }

        public Builder setBoundEw(double val){
            boundEw = val;
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

        public Ewakstar build() {
            return new Ewakstar(state, boundEw, maxSimultaneousMutations, printToConsole, logFile);
        }
    }


    public final double boundEw;
    public final int maxSimultaneousMutations;
    public final boolean printToConsole;
    public final File logFile;

    public final State state;
    public final SeqSpace seqSpace;

    private final Map<StateConfs.Key,StateConfs> stateConfsCache = new HashMap<>();

    private Ewakstar(State state, double boundEw, int maxSimultaneousMutations, boolean printToConsole, File logFile) {

        this.boundEw = boundEw;
        this.maxSimultaneousMutations = maxSimultaneousMutations;
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
    public List<SequenceInfo> findBestSequences(int numSequences, int orderMag, State PL) {

        // reset any previous state
        stateConfsCache.clear();

        state.checkConfig();

        List<SequenceInfo> bestPLSeqs = extractPLSeqsByLB(numSequences, orderMag, PL);

        return bestPLSeqs;

    }

    public List<SequenceInfo> extractPLSeqsByLB(int numSequences, int orderMag, State PL) {

        // start the A* search over sequences
        SeqAStarTree seqTree = new SeqAStarTree.Builder(new RTs(seqSpace))
                .setHeuristics(
                        new SequentialSeqAStarOrder(),
                        new NOPSeqAStarScorer(),
                        new SeqHScorer()
                )
                .setMaxSimultaneousMutations(maxSimultaneousMutations)
                .build();

        List<SequenceInfo> infos = new ArrayList<>();
        double pfUB;
        double wtPfLB = Double.POSITIVE_INFINITY;
        boolean didSeqMax = false;

        Boolean wtFound = false;
        double wtEnergy = 0.0;
        Boolean didEW = false;
        int wtSpot = 0;

        BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

        try (ConfDBs confDBs = new ConfDBs()) {

            for(int seqNum=0; seqNum<numSequences; seqNum++) {

                // get the next sequence from the tree
                SeqAStarNode node = seqTree.nextLeafNode();
                if (node == null) {
                    break;
                }

                SeqConfs confs = (SeqConfs) node.getData();
                if (confs == null) {

                    log("Estimated sequence: %s   objective lower bound: %12.6f",
                            node.makeSequence(seqSpace),
                            node.getScore()
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
                    System.out.println("Finding all sequences within " + boundEw + " kcal of WT minimized energy: "
                            + wtEnergy + " or until " + numSequences + " sequences are enumerated.");
                }

                Sequence newSequence = node.makeSequence(PL.confSpace.seqSpace);
                pfUB = Math.log10(bc.calc(node.getScore()).multiply(new BigDecimal(getNumConfsForSeq(newSequence, PL.confSpace))).doubleValue());
                SequenceInfo info = new SequenceInfo(node, confs, pfUB);
                if (!wtFound) {
                    infos.add(info);
                } else {
                    if(info.lowestConf.getScore() > wtEnergy+boundEw) { //if we have surpassed the energy window we want to quit
                        if(infos.size() >= numSequences)
                            didSeqMax = true;
                        didEW = true;
                        break;
                    } else if (infos.size() >= numSequences) { //if we have surpassed the maximum number of sequence we want to find, we want to quit
                        didSeqMax = true;
                        break;
                    } else if(pfUB > wtPfLB - orderMag) { //otherwise, add the sequence if it's partition function upper bound is within orderMag of the wt pf lower bound
                        infos.add(info);
                        reportSequence(infos.size() == 1, info);
                    }
                }
            }

            for(int i=0; i<wtSpot; i++) {
                if (infos.get(i).lowestConf.getScore() < wtPfLB - orderMag) {
                    infos.get(i).setCanPrune(true);
                }
            }

            for (SequenceInfo si : infos){
                if (si.canPrune)
                    infos.remove(si);
            }

            log("");
            if (infos.isEmpty()) {
                log("EWAK* didn't find any sequences within the window that satisfy all the constraints.");
            } else {
                log("EWAK* found the best %d sequences within %d orders of magnitude of the wild-type sequence partition function.", infos.size(), orderMag);
            }

            if(didEW){
                log("Found all sequences within %1.0f kcal/mol of the wild-type sequence.", boundEw);
            } else {
                log("The energy window of %1.0f kcal/mol was not fully enumerated.", boundEw);
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
    }

    private void log(String msg, Object ... args) {
        if (printToConsole) {
            edu.duke.cs.osprey.tools.Log.log(msg, args);
        }
    }

    public BigInteger getNumConfsForSeq(Sequence seq, SimpleConfSpace confSpace){
        return seq.makeRCs(confSpace).getNumConformations();
    }

    private void reportSequence(boolean isFirstSequence, SequenceInfo info) {

        if (printToConsole) {
            int cellSize = info.sequence.calcCellSize();
            log("\nFully calculated sequence: %s      %s\n                           %s",
                    info.sequence.toString(Sequence.Renderer.ResNum, cellSize),
                    info.sequence.isWildType() ? "Wild-type" : "",
                    info.sequence.toString(Sequence.Renderer.ResTypeMutations, cellSize)
            );
            log("\tState: %-3s    Lower Bound on Minimized Energy: %12.6f\n",
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
