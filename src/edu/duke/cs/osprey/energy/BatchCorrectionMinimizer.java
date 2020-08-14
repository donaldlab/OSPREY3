package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupE;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;
import edu.duke.cs.osprey.ematrix.UpdatingEnergyMatrix;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import org.apache.commons.math3.analysis.function.Cos;

import java.util.*;

import static edu.duke.cs.osprey.tools.Log.log;

public class BatchCorrectionMinimizer {

    public final ConfEnergyCalculator confEcalc;
    public final int CostThreshold = 100;
    public final int[] costs = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    private Queue<Batch> batches;
    private UpdatingEnergyMatrix.TupleTrie submittedConfs;
    private UpdatingEnergyMatrix correctionMatrix;
    public EnergyMatrix minimizingEnergyMatrix;
    private final Queue<PartialMinimizationTuple> waitingQueue;
    private int waitingCost = 0;
    private final Object lock = new Object();

    public BatchCorrectionMinimizer(ConfEnergyCalculator confEcalc, UpdatingEnergyMatrix correctionMatrix,
                                    EnergyMatrix minimizingEnergyMatrix) {
        this.confEcalc = confEcalc;
        this.correctionMatrix = correctionMatrix;
        this.minimizingEnergyMatrix = minimizingEnergyMatrix;
        this.submittedConfs = new UpdatingEnergyMatrix.TupleTrie(confEcalc.confSpace.positions);
        //TODO: redo the tupleTrie like was in the other branch
        this.waitingQueue = new LinkedList<>();
        this.batches = new LinkedList<>();
    }

    public synchronized void addTuple(PartialMinimizationTuple tuple) {
        if (submittedConfs.contains(tuple.tup))
            return;
        submittedConfs.insert(new TupE(tuple.tup, 0));
        int tupleSize = tuple.size();
        if(costs[tupleSize] < 0)
            costs[tupleSize] = confEcalc.makeFragInters(tuple.tup).size();
        waitingQueue.add(tuple);
        waitingCost += costs[tupleSize];
    }

    public synchronized boolean canBatch(){
        return waitingCost >= CostThreshold;
    }

    public void makeBatch(){
        Batch newBatch = new Batch();
        int batchCost = 0;
        while(batchCost < CostThreshold){
            synchronized(lock){
                PartialMinimizationTuple frag = waitingQueue.poll();
                batchCost += costs[frag.size()];
                waitingCost -= costs[frag.size()];
                newBatch.fragments.add(frag);
            }
        }
        synchronized(lock) {
            this.batches.add(newBatch);
        }
    }

    public Batch getBatch(){
        synchronized(lock){
            return this.batches.poll();
        }
    }

    public static class Batch {
        public List<PartialMinimizationTuple> fragments = new ArrayList<>();
        int cost = 0;
    }

    public boolean isParametricallyIncompatible(RCTuple tuple) {
        for (int i1=0; i1<tuple.size(); i1++) {
            SimpleConfSpace.ResidueConf rc1 = getRC(tuple, i1);
            for (int i2=0; i2<i1; i2++) {
                SimpleConfSpace.ResidueConf rc2 = getRC(tuple, i2);
                if (!isPairParametricallyCompatible(rc1, rc2)) {
                    return true;
                }
            }
        }
        return false;
    }

    private SimpleConfSpace.ResidueConf getRC(RCTuple tuple, int index) {
        return confEcalc.confSpace.positions.get(tuple.pos.get(index)).resConfs.get(tuple.RCs.get(index));
    }

    private boolean isPairParametricallyCompatible(SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.ResidueConf rc2) {
        for(String dofName : rc1.dofBounds.keySet()){
            if(rc2.dofBounds.containsKey(dofName)){
                //shared DOF between the RCs; make sure the interval matches
                double[] interval1 = rc1.dofBounds.get(dofName);
                double[] interval2 = rc2.dofBounds.get(dofName);
                for(int a=0; a<2; a++){
                    if( Math.abs(interval1[a] - interval2[a]) > 1e-8 ){
                        return false;
                    }
                }
            }
        }
        return true;//found no incompatibilities
    }

    /** PartialMinimizationTuple
     *
     * This class contains an RCTuple representing a conformation fragment, along with the minimized energy and
     * energy lowerbound of the full conformation that originally contained it.
     */
    public static class PartialMinimizationTuple extends TupE {
        public double parentConfLB;
        public RCTuple parentConf;

        public PartialMinimizationTuple(RCTuple tup, double E) {
            super(tup, E);
        }
        public PartialMinimizationTuple(RCTuple tup, double E, double parentConfLB, RCTuple parentConf) {
            super(tup, E);
            this.parentConfLB = parentConfLB;
            this.parentConf = parentConf;
        }

        public PartialMinimizationTuple(String repr) {
            super(repr);
        }
    }
}
