package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.confspace.*;
import edu.duke.cs.osprey.ematrix.Correctable;
import edu.duke.cs.osprey.ematrix.EnergyMatrix;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BatchCorrectionMinimizer {

    private final ConfEnergyCalculator confEcalc;
    public final int CostThreshold = 1000;
    public final int[] costs = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    private Batch batch = null;
    private final TupleTrie<TupE> submittedConfs;
    private final Correctable<TupE> correctionStorage;
    private final EnergyMatrix minimizingEnergyMatrix;
    public BatchCorrectionMinimizer(ConfEnergyCalculator confEcalc, Correctable<TupE> correctionStorage,
                                    EnergyMatrix minimizingEnergyMatrix) {
        this.confEcalc = confEcalc;
        this.correctionStorage = correctionStorage;
        this.minimizingEnergyMatrix = minimizingEnergyMatrix;
        submittedConfs = new TupleTrie<>(confEcalc.confSpace.positions);
    }

    public Batch getBatch() {
            if (batch == null) {
                batch = new Batch();
            }
            return batch;
        }

    public void submitIfFull() {
            if (batch != null && batch.cost >= CostThreshold) {
                submit();
            }
        }

    public void submit() {
            if (batch != null) {
                batch.submitTask();
                batch = null;
            }
        submittedConfs.clear();
    }

    public class Batch {

        List<RCTuple> fragments = new ArrayList<>();
        int cost = 0;

        public void addTuple(RCTuple tuple) {
            synchronized (this) {
                if (submittedConfs.contains(tuple))
                    return;
                submittedConfs.insert(new TupE(tuple, 0));
            }
            int tupleSize = tuple.size();
            if(costs[tupleSize] < 0)
                costs[tupleSize] = confEcalc.makeFragInters(tuple).size();
            fragments.add(tuple);
            cost += costs[tupleSize];
        }

        void submitTask() {
            confEcalc.tasks.submit(
                    () -> {

                        // calculate all the fragment energies
                        Map<RCTuple, Double> confs = new HashMap<>();
                        for (RCTuple frag : fragments) {

                            double energy;

                            // are there any RCs are from two different backbone states that can't connect?
                            if (isParametricallyIncompatible(frag)) {

                                // yup, give this frag an infinite energy so we never choose it
                                energy = Double.POSITIVE_INFINITY;

                            } else {

                                // nope, calculate the usual fragment energy
                                energy = confEcalc.calcEnergy(frag).energy;
                            }
                            // Record the filtered (possibly infinite) energy
                            confs.put(frag, energy);
                        }
                        System.out.println("Minimized "+fragments.size()+" tuples.");
                        return confs;
                    },
                    (Map<RCTuple, Double> confs) -> {
                        // update the energy matrix
                        for(RCTuple tuple : confs.keySet()) {
                            double lowerbound = minimizingEnergyMatrix.getInternalEnergy(tuple);
                            double tupleEnergy = confs.get(tuple);
                            if (tupleEnergy - lowerbound > 0) {
                                double correction = tupleEnergy - lowerbound;
                                RCTuple orderedTup = tuple.sorted();
                                correctionStorage.insertCorrection(new TupE(orderedTup, correction));
                            } else
                                System.err.println("Negative correction for " + tuple.stringListing());


                        }
                    }
            );
        }
    }



    protected boolean isParametricallyIncompatible(RCTuple tuple) {
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
}
