package edu.duke.cs.osprey.design.analysis;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;
import edu.duke.cs.osprey.confspace.ConfSearch;

import java.util.List;
import java.util.PriorityQueue;
import java.util.stream.Collectors;

public class EnergiedConfQueue {

    private final int numElements;

    private Ordering<ConfSearch.EnergiedConf> lowestEnergyFirstOrdering = new Ordering<>() {
        public int compare(ConfSearch.EnergiedConf left, ConfSearch.EnergiedConf right) {
            return Doubles.compare(left.getEnergy(), right.getEnergy());
        }
    };

    private final PriorityQueue<ConfSearch.EnergiedConf> confQueue = new PriorityQueue<>(lowestEnergyFirstOrdering.reverse());

    public EnergiedConfQueue(int maxElements) {
        this.numElements = maxElements;
    }

    public List<ConfSearch.EnergiedConf> toOrderedList() {
        return confQueue.stream().sorted(lowestEnergyFirstOrdering).collect(Collectors.toList());
    }

    public void add(ConfSearch.EnergiedConf conf) {
        if (confQueue.size() < numElements) {
            confQueue.add(conf);
            return;
        }

        assert confQueue.size() > 0;
        var maxConf = confQueue.peek();
        if (conf.getEnergy() < maxConf.getEnergy()) {
            confQueue.poll(); // throw out the conf with the maximum energy
            confQueue.add(conf);
        }
    }
}
