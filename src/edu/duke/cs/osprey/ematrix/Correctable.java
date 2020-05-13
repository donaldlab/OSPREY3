package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.RCTupleContainer;

public interface Correctable<T extends RCTupleContainer> {

    void insertCorrection(T t);

    boolean containsCorrectionFor(RCTuple tup);
}
