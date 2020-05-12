package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.TupE;

public interface Correctable<T extends TupE> {

    void insertCorrection(T t);

    boolean containsCorrectionFor(RCTuple tup);
}
