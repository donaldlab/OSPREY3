package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.ConfSpace;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;

public class DofMatrix extends TupleMatrixGeneric<double[]> {

	private static final long serialVersionUID = 7381812984847056950L;
	
	public DofMatrix(ConfSpace cSpace) {
		super(cSpace, 0, null);
	}
}
