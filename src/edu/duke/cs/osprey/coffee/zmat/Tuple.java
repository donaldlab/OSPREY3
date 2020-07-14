package edu.duke.cs.osprey.coffee.zmat;

import edu.duke.cs.osprey.confspace.TripleMatrix;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.DataOutput;
import java.io.IOException;


public interface Tuple {

	int type();
	ConfEnergyCalculator.MinimizationJob makeJob(ConfSpace confSpace, PosInterGen posInterGen);
	void setZ(BigExp z);
	boolean write(TupleMatrixGeneric<BigExp> singlesPairs, TripleMatrix<BigExp> triples);
	void write(DataOutput out) throws IOException;
}
