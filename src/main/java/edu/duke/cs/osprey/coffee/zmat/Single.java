package edu.duke.cs.osprey.coffee.zmat;

import edu.duke.cs.osprey.confspace.TripleMatrix;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.tools.BigExp;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;


public class Single implements Tuple {

	public static final int Type = 1;

	public int posi1;
	public int confi1;
	public BigExp z;

	public Single(int posi1, int confi1, BigExp z) {
		this.posi1 = posi1;
		this.confi1 = confi1;
		this.z = z;
	}

	public Single(int posi1, int confi1) {
		this(posi1, confi1, null);
	}

	@Override
	public int type() {
		return Type;
	}

	@Override
	public ConfEnergyCalculator.MinimizationJob makeJob(ConfSpace confSpace, PosInterGen posInterGen) {
		return new ConfEnergyCalculator.MinimizationJob(
			confSpace.assign(posi1, confi1),
			posInterGen.single(confSpace, posi1, confi1)
		);
	}

	@Override
	public void setZ(BigExp z) {
		this.z = z;
	}

	@Override
	public boolean write(TupleMatrixGeneric<BigExp> singlesPairs, TripleMatrix<BigExp> triples) {
		singlesPairs.setOneBody(posi1, confi1, z);
		return true;
	}

	@Override
	public void write(DataOutput out)
	throws IOException {
		out.writeInt(posi1);
		out.writeInt(confi1);
		out.writeDouble(z.fp);
		out.writeInt(z.exp);
	}

	public static Single read(DataInput in)
	throws IOException {
		int posi1 = in.readInt();
		int confi1 = in.readInt();
		double fp = in.readDouble();
		int exp = in.readInt();
		var z = new BigExp(fp, exp);
		return new Single(posi1, confi1, z);
	}
}
