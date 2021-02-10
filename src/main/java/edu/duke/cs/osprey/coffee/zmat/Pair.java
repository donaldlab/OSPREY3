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


public class Pair implements Tuple {

	public static final int Type = 2;

	public int posi1;
	public int confi1;
	public int posi2;
	public int confi2;
	public BigExp z;

	public Pair(int posi1, int confi1, int posi2, int confi2, BigExp z) {
		this.posi1 = posi1;
		this.confi1 = confi1;
		this.posi2 = posi2;
		this.confi2 = confi2;
		this.z = z;
	}

	public Pair(int posi1, int confi1, int posi2, int confi2) {
		this(posi1, confi1, posi2, confi2, null);
	}

	@Override
	public int type() {
		return Type;
	}

	@Override
	public ConfEnergyCalculator.MinimizationJob makeJob(ConfSpace confSpace, PosInterGen posInterGen) {
		return new ConfEnergyCalculator.MinimizationJob(
			confSpace.assign(posi1, confi1, posi2, confi2),
			posInterGen.pair(confSpace, posi1, confi1, posi2, confi2)
		);
	}

	@Override
	public void setZ(BigExp z) {
		this.z = z;
	}

	@Override
	public boolean write(TupleMatrixGeneric<BigExp> singlesPairs, TripleMatrix<BigExp> triples) {
		singlesPairs.setPairwise(posi1, confi1, posi2, confi2, z);
		return true;
	}

	@Override
	public void write(DataOutput out)
	throws IOException {
		out.writeInt(posi1);
		out.writeInt(confi1);
		out.writeInt(posi2);
		out.writeInt(confi2);
		out.writeDouble(z.fp);
		out.writeInt(z.exp);
	}

	public static Pair read(DataInput in)
	throws IOException {
		int posi1 = in.readInt();
		int confi1 = in.readInt();
		int posi2 = in.readInt();
		int confi2 = in.readInt();
		double fp = in.readDouble();
		int exp = in.readInt();
		var z = new BigExp(fp, exp);
		return new Pair(posi1, confi1, posi2, confi2, z);
	}
}
