package edu.duke.cs.osprey.coffee.zmat;

import edu.duke.cs.osprey.confspace.TripleMatrix;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.confspace.compiled.ConfSpace;
import edu.duke.cs.osprey.energy.compiled.ConfEnergyCalculator;
import edu.duke.cs.osprey.energy.compiled.PosInterGen;
import edu.duke.cs.osprey.tools.BigExp;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;


public class Triple implements Tuple {

	public static final int Type = 3;

	public int posi1;
	public int confi1;
	public int posi2;
	public int confi2;
	public int posi3;
	public int confi3;
	public BigExp z;

	public Triple(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3, BigExp z) {
		this.posi1 = posi1;
		this.confi1 = confi1;
		this.posi2 = posi2;
		this.confi2 = confi2;
		this.posi3 = posi3;
		this.confi3 = confi3;
		this.z = z;
	}

	public Triple(int posi1, int confi1, int posi2, int confi2, int posi3, int confi3) {
		this(posi1, confi1, posi2, confi2, posi3, confi3, null);
	}

	@Override
	public int type() {
		return Type;
	}

	@Override
	public ConfEnergyCalculator.MinimizationJob makeJob(ConfSpace confSpace, PosInterGen posInterGen) {
		return new ConfEnergyCalculator.MinimizationJob(
			confSpace.assign(posi1, confi1, posi2, confi2, posi3, confi3),
			posInterGen.tripleCorrection(confSpace, posi1, confi1, posi2, confi2, posi3, confi3)
		);
	}

	@Override
	public void setZ(BigExp z) {
		this.z = z;
	}

	@Override
	public boolean write(TupleMatrixGeneric<BigExp> singlesPairs, TripleMatrix<BigExp> triples) {

		// convert the triple energy into a correction
		BigExp divisor = new BigExp(singlesPairs.getPairwise(posi1, confi1, posi2, confi2));
		divisor.mult(singlesPairs.getPairwise(posi1, confi1, posi3, confi3));
		divisor.mult(singlesPairs.getPairwise(posi2, confi2, posi3, confi3));
		divisor.pow(1.0/MathTools.numTriplesPerPair(singlesPairs.getNumPos()));

		// filter out unhelpful corrections
		// ie, if the correction factor is >= 1
		if (z.greaterThanOrEqual(divisor)) {
			return false;
		}

		var correction = new BigExp(z);
		correction.div(divisor);
		triples.set(posi1, confi1, posi2, confi2, posi3, confi3, correction);
		return true;
	}

	@Override
	public void write(DataOutput out)
	throws IOException {
		out.writeInt(posi1);
		out.writeInt(confi1);
		out.writeInt(posi2);
		out.writeInt(confi2);
		out.writeInt(posi3);
		out.writeInt(confi3);
		out.writeDouble(z.fp);
		out.writeInt(z.exp);
	}

	public static Triple read(DataInput in)
	throws IOException {
		int posi1 = in.readInt();
		int confi1 = in.readInt();
		int posi2 = in.readInt();
		int confi2 = in.readInt();
		int posi3 = in.readInt();
		int confi3 = in.readInt();
		double fp = in.readDouble();
		int exp = in.readInt();
		var z = new BigExp(fp, exp);
		return new Triple(posi1, confi1, posi2, confi2, posi3, confi3, z);
	}
}
