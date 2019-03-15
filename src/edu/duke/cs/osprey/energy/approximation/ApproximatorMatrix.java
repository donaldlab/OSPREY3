package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction.Approximator;
import edu.duke.cs.osprey.tools.IOable;

import java.io.*;
import java.util.*;


public class ApproximatorMatrix implements IOable {

	public final SimpleConfSpace confSpace;

	private final int[] offsets;
	private final Map<String,Approximator.Addable[]> approximators;

	public ApproximatorMatrix(SimpleConfSpace confSpace) {

		this.confSpace = confSpace;

		// compute the offsets
		int offset = 0;
		offsets = new int[confSpace.positions.size()];
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			offsets[pos.index] = offset;
			offset += pos.resConfs.size();
		}

		approximators = new LinkedHashMap<>();
		for (String fixedResNum : confSpace.shellResNumbers) {
			approximators.put(fixedResNum, new Approximator.Addable[offset]);
		}
	}

	private int getIndex(int pos, int rc) {
		return offsets[pos] + rc;
	}

	public Approximator.Addable get(String fixedResNum, int pos, int rc) {
		Approximator.Addable[] approximators = this.approximators.get(fixedResNum);
		if (approximators == null) {
			return null;
		}
		return approximators[getIndex(pos, rc)];
	}

	public void set(String fixedResNum, int pos, int rc, Approximator.Addable approximator) {
		approximators.get(fixedResNum)[getIndex(pos, rc)] = approximator;
	}

	private static class ResPairApproximator {

		public final ResidueInteractions.Pair pair;
		public final Approximator.Addable approximator;

		public ResPairApproximator(ResidueInteractions.Pair pair, Approximator.Addable approximator) {
			this.pair = pair;
			this.approximator = approximator;
		}
	}

	public ResidueInteractionsApproximator get(RCTuple tuple, ResidueInteractions inters, double errorBudget) {

		double errorBudgetPerResidue = errorBudget/tuple.size();

		// TODO: optimize this?

		// TODO: does this match the DOF order that the conf space picks?
		Approximator[] approximators = new Approximator[tuple.size()];
		ResidueInteractions approxInters = new ResidueInteractions();

		for (int i=0; i<tuple.size(); i++) {
			int pos = tuple.pos.get(i);
			int rc = tuple.RCs.get(i);

			String resNum = confSpace.positions.get(pos).resNum;

			// collect the approximators for the fixed residues
			List<ResPairApproximator> resPairApproximators = new ArrayList<>();
			for (ResidueInteractions.Pair pair : inters) {
				 String fixedResNum = pair.getOtherResNum(resNum);
				 if (fixedResNum != null) {
				 	Approximator.Addable approximator = get(fixedResNum, pos, rc);
					if (approximator != null) {
				 		resPairApproximators.add(new ResPairApproximator(pair, approximator));
					}
				 }
			}

			// no approximators good enough? make a no-op approximator
			if (resPairApproximators.isEmpty()) {
				int d = confSpace.makeMolecule(new RCTuple(pos, rc)).dofs.size();
				approximators[i] = new NOPApproximator(d);
				continue;
			}

			// add as many approximators as we can and stay under the error budget
			resPairApproximators.sort(Comparator.comparing(a -> a.approximator.error()));
			Approximator.Addable approximator = resPairApproximators.get(0).approximator.makeIdentity();
			for (ResPairApproximator pair : resPairApproximators) {
				if (approximator.error() + pair.approximator.error() > errorBudgetPerResidue) {
					break;
				}
				approxInters.add(pair.pair);
				approximator.add(pair.approximator, pair.pair.weight, pair.pair.offset);
			}

			approximators[i] = approximator;
		}

		// combine the approximators for all the residues
		Approximator approximator = new ApproximatedObjectiveFunction.Approximators(approximators);

		/* DEBUG
		log("approximated %d/%d (%.1f%%) inters   error = %12.6f/%12.6f    for tuple %s",
			approxInters.size(), inters.size(),
			100.0*approxInters.size()/inters.size(),
			approximator.error(), errorBudget,
			tuple
		);
		*/

		return new ResidueInteractionsApproximator(approxInters, approximator);
	}

	private byte getType(Approximator approximator) {
		if (approximator == null) {
			return 0;
		} else if (approximator instanceof QuadraticApproximator) {
			return 1;
		} else {
			throw new IllegalArgumentException("unrecognized approximator type: " + approximator.getClass().getName());
		}
	}

	private Approximator alloc(byte type, int d) {
		switch (type) {
			case 0: return null;
			case 1: return new QuadraticApproximator(d);
			default: throw new IllegalArgumentException("unrecognized approximator type: " + type);
		}
	}

	@Override
	public void writeTo(DataOutput out)
	throws IOException {

		for (String fixedResNum : confSpace.shellResNumbers) {
			for (Approximator approximator : approximators.get(fixedResNum)) {
				out.writeByte(getType(approximator));
				if (approximator != null) {
					out.writeInt(approximator.numDofs());
					IOable.of(approximator).writeTo(out);
				}
			}
		}
	}

	@Override
	public void readFrom(DataInput in)
	throws IOException {

		for (String fixedResNum : confSpace.shellResNumbers) {
			Approximator[] approximators = this.approximators.get(fixedResNum);
			for (int i=0; i<approximators.length; i++) {
				byte type = in.readByte();
				if (type != 0) {
					approximators[i] = alloc(type, in.readInt());
					IOable.of(approximators[i]).readFrom(in);
				} else {
					approximators[i] = null;
				}
			}
		}
	}
}
