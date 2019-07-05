/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

package edu.duke.cs.osprey.energy.approximation;

import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.TupleMatrixGeneric;
import edu.duke.cs.osprey.dof.DofInfo;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction.Approximator;
import edu.duke.cs.osprey.tools.IOable;

import java.io.*;
import java.util.*;


/**
 * Forcefield approximators for RC,residue pairs in a conf space
 */
public class ApproximatorMatrix implements IOable {

	public final SimpleConfSpace confSpace;

	private final int[] offsets;
	private final List<Map<String,Approximator.Addable>> fixedApproximators;
	private final TupleMatrixGeneric<Approximator.Addable> tupleApproximators;

	public ApproximatorMatrix(SimpleConfSpace confSpace) {

		this.confSpace = confSpace;

		// compute the offsets
		int offset = 0;
		offsets = new int[confSpace.positions.size()];
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			offsets[pos.index] = offset;
			offset += pos.resConfs.size();
		}

		fixedApproximators = new ArrayList<>(offset);
		for (int i=0; i<offset; i++) {
			fixedApproximators.add(new HashMap<>());
		}

		tupleApproximators = new TupleMatrixGeneric<>(confSpace);
	}

	public Approximator.Addable get(int pos1, int rc1) {
		return tupleApproximators.getOneBody(pos1, rc1);
	}
	public void set(int pos1, int rc1, Approximator.Addable approximator) {
		tupleApproximators.setOneBody(pos1, rc1, approximator);
	}
	public Approximator.Addable get(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1) {
		return tupleApproximators.getOneBody(pos1, rc1);
	}
	public void set(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, Approximator.Addable approximator) {
		tupleApproximators.setOneBody(pos1, rc1, approximator);
	}

	public Approximator.Addable get(int pos1, int rc1, int pos2, int rc2) {
		return tupleApproximators.getPairwise(pos1, rc1, pos2, rc2);
	}
	public void set(int pos1, int rc1, int pos2, int rc2, Approximator.Addable approximator) {
		tupleApproximators.setPairwise(pos1, rc1, pos2, rc2, approximator);
	}
	public Approximator.Addable get(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf rc2) {
		return tupleApproximators.getPairwise(pos1, rc1, pos2, rc2);
	}
	public void set(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, SimpleConfSpace.Position pos2, SimpleConfSpace.ResidueConf rc2, Approximator.Addable approximator) {
		tupleApproximators.setPairwise(pos1, rc1, pos2, rc2, approximator);
	}

	private int getIndex(int pos, int rc) {
		return offsets[pos] + rc;
	}

	public Approximator.Addable get(int pos1, int rc1, String fixedResNum) {
		return fixedApproximators.get(getIndex(pos1, rc1)).get(fixedResNum);
	}
	public void set(int pos1, int rc1, String fixedResNum, Approximator.Addable approximator) {
		fixedApproximators.get(getIndex(pos1, rc1)).put(fixedResNum, approximator);
	}
	public Approximator.Addable get(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, String fixedResNum) {
		return get(pos1.index, rc1.index, fixedResNum);
	}
	public void set(SimpleConfSpace.Position pos1, SimpleConfSpace.ResidueConf rc1, String fixedResNum, Approximator.Addable approximator) {
		set(pos1.index, rc1.index, fixedResNum, approximator);
	}

	private static class InteractionApproximator {

		final ResidueInteractions.Pair inter;
		final Approximator.Addable approximator;

		InteractionApproximator(ResidueInteractions.Pair inter, Approximator.Addable approximator) {
			this.inter = inter;
			this.approximator = approximator;
		}
	}

	public ResidueInteractionsApproximator get(RCTuple tuple, ResidueInteractions inters, double errorBudget) {

		// TODO: optimize me!

		DofInfo dofInfo = confSpace.makeDofInfo(tuple);

		// DEBUG
		//log("%s", dofInfo);

		double errorBudgetPerInter = errorBudget/inters.size();

		ResidueInteractionsApproximator.Builder builder = new ResidueInteractionsApproximator.Builder(dofInfo);

		// which residue interactions can be approximated?
		List<InteractionApproximator> leftovers = new ArrayList<>();
		for (ResidueInteractions.Pair inter : inters) {

			// find out if any design positions are involved in this residue interaction
			Integer blockIndex1 = dofInfo.getBlockIndex(inter.resNum1);
			Integer blockIndex2 = dofInfo.getBlockIndex(inter.resNum2);

			// get the approximator for this residue interaction
			Approximator.Addable approximator = null;
			if (blockIndex1 != null && blockIndex2 != null) {

				if (blockIndex1.equals(blockIndex2)) {

					// single
					SimpleConfSpace.Position pos = dofInfo.positions.get(blockIndex1);
					SimpleConfSpace.ResidueConf rc = dofInfo.resConfs.get(blockIndex1);
					approximator = get(pos, rc);

				} else {

					// pair
					SimpleConfSpace.Position pos1 = dofInfo.positions.get(blockIndex1);
					SimpleConfSpace.ResidueConf rc1 = dofInfo.resConfs.get(blockIndex1);
					SimpleConfSpace.Position pos2 = dofInfo.positions.get(blockIndex2);
					SimpleConfSpace.ResidueConf rc2 = dofInfo.resConfs.get(blockIndex2);
					approximator = get(pos1, rc1, pos2, rc2);
				}

			} else if (blockIndex1 != null || blockIndex2 != null) {

				// fixed residue interaction
				int blockIndex = blockIndex1 != null ? blockIndex1 : blockIndex2;
				String fixedResNum = blockIndex1 != null ? inter.resNum2 : inter.resNum1;

				SimpleConfSpace.Position pos = dofInfo.positions.get(blockIndex);
				SimpleConfSpace.ResidueConf rc = dofInfo.resConfs.get(blockIndex);
				approximator = get(pos, rc, fixedResNum);

			} else {

				// just in case...
				assert (false) : String.format("residue interaction %s:%s appears unrelated to tuple %s",
					inter.resNum1, inter.resNum2, tuple
				);
			}

			if (approximator == null) {
				builder.dontApproximate(inter);
			} else if (approximator.error() <= errorBudgetPerInter) {
				builder.approximate(inter, approximator);
			} else {
				leftovers.add(new InteractionApproximator(inter, approximator));
			}
		}

		assert (builder.error() <= errorBudget);

		if (!leftovers.isEmpty()) {

			// sort the leftover approximators by error, we'll add whatever else we can and stay under budget
			leftovers.sort(Comparator.comparing(interApproximator -> interApproximator.approximator.error()));

			for (InteractionApproximator interApproximator : leftovers) {

				/* DEBUG
				log("leftover:  %3s:%3s  err=%8.6f   accept? %b   coefficients=%s",
					interApproximator.inter.resNum1, interApproximator.inter.resNum2,
					interApproximator.approximator.error(),
					builder.error() + interApproximator.approximator.error() <= errorBudget,
					Arrays.toString(((QuadraticApproximator)interApproximator.approximator).coefficients.toArray())
				);
				*/

				if (builder.error() + interApproximator.approximator.error() <= errorBudget) {
					builder.approximate(interApproximator.inter, interApproximator.approximator);
				} else {
					// we won't approximate this interaction, let the forcefield compute it as usual
					builder.dontApproximate(interApproximator.inter);
				}
			}
		}

		return builder.build();

		/* DEBUG
		log("approximated %2d + %2d of %2d (%5.1f%%) inters   error = %12.6f/%12.6f    for tuple %s",
			out.approxInters.size(), out.ffInters.size(), inters.size(),
			100.0*out.approxInters.size()/inters.size(),
			out.approximator.error(), errorBudget,
			tuple
		);
		*/
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

	private Approximator.Addable alloc(byte type, RCTuple tuple) {
		DofInfo dofInfo = confSpace.makeDofInfo(tuple);
		switch (type) {
			case 0: return null;
			case 1: return new QuadraticApproximator(dofInfo.ids, dofInfo.counts);
			default: throw new IllegalArgumentException("unrecognized approximator type: " + type);
		}
	}

	private void write(Approximator.Addable approximator, DataOutput out)
	throws IOException {
		out.writeByte(getType(approximator));
		if (approximator != null) {
			IOable.of(approximator).writeTo(out);
		}
	}

	@Override
	public void writeTo(DataOutput out)
	throws IOException {

		// singles
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
				write(get(pos1, rc1), out);
			}
		}

		// pairs
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
				for (SimpleConfSpace.Position pos2 : confSpace.positions.subList(0, pos1.index)) {
					for (SimpleConfSpace.ResidueConf rc2: pos2.resConfs) {
						write(get(pos1, rc1, pos2, rc2), out);
					}
				}
			}
		}

		// inters with fixed residues
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
				for (String resNum : confSpace.shellResNumbers) {
					write(get(pos1, rc1, resNum), out);
				}
			}
		}
	}

	private Approximator.Addable read(DataInput in, RCTuple tuple)
	throws IOException {
		byte type = in.readByte();
		if (type == 0) {
			return null;
		}
		Approximator.Addable approximator = alloc(type, tuple);
		IOable.of(approximator).readFrom(in);
		return approximator;
	}

	@Override
	public void readFrom(DataInput in)
	throws IOException {

		// singles
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
				set(pos1, rc1, read(in, new RCTuple(pos1.index, rc1.index)));
			}
		}

		// pairs
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
				for (SimpleConfSpace.Position pos2 : confSpace.positions.subList(0, pos1.index)) {
					for (SimpleConfSpace.ResidueConf rc2: pos2.resConfs) {
						set(pos1, rc1, pos2, rc2, read(in, new RCTuple(pos1.index, rc1.index, pos2.index, rc2.index)));
					}
				}
			}
		}

		// inters with fixed residues
		for (SimpleConfSpace.Position pos1 : confSpace.positions) {
			for (SimpleConfSpace.ResidueConf rc1 : pos1.resConfs) {
				for (String resNum : confSpace.shellResNumbers) {
					set(pos1, rc1, resNum, read(in, new RCTuple(pos1.index, rc1.index)));
				}
			}
		}
	}
}
