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
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.energy.approximation.ApproximatedObjectiveFunction.Approximator;
import edu.duke.cs.osprey.tools.IOable;

import java.io.*;
import java.util.*;


/**
 * Forcefield approximators for RC,residue pairs in a conf space
 */
public class ApproximatorMatrix implements IOable {

	public static class Entry {

		public final String resNum;
		public final Approximator.Addable approximator;

		public Entry(String resNum, Approximator.Addable approximator) {
			this.resNum = resNum;
			this.approximator = approximator;
		}
	}

	public final SimpleConfSpace confSpace;

	private final int[] offsets;
	private final List<List<Entry>> entries;

	public ApproximatorMatrix(SimpleConfSpace confSpace) {

		this.confSpace = confSpace;

		// compute the offsets
		int offset = 0;
		offsets = new int[confSpace.positions.size()];
		for (SimpleConfSpace.Position pos : confSpace.positions) {
			offsets[pos.index] = offset;
			offset += pos.resConfs.size();
		}

		entries = new ArrayList<>(offset);
		for (int i=0; i<offset; i++) {
			entries.add(new ArrayList<>());
		}
	}

	private int getIndex(int pos, int rc) {
		return offsets[pos] + rc;
	}

	public List<Entry> get(int pos, int rc) {
		return Collections.unmodifiableList(entries.get(getIndex(pos, rc)));
	}

	public Approximator.Addable get(int pos, int rc, String resNum) {
		for (Entry entry : entries.get(getIndex(pos, rc))) {
			if (entry.resNum.equals(resNum)) {
				return entry.approximator;
			}
		}
		return null;
	}

	public void set(int pos, int rc, String resNum, Approximator.Addable approximator) {

		List<Entry> entries = this.entries.get(getIndex(pos, rc));

		// if the list already has an entry for this res num, remove it
		for (int i=0; i<entries.size(); i++) {
			if (entries.get(i).resNum.equals(resNum)) {
				entries.remove(i);
				break;
			}
		}

		// add the new approximator in order of weakly increasing error
		int i = 0;
		for (; i<entries.size(); i++) {
			if (approximator.error() < entries.get(i).approximator.error()) {
				break;
			}
		}
		entries.add(i, new Entry(resNum, approximator));
	}

	public ResidueInteractionsApproximator get(RCTuple tuple, ResidueInteractions inters, double errorBudget) {

		double errorBudgetPerResidue = errorBudget/tuple.size();

		Approximator[] approximators = new Approximator[tuple.size()];
		ResidueInteractions approxInters = new ResidueInteractions();

		for (int i=0; i<tuple.size(); i++) {

			int pos = tuple.pos.get(i);
			int rc = tuple.RCs.get(i);

			String resNum = confSpace.positions.get(pos).resNum;

			Approximator.Addable approximator = null;
			for (Entry entry : entries.get(getIndex(pos, rc))) {
				ResidueInteractions.Pair pair = inters.get(resNum, entry.resNum);
				if (pair != null) {

					// check the budget
					double error = approximator != null ? approximator.error() : 0.0;
					if (error + entry.approximator.error() > errorBudgetPerResidue) {
						break;
					}

					approxInters.add(pair);

					if (approximator == null) {
						approximator = entry.approximator.makeIdentity();
					}
					approximator.add(entry.approximator, pair.weight, pair.offset);
				}
			}

			// didn't find anything within budget? make a no-op approximator
			if (approximator == null) {
				int d = confSpace.countDofs(new RCTuple(pos, rc));
				approximators[i] = new NOPApproximator(d);
			} else {
				approximators[i] = approximator;
			}
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
		if (approximator instanceof QuadraticApproximator) {
			return 1;
		} else {
			throw new IllegalArgumentException("unrecognized approximator type: " + approximator.getClass().getName());
		}
	}

	private Approximator.Addable alloc(byte type, int d) {
		switch (type) {
			case 1: return new QuadraticApproximator(d);
			default: throw new IllegalArgumentException("unrecognized approximator type: " + type);
		}
	}

	@Override
	public void writeTo(DataOutput out)
	throws IOException {

		for (List<Entry> entries : this.entries) {
			out.writeInt(entries.size());
			for (Entry entry : entries) {
				out.writeUTF(entry.resNum);
				out.writeByte(getType(entry.approximator));
				out.writeInt(entry.approximator.numDofs());
				IOable.of(entry.approximator).writeTo(out);
			}
		}
	}

	@Override
	public void readFrom(DataInput in)
	throws IOException {

		for (List<Entry> entries : this.entries) {

			entries.clear();
			int numEntries = in.readInt();
			for (int i=0; i<numEntries; i++) {

				String resNum = in.readUTF();
				byte type = in.readByte();
				int d = in.readInt();
				Approximator.Addable approximator = alloc(type, d);
				IOable.of(approximator).readFrom(in);

				entries.add(new Entry(resNum, approximator));
			}
		}
	}
}
