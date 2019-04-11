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

package edu.duke.cs.osprey.pruning;


import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.parallelism.TaskExecutor;
import edu.duke.cs.osprey.tools.Progress;

import java.util.concurrent.atomic.AtomicLong;


/**
 * Checks if each tuple can possibly be part of a full unpruned conformation.
 *
 * This pruning essentially checks for the same thing as the LUTE conf sampler,
 * (namely, can we map a tuple into a set of confs)
 * but using a polynomial time algorithm instead of an exponential time one
 */
public class TransitivePruner {

	public final SimpleConfSpace confSpace;

	public TransitivePruner(SimpleConfSpace confSpace) {
		this.confSpace = confSpace;
	}

	public void pruneSingles(PruningMatrix pmat, TaskExecutor tasks) {

		// this can take a while, so track progress
		AtomicLong numSingles = new AtomicLong(0);
		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			numSingles.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numSingles.get());

		pmat.forEachUnprunedSingle((pos1, rc1) -> {
			tasks.submit(
				() -> isPrunedTransitively(pmat, pos1, rc1),
				(isPruned) -> {
					if (isPruned) {
						pmat.pruneSingle(pos1, rc1);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	private boolean isPrunedTransitively(PruningMatrix pmat, int pos1, int rc1) {

		for (int pos2=0; pos2<pmat.getNumPos(); pos2++) {

			// skip assigned positions
			if (pos2 == pos1) {
				continue;
			}

			// check pairs
			if (!isPos2Assignable(pmat, pos1, rc1, pos2)) {
				return true;
			}

			for (int pos3=0; pos3<pos2; pos3++) {

				// skip assigned positions
				if (pos3 == pos1) {
					continue;
				}

				// check triples
				if (!arePos23Assignable(pmat, pos1, rc1, pos2, pos3)) {
					return true;
				}
			}
		}

		return false;
	}

	private boolean isPos2Assignable(PruningMatrix pmat, int pos1, int rc1, int pos2) {

		// is there at least one possible assignment to pos2, given rc1?
		for (int rc2=0; rc2<pmat.getNumConfAtPos(pos2); rc2++) {
			if (!pmat.isPairPruned(pos1, rc1, pos2, rc2)) {
				return true;
			}
		}

		return false;
	}

	private boolean arePos23Assignable(PruningMatrix pmat, int pos1, int rc1, int pos2, int pos3) {

		// is there at least one possible assignment to pos2,pos3, given rc1?
		for (int rc2=0; rc2<pmat.getNumConfAtPos(pos2); rc2++) {
			for (int rc3=0; rc3<pmat.getNumConfAtPos(pos3); rc3++) {
				if (!pmat.isTriplePruned(pos1, rc1, pos2, rc2, pos3, rc3)) {
					return true;
				}
			}
		}

		return false;
	}

	public void prunePairs(PruningMatrix pmat, TaskExecutor tasks) {

		// this can take a while, so track progress
		AtomicLong numPairs = new AtomicLong(0);
		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			numPairs.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numPairs.get());

		pmat.forEachUnprunedPair((pos1, rc1, pos2, rc2) -> {
			tasks.submit(
				() -> isPrunedTransitively(pmat, pos1, rc1, pos2, rc2),
				(isPruned) -> {
					if (isPruned) {
						pmat.prunePair(pos1, rc1, pos2, rc2);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	private boolean isPrunedTransitively(PruningMatrix pmat, int pos1, int rc1, int pos2, int rc2) {

		for (int pos3=0; pos3<pmat.getNumPos(); pos3++) {

			// skip assigned positions
			if (pos3 == pos1 || pos3 == pos2) {
				continue;
			}

			// check pairs
			if (!isPos3Assignable(pmat, pos1, rc1, pos2, rc2, pos3)) {
				return true;
			}

			for (int pos4=0; pos4<pos3; pos4++) {

				// skip assigned positions
				if (pos4 == pos1 || pos4 == pos2) {
					continue;
				}

				// check triples
				if (!arePos34Assignable(pmat, pos1, rc1, pos2, rc2, pos3, pos4)) {
					return true;
				}
			}
		}

		return false;
	}

	private boolean isPos3Assignable(PruningMatrix pmat, int pos1, int rc1, int pos2, int rc2, int pos3) {

		// is there at least one possible assignment to pos3, given rc1,rc2?
		for (int rc3=0; rc3<pmat.getNumConfAtPos(pos3); rc3++) {
			if (!pmat.isTriplePruned(pos1, rc1, pos2, rc2, pos3, rc3)) {
				return true;
			}
		}

		return false;
	}

	private boolean arePos34Assignable(PruningMatrix pmat, int pos1, int rc1, int pos2, int rc2, int pos3, int pos4) {

		// is there at least one possible assignment to pos3,pos4, given rc1,rc2?
		for (int rc3=0; rc3<pmat.getNumConfAtPos(pos3); rc3++) {
			for (int rc4=0; rc4<pmat.getNumConfAtPos(pos4); rc4++) {
				if (!pmat.isQuadruplePruned(pos1, rc1, pos2, rc2, pos3, rc3, pos4, rc4)) {
					return true;
				}
			}
		}

		return false;
	}

	public void pruneTriples(PruningMatrix pmat, TaskExecutor tasks) {

		// this can take a while, so track progress
		AtomicLong numTriples = new AtomicLong(0);
		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			numTriples.incrementAndGet();
			return PruningMatrix.IteratorCommand.Continue;
		});
		Progress progress = new Progress(numTriples.get());

		pmat.forEachUnprunedTriple((pos1, rc1, pos2, rc2, pos3, rc3) -> {
			tasks.submit(
				() -> isPrunedTransitively(pmat, pos1, rc1, pos2, rc2, pos3, rc3),
				(isPruned) -> {
					if (isPruned) {
						pmat.pruneTriple(pos1, rc1, pos2, rc2, pos3, rc3);
					}
					progress.incrementProgress();
				}
			);
			return PruningMatrix.IteratorCommand.Continue;
		});
	}

	private boolean isPrunedTransitively(PruningMatrix pmat, int pos1, int rc1, int pos2, int rc2, int pos3, int rc3) {

		for (int pos4=0; pos4<pmat.getNumPos(); pos4++) {

			// skip assigned positions
			if (pos4 == pos1 || pos4 == pos2 || pos4 == pos3) {
				continue;
			}

			// check pairs
			if (!isPos4Assignable(pmat, pos1, rc1, pos2, rc2, pos3, rc3, pos4)) {
				return true;
			}

			for (int pos5=0; pos5<pos4; pos5++) {

				// skip assigned positions
				if (pos5 == pos1 || pos5 == pos2 || pos5 == pos3) {
					continue;
				}

				// check triples
				if (!arePos45Assignable(pmat, pos1, rc1, pos2, rc2, pos3, rc3, pos4, pos5)) {
					return true;
				}
			}
		}

		return false;
	}

	private boolean isPos4Assignable(PruningMatrix pmat, int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4) {

		// is there at least one possible assignment to pos4, given rc1,rc2,rc3?
		for (int rc4=0; rc4<pmat.getNumConfAtPos(pos4); rc4++) {
			if (!pmat.isQuadruplePruned(pos1, rc1, pos2, rc2, pos3, rc3, pos4, rc4)) {
				return true;
			}
		}

		return false;
	}

	private boolean arePos45Assignable(PruningMatrix pmat, int pos1, int rc1, int pos2, int rc2, int pos3, int rc3, int pos4, int pos5) {

		// is there at least one possible assignment to pos4,pos5, given rc1,rc2,rc3?
		for (int rc4=0; rc4<pmat.getNumConfAtPos(pos4); rc4++) {
			for (int rc5=0; rc5<pmat.getNumConfAtPos(pos5); rc5++) {
				if (!pmat.isQuintuplePruned(pos1, rc1, pos2, rc2, pos3, rc3, pos4, rc4, pos5, rc5)) {
					return true;
				}
			}
		}

		return false;
	}
}
