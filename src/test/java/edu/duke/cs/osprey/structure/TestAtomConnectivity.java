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

package edu.duke.cs.osprey.structure;

import static org.hamcrest.Matchers.*;
import static org.hamcrest.MatcherAssert.*;
import static org.junit.jupiter.api.Assertions.fail;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPairs;


public class TestAtomConnectivity {
	
	private static class Case {

		Strand strand;
		SimpleConfSpace confSpace;
		AtomConnectivity connectivity;

		Case(Strand strand) {
			this.strand = strand;
			this.confSpace = new SimpleConfSpace.Builder()
				.addStrand(strand)
				.build();
			this.connectivity = new AtomConnectivity.Builder()
				.set15HasNonBonded(true)
				.build();
		}

		List<Residue> flexibleResidues() {
			return strand.flexibility.getFlexibleResidueNumbers().stream()
				.map(resNum -> strand.mol.residues.getOrThrow(resNum))
				.collect(Collectors.toList());
		}

		void assertResiduePair(Residue res1, Residue res2) {

			AtomPairs pairs = connectivity.getAtomPairs(res1, res2);

			if (pairs == null) {
				fail("no pairs for residue types: " + res1.template.name + ", " + res2.template.name);
			}

			assertThat(pairs.templ1, is(res1.template));
			assertThat(pairs.templ2, is(res2.template));

			for (AtomNeighbors.Type neighborType : AtomNeighbors.Type.values()) {

				List<int[]> typedPairs = AtomNeighbors.getPairIndicesByType(res1.atoms, res2.atoms, res1 == res2, neighborType);

				String desc = String.format("%d:%-14s - %d:%-14s - %s",
					res1.indexInMolecule, res1.template,
					res2.indexInMolecule, res2.template,
					neighborType
				);

				// check that the pairs from the connectivity lookup match the brute force search
				assertAtomPairs(desc, neighborType, pairs, typedPairs);
			}
		}

		void assertAtomPairs(String desc, AtomNeighbors.Type neighborType, AtomPairs pairs, List<int[]> typedPairs) {

			// try to match each typed pair
			int matched = 0;
			for (int[] typedPair : typedPairs) {
				assertThat(
					desc + "\natom pair not found: " + Arrays.toString(typedPair) + "\nin: " + pairs.dumpPairs(neighborType),
					findAtomPair(pairs, neighborType, typedPair),
					is(true)
				);
				matched++;
			}

			assertThat(desc, matched, is(pairs.getNumPairs(neighborType)));
		}

		boolean findAtomPair(AtomPairs pairs, AtomNeighbors.Type neighborType, int[] typedPair) {
			for (int[] pair : pairs.getPairs(neighborType)) {
				if (pair[0] == typedPair[0] && pair[1] == typedPair[1]) {
					return true;
				}
			}
			return false;
		}
	}

	private static Strand make1CC8() {
		return new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
	}

	private static Strand make1CC8Cyx() {
		return new Strand.Builder(PDBIO.readResource("/1CC8.CYX.pdb")).build();
	}

	private static Strand setAllFlexible(Strand strand) {
		for (Residue res : strand.mol.residues) {
			strand.flexibility.get(res.getPDBResNumber()).addWildTypeRotamers();
		}
		return strand;
	}

	private static Strand setNonCyxFlexible(Strand strand) {
		for (Residue res : strand.mol.residues) {
			if (res.template.name.equals("CYX")) {
				continue;
			}
			strand.flexibility.get(res.getPDBResNumber()).addWildTypeRotamers();
		}
		return strand;
	}

	private static void allSelfPairs(Case c) {

		List<Residue> residues = c.flexibleResidues();

		for (int i=0; i<residues.size(); i++) {
			Residue res1 = c.strand.mol.residues.get(i);
			c.assertResiduePair(res1, res1);
		}
	}
	@Test public void allSelfPairs_1CC8_allFlexible() { allSelfPairs(new Case(setAllFlexible(make1CC8()))); }
	@Test public void allSelfPairs_1CC8Cyx_nonCyxFlexible() { allSelfPairs(new Case(setNonCyxFlexible(make1CC8Cyx()))); }

	private static void allAdjacentPairsForward(Case c) {

		List<Residue> residues = c.flexibleResidues();

		for (int i=0; i<residues.size(); i++) {
			Residue res1 = residues.get(i);
			
			if (i > 0) {
				Residue res2 = residues.get(i - 1);
				c.assertResiduePair(res2, res1);
			} else if (i < residues.size() - 1) {
				Residue res2 = residues.get(i + 1);
				c.assertResiduePair(res1, res2);
			}
		}
	}
	@Test public void allAdjacentPairsForward_1CC8_allFlexible() { allAdjacentPairsForward(new Case(setAllFlexible(make1CC8()))); }
	@Test public void allAdjacentPairsForward_1CC8Cyx_nonCyxFlexible() { allAdjacentPairsForward(new Case(setNonCyxFlexible(make1CC8Cyx()))); }

	private static void allAdjacentPairsReverse(Case c) {

		List<Residue> residues = c.flexibleResidues();

		for (int i=0; i<residues.size(); i++) {
			Residue res1 = residues.get(i);
			
			if (i > 0) {
				Residue res2 = residues.get(i - 1);
				c.assertResiduePair(res1, res2);
			} else if (i < residues.size() - 1) {
				Residue res2 = residues.get(i + 1);
				c.assertResiduePair(res2, res1);
			}
		}
	}
	@Test public void allAdjacentPairsReverse_1CC8_allFlexible() { allAdjacentPairsReverse(new Case(setAllFlexible(make1CC8()))); }
	@Test public void allAdjacentPairsReverse_1CC8Cyx_nonCyxFlexible() { allAdjacentPairsReverse(new Case(setNonCyxFlexible(make1CC8Cyx()))); }

	private static void allSeparatePairs(Case c) {

		List<Residue> residues = c.flexibleResidues();

		for (int i=0; i<residues.size(); i++) {
			Residue res1 = residues.get(i);
			for (int j=0; j<residues.size(); j++) {
				
				if (Math.abs(i - j) <= 1) {
					continue;
				}
				
				Residue res2 = residues.get(j);
				
				c.assertResiduePair(res1, res2);
			}
		}
	}
	@Test public void allSeparatePairs_1CC8_allFlexible() { allSeparatePairs(new Case(setAllFlexible(make1CC8()))); }
	@Test public void allSeparatePairs_1CC8Cyx_nonCyxFlexible() { allSeparatePairs(new Case(setNonCyxFlexible(make1CC8Cyx()))); }

	// flexibiility across disulfide bonds is not yet supported, if the user tries this, they should see an error
	@Test
	public void allSeparatePairs_1CC8Cyx_allFlexible() {
		Assertions.assertThrows(UnsupportedOperationException.class, () -> {
			allSeparatePairs(new Case(setAllFlexible(make1CC8Cyx())));
		});
	}
}
