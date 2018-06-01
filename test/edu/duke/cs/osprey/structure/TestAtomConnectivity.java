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
import static org.junit.Assert.*;

import java.util.Arrays;
import java.util.List;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.parallelism.Parallelism;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPairs;


public class TestAtomConnectivity {
	
	private static Strand strand;
	private static SimpleConfSpace confSpace;
	private static AtomConnectivity connectivity;
	
	@BeforeClass
	public static void beforeClass() {
		
		// get an arbitrary conf space
		strand = new Strand.Builder(PDBIO.readFile("examples/python.GMEC/1CC8.ss.pdb")).build();
		for (Residue res : strand.mol.residues) {
			strand.flexibility.get(res.getPDBResNumber()).addWildTypeRotamers();
		}
		confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		connectivity = new AtomConnectivity.Builder()
			.addTemplates(confSpace)
			.setParallelism(Parallelism.makeCpu(1))
			.build();
	}
	
	@Test
	public void allSelfPairs() {
		
		for (int i=0; i<strand.mol.residues.size(); i++) {
			Residue res1 = strand.mol.residues.get(i);
			assertResiduePair(res1, res1);
		}
	}
	
	@Test
	public void allAdjacentPairsForward() {
		
		for (int i=0; i<strand.mol.residues.size(); i++) {
			Residue res1 = strand.mol.residues.get(i);
			
			if (i > 0) {
				Residue res2 = strand.mol.residues.get(i - 1);
				assertResiduePair(res2, res1);
			} else if (i < strand.mol.residues.size() - 1) {
				Residue res2 = strand.mol.residues.get(i + 1);
				assertResiduePair(res1, res2);
			}
		}
	}

	@Test
	public void allAdjacentPairsReverse() {
		
		for (int i=0; i<strand.mol.residues.size(); i++) {
			Residue res1 = strand.mol.residues.get(i);
			
			if (i > 0) {
				Residue res2 = strand.mol.residues.get(i - 1);
				assertResiduePair(res1, res2);
			} else if (i < strand.mol.residues.size() - 1) {
				Residue res2 = strand.mol.residues.get(i + 1);
				assertResiduePair(res2, res1);
			}
		}
	}
	
	// this test takes a really long time (~10 minutes)
	// you probably only want to run it if you're tracking a specific bug
	// NOPE: it takes about 30 seconds now after the latest optimizations, but meh
	//@Test
	public void allSeparatePairs() {
		
		for (int i=0; i<strand.mol.residues.size(); i++) {
			Residue res1 = strand.mol.residues.get(i);
			for (int j=0; j<strand.mol.residues.size(); j++) {
				
				if (Math.abs(i - j) <= 1) {
					continue;
				}
				
				Residue res2 = strand.mol.residues.get(j);
				
				assertResiduePair(res1, res2);
			}
		}
	}
	
	@Test
	public void someSeparatePairs() {
		
		int a = 5;
		int b = 12;
		
		for (int i=a; i<b; i++) {
			Residue res1 = strand.mol.residues.get(i);
			for (int j=a; j<b; j++) {
				
				if (Math.abs(i - j) <= 1) {
					continue;
				}
				
				Residue res2 = strand.mol.residues.get(j);
				
				assertResiduePair(res1, res2);
			}
		}
	}
	
	private void assertResiduePair(Residue res1, Residue res2) {
		
		AtomPairs pairs = connectivity.getAtomPairs(res1, res2);
		
		if (pairs == null) {
			fail("no pairs for residue types: " + res1.template.name + ", " + res2.template.name);
		}
		
		assertThat(pairs.res1.template, is(res1.template));
		assertThat(pairs.res2.template, is(res2.template));
		
		for (AtomNeighbors.Type neighborType : AtomNeighbors.Type.values()) {
			
			List<int[]> typedPairs = AtomNeighbors.getPairIndicesByType(res1.atoms, res2.atoms, res1 == res2, neighborType);
			
			String desc = String.format("%d:%-14s - %d:%-14s - %s",
				res1.indexInMolecule, res1.template,
				res2.indexInMolecule, res2.template,
				neighborType
			);
			
			// check that the pairs from the connectivity lookup match the brute force search
			assertAtomPairs(desc, res1, res2, neighborType, pairs, typedPairs);
		}
	}

	private void assertAtomPairs(String desc, Residue res1, Residue res2, AtomNeighbors.Type neighborType, AtomPairs pairs, List<int[]> typedPairs) {
		
		// try to match each typed pair
		int matched = 0;
		for (int[] typedPair : typedPairs) {
			assertThat(
				desc + "\natom pair not found: " + Arrays.toString(typedPair) + "\nin: " + pairs.dumpPairs(neighborType),
				findAtomPair(res1, res2, pairs, neighborType, typedPair),
				is(true)
			);
			matched++;
		}
		
		assertThat(desc, matched, is(pairs.getNumPairs(neighborType)));
	}

	private boolean findAtomPair(Residue res1, Residue res2, AtomPairs pairs, AtomNeighbors.Type neighborType, int[] typedPair) {
		for (int[] pair : pairs.getPairs(neighborType)) {
			if (pair[0] == typedPair[0] && pair[1] == typedPair[1]) {
				return true;
			}
		}
		return false;
	}
}
