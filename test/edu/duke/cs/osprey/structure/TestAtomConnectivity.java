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
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPair;
import edu.duke.cs.osprey.structure.AtomConnectivity.AtomPairList;


public class TestAtomConnectivity {
	
	private static Strand strand;
	private static SimpleConfSpace confSpace;
	private static AtomConnectivity connectivity;
	
	@BeforeClass
	public static void beforeClass() {
		
		// get an arbitrary conf space
		strand = new Strand.Builder(PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb")).build();
		for (Residue res : strand.mol.residues) {
			strand.flexibility.get(res.getPDBResNumber()).addWildTypeRotamers();
		}
		confSpace = new SimpleConfSpace.Builder().addStrand(strand).build();
		
		connectivity = new AtomConnectivity.Builder()
			.setConfSpace(confSpace)
			.setParallelism(Parallelism.makeCpu(4))
			.build();
	}
	
	@Test
	public void pairDirectionsSame() {
		
		// the same adjacent residues in different orders should get the same atom pairs
		
		// VAL->LEU pair
		Residue val25 = strand.mol.getResByPDBResNumber("25");
		assertThat(val25.template.name, is("VAL"));
		Residue leu26 = strand.mol.getResByPDBResNumber("26");
		assertThat(leu26.template.name, is("LEU"));
		AtomPairList pairs1 = connectivity.getAtomPairs(strand.mol, val25, leu26);
		AtomPairList pairs2 = connectivity.getAtomPairs(strand.mol, leu26, val25);
		assertThat(pairs1 == pairs2, is(true));
	}
	
	@Test
	public void pairDirectionsDifferent() {
		
		// different adjacent residues should get different atom pairs
		
		// VAL->LEU pair
		Residue val25 = strand.mol.getResByPDBResNumber("25");
		assertThat(val25.template.name, is("VAL"));
		Residue leu26 = strand.mol.getResByPDBResNumber("26");
		assertThat(leu26.template.name, is("LEU"));
		AtomPairList pairs1 = connectivity.getAtomPairs(strand.mol, val25, leu26);
		
		// LEU->VAL pair
		Residue leu44 = strand.mol.getResByPDBResNumber("44");
		assertThat(leu44.template.name, is("LEU"));
		Residue val45 = strand.mol.getResByPDBResNumber("45");
		assertThat(val45.template.name, is("VAL"));
		AtomPairList pairs2 = connectivity.getAtomPairs(strand.mol, leu44, val45);
		
		assertThat(pairs1 == pairs2, is(false));
	}
	
	@Test
	public void separateDirections() {
		
		// non-adjacent residues in different orders should get the same atom pairs
		
		// VAL->LEU pair
		Residue val25 = strand.mol.getResByPDBResNumber("25");
		assertThat(val25.template.name, is("VAL"));
		Residue gln43 = strand.mol.getResByPDBResNumber("43");
		assertThat(gln43.template.name, is("GLN"));
		AtomPairList pairs1 = connectivity.getAtomPairs(strand.mol, val25, gln43);
		AtomPairList pairs2 = connectivity.getAtomPairs(strand.mol, gln43, val25);
		
		assertThat(pairs1 == pairs2, is(true));
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
		
		AtomPairList pairs = connectivity.getAtomPairs(strand.mol, res1, res2);
		
		if (pairs == null) {
			fail("no pairs for residue types: " + res1.template.name + ", " + res2.template.name);
		}
		
		for (AtomNeighbors.Type neighborType : AtomNeighbors.Type.values()) {
			
			List<int[]> typedPairs = AtomNeighbors.getPairIndicesByType(res1.atoms, res2.atoms, res1 == res2, neighborType);
			
			String desc = String.format("%d:%-14s - %d:%-14s - %s --> %s - %s",
				res1.indexInMolecule, res1.template,
				res2.indexInMolecule, res2.template,
				neighborType,
				pairs.templa, pairs.templb
			);
			
			// check that the pairs from the connectivity lookup match the brute force search
			assertAtomPairs(desc, res1, res2, neighborType, pairs, typedPairs);
		}
	}

	private void assertAtomPairs(String desc, Residue res1, Residue res2, AtomNeighbors.Type neighborType, AtomPairList pairs, List<int[]> typedPairs) {
		
		// count the pairs of this type
		long numPairsOfType = pairs.pairs.stream()
			.filter((AtomPair pair) -> pair.type == neighborType)
			.count();
		
		// try to match each typed pair
		long matched = 0;
		for (int[] typedPair : typedPairs) {
			assertThat(
				desc + "\natom pair not found: " + Arrays.toString(typedPair) + "\nin: " + pairs.pairs,
				findAtomPair(res1, res2, pairs, typedPair), is(true)
			);
			matched++;
		}
		
		assertThat(desc, matched, is(numPairsOfType));
	}

	private boolean findAtomPair(Residue res1, Residue res2, AtomPairList pairs, int[] typedPair) {
		for (int i=0; i<pairs.size(); i++) {
			if (pairs.getIndex1(res1, res2, i) == typedPair[0] && pairs.getIndex2(res1, res2, i) == typedPair[1]) {
				return true;
			}
		}
		return false;
	}
}
