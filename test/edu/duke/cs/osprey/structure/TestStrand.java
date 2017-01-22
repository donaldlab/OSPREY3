package edu.duke.cs.osprey.structure;

import static org.junit.Assert.*;
import static org.hamcrest.Matchers.*;

import org.junit.Test;

import edu.duke.cs.osprey.confspace.Strand;

public class TestStrand {
	
	@Test
	public void full() {
		
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		Strand strand = Strand.builder(mol)
			.setErrorOnNonTemplateResidues(true)
			.build();
		
		assertThat(strand.mol.residues.get(0).getPDBResNumber(), is("2"));
		assertThat(strand.mol.residues.get(strand.mol.residues.size() - 1).getPDBResNumber(), is("73"));
	}
	
	@Test
	public void subsequence() {
		
		Molecule mol = PDBIO.readFile("examples/1CC8.python/1CC8.ss.pdb");
		Strand strand = Strand.builder(mol)
			.setResidues(5, 70)
			.setErrorOnNonTemplateResidues(true)
			.build();
		
		assertThat(strand.mol.residues.get(0).getPDBResNumber(), is("5"));
		assertThat(strand.mol.residues.get(strand.mol.residues.size() - 1).getPDBResNumber(), is("70"));
	}
}
