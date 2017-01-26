package edu.duke.cs.osprey.structure;

import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.confspace.Strand;

public class TestStrand {
	
	private static Molecule mol;
	
	@BeforeClass
	public static void beforeClass() {
		mol = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
	}
	
	@Test
	public void fullResidues() {
		
		Strand strand = Strand.builder(mol)
			.setErrorOnNonTemplateResidues(true)
			.build();
		
		assertThat(strand.mol.residues.get(0).getPDBResNumber(), is("2"));
		assertThat(strand.mol.residues.get(strand.mol.residues.size() - 1).getPDBResNumber(), is("73"));
	}
	
	@Test
	public void subsequenceResidues() {
		
		Strand strand = Strand.builder(mol)
			.setResidues(5, 70)
			.build();
		
		assertThat(strand.mol.residues.get(0).getPDBResNumber(), is("5"));
		assertThat(strand.mol.residues.get(strand.mol.residues.size() - 1).getPDBResNumber(), is("70"));
	}
	
	@Test
	public void structureReferences() {
		
		Strand strand = Strand.builder(mol).build();
		
		// make sure various structure the references are all set correctly
		for (int i=0; i<strand.mol.residues.size(); i++) {
			
			Residue res = strand.mol.residues.get(i);
			assert (res.molec == strand.mol);
			
			// make sure all the bonds are marked
			assert (res.intraResBondsMarked);
			assert (res.interResBondsMarked);
			
			// make sure all the atoms point to the right residues
			for (Atom atom : res.atoms) {
				assert (atom.res == res);
			}
			
			// check the alternates too
			for (Residue altRes : strand.mol.getAlternates(i)) {
				for (Atom atom : altRes.atoms) {
					assert (atom.res == altRes);
				}
			}
			
			// make sure every residue has a template
			assert (res.template != null);
			for (Residue altRes : strand.mol.getAlternates(i)) {
				assert (altRes.template != null);
			}
		}
	}
}
