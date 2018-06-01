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

import org.junit.BeforeClass;
import org.junit.Test;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.confspace.Strand.ResidueFlex;

public class TestStrand {
	
	private static Molecule mol;
	
	@BeforeClass
	public static void beforeClass() {
		mol = PDBIO.readFile("examples/1CC8/1CC8.ss.pdb");
	}
	
	@Test
	public void fullResidues() {
		
		Strand strand = new Strand.Builder(mol)
			.setErrorOnNonTemplateResidues(true)
			.build();
		
		assertThat(strand.mol.residues.get(0).getPDBResNumber(), is("A2"));
		assertThat(strand.mol.residues.get(strand.mol.residues.size() - 1).getPDBResNumber(), is("A73"));
	}
	
	@Test
	public void subsequenceResidues() {
		
		Strand strand = new Strand.Builder(mol)
			.setResidues("A5", "A70")
			.build();
		
		assertThat(strand.mol.residues.get(0).getPDBResNumber(), is("A5"));
		assertThat(strand.mol.residues.get(strand.mol.residues.size() - 1).getPDBResNumber(), is("A70"));
	}
	
	@Test
	public void structureReferences() {
		
		Strand strand = new Strand.Builder(mol).build();
		
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
	
	@Test
	public void defaultFlexibilty() {
		
		Strand strand = new Strand.Builder(mol)
			.setResidues("A2", "A5")
			.build();
		
		for (Residue res : strand.mol.residues) {
			ResidueFlex resFlex = strand.flexibility.get(res.getPDBResNumber());
			assertThat(resFlex.isFlexible(), is(false));
		}
		
		assertThat(strand.flexibility.getFlexibleResidueNumbers().isEmpty(), is(true));
		assertThat(strand.flexibility.getStaticResidueNumbers(), contains("A2", "A3", "A4", "A5"));
	}
	
	@Test
	public void flexibleResidueNumbers() {
		
		Strand strand = new Strand.Builder(mol)
			.setResidues("A2", "A5")
			.build();
		
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		
		assertThat(strand.flexibility.getFlexibleResidueNumbers(), contains("A2", "A4"));
	}
	
	@Test
	public void staticResidueNumbers() {
		
		Strand strand = new Strand.Builder(mol)
			.setResidues("A2", "A5")
			.build();
		
		strand.flexibility.get("A2").setLibraryRotamers(Strand.WildType);
		strand.flexibility.get("A4").setLibraryRotamers(Strand.WildType);
		
		assertThat(strand.flexibility.getStaticResidueNumbers(), contains("A3", "A5"));
	}
}
