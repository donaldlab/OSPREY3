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

import static edu.duke.cs.osprey.tools.Log.log;
import static org.hamcrest.Matchers.*;
import static org.junit.Assert.*;

import edu.duke.cs.osprey.confspace.ParametricMolecule;
import edu.duke.cs.osprey.confspace.RCTuple;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.energy.ResInterGen;
import edu.duke.cs.osprey.energy.ResidueInteractions;
import edu.duke.cs.osprey.structure.*;
import org.junit.Test;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;


public class TestPLUG {

	private static final double tolerance = 0.4;

	public static void dumpContacts(SimpleConfSpace confSpace, RCTuple tuple) {

		Probe probe = new Probe();
		probe.matchTemplates(confSpace);

		AtomConnectivity connectivity = new AtomConnectivity.Builder()
			.addTemplates(confSpace)
			.build();

		// make full inters
		ResidueInteractions inters = ResInterGen.of(confSpace)
			.addIntras(tuple)
			.addInters(tuple)
			.addShell(tuple)
			.make();

		ParametricMolecule pmol = confSpace.makeMolecule(tuple);

		// get the contacts sorted by violation
		List<Probe.AtomPair.Interaction> contacts = probe.getInteractions(pmol.mol.residues, inters, connectivity).stream()
			.filter(interaction -> interaction.contact.isContact)
			.sorted(Comparator.comparingDouble(contact -> contact.getViolation(tolerance)))
			.collect(Collectors.toList());

		for (Probe.AtomPair.Interaction contact : contacts) {
			log("v=%5.2f %s - %s", contact.getViolation(tolerance), contact.atomPair, contact);
		}
	}

	public static void writeViolations2D(SimpleConfSpace confSpace, RCTuple tuple, String resNum1, String atomName1, String resNum2, String atomName2) {

		ParametricMolecule pmol = confSpace.makeMolecule(tuple);

		Residue res1 = pmol.mol.residues.getOrThrow(resNum1);
		Atom atom1 = res1.getAtomByName(atomName1);
		Residue res2 = pmol.mol.residues.getOrThrow(resNum2);
		Atom atom2 = res2.getAtomByName(atomName2);

		Probe probe = new Probe();
		probe.matchTemplates(confSpace);
		Probe.AtomPair pair = probe.new AtomPair(atom1, atom2);

		File file = new File("PLUG.vtk");
		PLUGLab.writeViolations2D(pmol, 0, 1, pair, tolerance, file);
	}

	@Test
	public void singleLeucine() {

		Strand strand = new Strand.Builder(PDBIO.readResource("/1CC8.ss.pdb")).build();
		strand.flexibility.get("A23").setLibraryRotamers("LEU").setContinuous();
		SimpleConfSpace confSpace = new SimpleConfSpace.Builder()
			.addStrand(strand)
			.build();

		PLUG plug = new PLUG(confSpace);

		RCTuple tuple;

		tuple = new RCTuple(0, 0);
		//dumpContacts(confSpace, tuple);
		//writeViolations2D(confSpace, new RCTuple(0, 0), "A23", "HD12", "A24", "N");
		// Unavoidable: v= 0.99 A23:HD12 <-> A24:N    -            overlap=   1.391  >0.000 BadClash
		assertThat(plug.shouldPruneTuple(tuple, tolerance), is(true));

		tuple = new RCTuple(0, 1);
		//dumpContacts(confSpace, tuple);
		//writeViolations2D(confSpace, new RCTuple(0, 1), "A23", "HD13", "A36", "HD11");
		// Unavoidable: v= 0.67 A23:HD13 <-> A36:HD11 -            overlap=   1.074  >0.000 BadClash
		assertThat(plug.shouldPruneTuple(tuple, tolerance), is(true));

		tuple = new RCTuple(0, 2);
		//dumpContacts(confSpace, tuple);
		//writeViolations2D(confSpace, tuple, "A23", "C", "A23", "HD23");
		// Avoidable: v= 0.10 A23:C    <-> A23:HD23 -            overlap=   0.499  >0.000 BadClash
		assertThat(plug.shouldPruneTuple(tuple, tolerance), is(false));

		tuple = new RCTuple(0, 3);
		//dumpContacts(confSpace, tuple);
		//writeViolations2D(confSpace, tuple, "A23", "HD11", "A38", "HD11");
		// Unavoidable: v= 1.24 A23:HD11 <-> A38:HD11 -            overlap=   1.643  >0.000 BadClash
		assertThat(plug.shouldPruneTuple(tuple, tolerance), is(true));

		tuple = new RCTuple(0, 4);
		//dumpContacts(confSpace, tuple);
		//writeViolations2D(confSpace, tuple, "A23", "CD2", "A36", "HD11");
		// Unavoidable: v= 1.13 A23:CD2  <-> A36:HD11 -            overlap=   1.526  >0.000 BadClash
		assertThat(plug.shouldPruneTuple(tuple, tolerance), is(true));
	}
}
