package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.SharedSpec
import edu.duke.cs.osprey.gui.io.fromPDB
import edu.duke.cs.osprey.gui.io.withService
import io.kotlintest.shouldBe


class TestMissingAtoms : SharedSpec({

	fun List<MissingAtom>.shouldHave(expAtom: Atom, expRes: Polymer.Residue?) {
		any { obs ->
			obs.atom.name == expAtom.name && obs.res?.id == expRes?.id
		} shouldBe true
	}

	test("1cc8") {
		withService {

			val mol = Molecule.fromPDB(OspreyGui.getResourceAsString("1cc8.protein.pdb"))

			// 1cc8 should have all its heavy atoms already
			mol.inferMissingAtomsAmber().size shouldBe 0

			// but let's delete a couple atoms
			val A = (mol as Polymer).chains.find { it.id == "A" }!!
			val A23 = A.residues.find { it.id == "23" }!!
			val A23N = A23.atoms.find { it.name == "N" }!!
			mol.atoms.remove(A23N)
			A23.atoms.remove(A23N)
			val A45 = A.residues.find { it.id == "45" }!!
			val A45CA = A45.atoms.find { it.name == "CA" }!!
			mol.atoms.remove(A45CA)
			A45.atoms.remove(A45CA)

			// they should show up again in the missing atoms list
			val missingAtoms = mol.inferMissingAtomsAmber()
			missingAtoms.shouldHave(A23N, A23)
			missingAtoms.shouldHave(A45CA, A45)
			missingAtoms.size shouldBe 2
		}
	}
})
