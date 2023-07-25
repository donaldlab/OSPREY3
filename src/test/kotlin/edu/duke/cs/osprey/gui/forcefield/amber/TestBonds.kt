package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.fromMol2
import edu.duke.cs.osprey.gui.io.fromPDB
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.shouldBe
import io.kotest.matchers.types.shouldBeSameInstanceAs


/**
 * Test that the bond generation codes match expected bonds for proteins and small molecules.
 *  Under the hood OSPREY uses LEaP or Antechamber to generate a MOL2 file, which has a list of bonds in it.
 **/
class TestBonds : FunSpec({

	fun Molecule.Bonds.toContentSet() = toSet().map { it.toContent() }.toSet()

	context("1cc8") {

		test("protein") {
			withService {

				val mol = Molecule.fromPDB(OspreyGui.getResourceAsString("1cc8.protein.pdb"))

				for ((a1, a2) in mol.inferBondsAmber()) {
					mol.bonds.add(a1, a2)
				}

				val bondedMol = Molecule.fromMol2(OspreyGui.getResourceAsString("1cc8.protein.sybyl.mol2"))
				mol.bonds.toContentSet() shouldBe bondedMol.bonds.toContentSet()
			}
		}

		test("benzamidine") {
			withService {

				val mol = Molecule.fromPDB(OspreyGui.getResourceAsString("benzamidine.pdb"))

				for ((a1, a2) in mol.inferBondsAmber()) {
					mol.bonds.add(a1, a2)
				}

				val bondedMol = Molecule.fromMol2(OspreyGui.getResourceAsString("benzamidine.sybyl.mol2"))
				mol.bonds.toContentSet() shouldBe bondedMol.bonds.toContentSet()
			}
		}

		test("no cross-partition bonds") {
			withService {

				val mol = Molecule.fromPDB(OspreyGui.getResourceAsString("1cc8.pdb"))

				val partition = mol.partition(combineSolvent = true)

				fun Atom.findMol() =
					partition
						.find { (_, mol) -> this in mol.atoms }

				for ((a, b) in mol.inferBondsAmber()) {
					val mola = a.findMol()!!
					val molb = b.findMol()!!
					mola shouldBeSameInstanceAs molb
				}
			}
		}
	}
})
