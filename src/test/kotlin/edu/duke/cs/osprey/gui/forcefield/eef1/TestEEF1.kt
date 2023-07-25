package edu.duke.cs.osprey.gui.forcefield.eef1

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.fromOMOL
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe


class TestEEF1 : FunSpec({

	context("1cc8") {

		val mol = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0] as Polymer

		fun Polymer.Residue.type(atomName: String) =
			findAtomOrThrow(atomName).atomTypeEEF1(mol)

		test("no missing heavy atoms") {

			mol.atoms
				.filter { it.element != Element.Hydrogen }
				.forEach {
					it.atomTypeEEF1OrThrow(mol)
				}
		}

		// spot-check a few amino acids, try to cover all the atom types

		test("arginine") {
			mol.findChainOrThrow("A").findResidueOrThrow("68").run {
				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("CG") shouldBe EEF1.AtomType.CH2E
				type("CD") shouldBe EEF1.AtomType.CH2E
				type("NE") shouldBe EEF1.AtomType.NH1
				type("CZ") shouldBe EEF1.AtomType.CR
				type("NH1") shouldBe EEF1.AtomType.NC2
				type("NH2") shouldBe EEF1.AtomType.NC2
			}
		}

		test("proline") {
			mol.findChainOrThrow("A").findResidueOrThrow("31").run {
				type("N") shouldBe EEF1.AtomType.N
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("CG") shouldBe EEF1.AtomType.CH2E
				type("CD") shouldBe EEF1.AtomType.CH2E
			}
		}

		test("histidine-dHD1") {
			mol.findChainOrThrow("A").findResidueOrThrow("6").run {

				// should be dHD1
				atoms.find { it.name == "HD1" } shouldBe null
				atoms.find { it.name == "HE2" } shouldNotBe null

				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("CG") shouldBe EEF1.AtomType.CR
				type("ND1") shouldBe EEF1.AtomType.NR
				type("CD2") shouldBe EEF1.AtomType.CR1E
				type("CE1") shouldBe EEF1.AtomType.CR1E
				type("NE2") shouldBe EEF1.AtomType.NH1
			}
		}

		test("glycine") {
			mol.findChainOrThrow("A").findResidueOrThrow("17").run {
				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH2E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
			}
		}

		test("alanine-Nterm") {
			mol.findChainOrThrow("A").findResidueOrThrow("2").run {
				type("N") shouldBe EEF1.AtomType.NH3
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH3E
			}
		}

		test("leucine-Cterm") {
			mol.findChainOrThrow("A").findResidueOrThrow("73").run {
				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.OC
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("CG") shouldBe EEF1.AtomType.CH1E
				type("CD1") shouldBe EEF1.AtomType.CH3E
				type("CD2") shouldBe EEF1.AtomType.CH3E
				type("OXT") shouldBe EEF1.AtomType.OC
			}
		}

		test("serine") {
			mol.findChainOrThrow("A").findResidueOrThrow("16").run {
				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("OG") shouldBe EEF1.AtomType.OH1
			}
		}

		test("cysteine") {
			mol.findChainOrThrow("A").findResidueOrThrow("18").run {
				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("SG") shouldBe EEF1.AtomType.SH1E
			}
		}

		test("methionine") {
			mol.findChainOrThrow("A").findResidueOrThrow("13").run {
				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("CG") shouldBe EEF1.AtomType.CH2E
				type("SD") shouldBe EEF1.AtomType.S
				type("CE") shouldBe EEF1.AtomType.CH3E
			}
		}

		test("asparagine") {
			mol.findChainOrThrow("A").findResidueOrThrow("10").run {
				type("N") shouldBe EEF1.AtomType.NH1
				type("CA") shouldBe EEF1.AtomType.CH1E
				type("C") shouldBe EEF1.AtomType.C
				type("O") shouldBe EEF1.AtomType.O
				type("CB") shouldBe EEF1.AtomType.CH2E
				type("CG") shouldBe EEF1.AtomType.C
				type("OD1") shouldBe EEF1.AtomType.O
				type("ND2") shouldBe EEF1.AtomType.NH2
			}
		}
	}
})
