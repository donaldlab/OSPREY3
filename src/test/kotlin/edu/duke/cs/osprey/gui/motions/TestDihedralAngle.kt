package edu.duke.cs.osprey.gui.motions

import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.SharedSpec
import edu.duke.cs.osprey.gui.io.fromOMOL
import io.kotlintest.matchers.doubles.plusOrMinus
import io.kotlintest.shouldBe


class TestDihedralAngle : SharedSpec({

	fun Double.shouldBeAbsolutely(expected: Double, epsilon: Double = 1e-3) {
		this shouldBe expected.plusOrMinus(epsilon)
	}

	val mol = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0] as Polymer

	test("measure") {

		val lys24 = mol.findChainOrThrow("A").findResidueOrThrow("24")

		DihedralAngle(
			mol,
			lys24.findAtomOrThrow("N"),
			lys24.findAtomOrThrow("CA"),
			lys24.findAtomOrThrow("CB"),
			lys24.findAtomOrThrow("CG")
		).measureDegrees().shouldBeAbsolutely(-175.368)

		DihedralAngle(
			mol,
			lys24.findAtomOrThrow("CA"),
			lys24.findAtomOrThrow("CB"),
			lys24.findAtomOrThrow("CG"),
			lys24.findAtomOrThrow("CD")
		).measureDegrees().shouldBeAbsolutely(-148.553)

		DihedralAngle(
			mol,
			lys24.findAtomOrThrow("CB"),
			lys24.findAtomOrThrow("CG"),
			lys24.findAtomOrThrow("CD"),
			lys24.findAtomOrThrow("CE")
		).measureDegrees().shouldBeAbsolutely(-62.106)

		DihedralAngle(
			mol,
			lys24.findAtomOrThrow("CG"),
			lys24.findAtomOrThrow("CD"),
			lys24.findAtomOrThrow("CE"),
			lys24.findAtomOrThrow("NZ")
		).measureDegrees().shouldBeAbsolutely(-85.934)
	}

	test("set") {

		val ser39 = mol.findChainOrThrow("A").findResidueOrThrow("39")

		fun DihedralAngle.check(degrees: Double) {
			setDegrees(degrees)
			measureDegrees().shouldBeAbsolutely(degrees)
		}

		DihedralAngle(
			mol,
			ser39.findAtomOrThrow("N"),
			ser39.findAtomOrThrow("CA"),
			ser39.findAtomOrThrow("CB"),
			ser39.findAtomOrThrow("OG")
		).run {

			rotatedAtoms.map { it.name }.toSet() shouldBe setOf("HB2", "HB3", "OG", "HG")
			check(-180.0)
			check(-150.0)
			check(-120.0)
			check(-90.0)
			check(-60.0)
			check(-30.0)
			check(0.0)
			check(30.0)
			check(60.0)
			check(90.0)
			check(120.0)
			check(150.0)
			check(180.0)
		}
	}
})
