package edu.duke.cs.osprey.molscope.molecule

import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe


class TestMoleculeHashEquals : FunSpec({

	val mol1 = Molecule("test1").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.0, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1Same = Molecule("test1").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.0, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1Copy = mol1.copy()

	val mol1DifName = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.0, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1DifA1Name = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C3", 1.0, 2.0, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1DifA1Element = Molecule("test1").apply {
		val a1 = Atom(Element.Oxygen, "C1", 1.0, 2.0, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1DifA1x = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.5, 2.0, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1DifA1y = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.5, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1DifA1z = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.0, 3.5).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
	}

	val mol1MissingBonds = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.0, 3.5).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
	}

	val mol1AddedBonds = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.0, 3.5).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a1, a3)
		bonds.add(a2, a3)
	}

	val mol1DifBonds = Molecule("test1 different").apply {
		val a1 = Atom(Element.Carbon, "C1", 1.0, 2.0, 3.5).also { atoms.add(it) }
		val a2 = Atom(Element.Carbon, "C2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		val a3 = Atom(Element.Nitrogen, "N", 5.0, 6.0, 7.0).also { atoms.add(it) }
		bonds.add(a1, a2)
		bonds.add(a2, a3)
	}

	val mol2 = Molecule("test2").apply {
		val a1 = Atom(Element.Oxygen, "O1", 1.0, 2.0, 3.0).also { atoms.add(it) }
		val a2 = Atom(Element.Oxygen, "O2", 2.0, 3.0, 4.0).also { atoms.add(it) }
		bonds.add(a1, a2)
	}

	test("hashCode") {

		// hash codes should be deterministic
		mol1.hashCode() shouldBe mol1.hashCode()
		mol1Same.hashCode() shouldBe mol1Same.hashCode()
		mol2.hashCode() shouldBe mol2.hashCode()

		// identical molecules should have the same hash code
		mol1.hashCode() shouldBe mol1Same.hashCode()
		mol1.hashCode() shouldBe mol1Copy.hashCode()

		// different molecules should (probably) have different hash codes
		mol1.hashCode() shouldNotBe mol1DifName.hashCode()
		mol1.hashCode() shouldNotBe mol1DifA1Name.hashCode()
		mol1.hashCode() shouldNotBe mol1DifA1Element.hashCode()
		mol1.hashCode() shouldNotBe mol1DifA1x.hashCode()
		mol1.hashCode() shouldNotBe mol1DifA1y.hashCode()
		mol1.hashCode() shouldNotBe mol1DifA1z.hashCode()
		mol1.hashCode() shouldNotBe mol1MissingBonds.hashCode()
		mol1.hashCode() shouldNotBe mol1AddedBonds.hashCode()
		mol1.hashCode() shouldNotBe mol1DifBonds.hashCode()
	}

	test("equals") {

		mol1 shouldBe mol1
		mol1 shouldBe mol1Same
		mol1 shouldBe mol1Copy
		mol2 shouldBe mol2

		mol1 shouldNotBe mol1DifName
		mol1 shouldNotBe mol1DifA1Name
		mol1 shouldNotBe mol1DifA1Element
		mol1 shouldNotBe mol1DifA1x
		mol1 shouldNotBe mol1DifA1y
		mol1 shouldNotBe mol1DifA1z
		mol1 shouldNotBe mol1MissingBonds
		mol1 shouldNotBe mol1AddedBonds
		mol1 shouldNotBe mol1DifBonds

		mol1 shouldNotBe mol2
	}
})
