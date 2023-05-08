package edu.duke.cs.osprey.molscope.molecule

import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldContainExactly
import io.kotest.matchers.collections.shouldContainExactlyInAnyOrder
import io.kotest.matchers.shouldBe
import io.kotest.matchers.types.shouldBeSameInstanceAs
import io.kotest.matchers.types.shouldNotBeSameInstanceAs


class TestMoleculeCombine : FunSpec({
	
	fun AtomMap.assert(a: Atom, b: Atom) {
		getB(a) shouldBeSameInstanceAs b
		getA(b) shouldBeSameInstanceAs a
	}

	test("one small") {

		val mol = Molecule("small")
		val srcC = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol.atoms.add(it) }
		val srcN = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol.atoms.add(it) }
		mol.bonds.add(srcC, srcN)

		val (combined, atomMap) = listOf(mol).combine("combined")

		combined.atoms.size shouldBe 2
		val dstC = combined.atoms.find { it.name == "C" }!!
		val dstN = combined.atoms.find { it.name == "N" }!!
		srcC shouldNotBeSameInstanceAs dstC
		srcN shouldNotBeSameInstanceAs dstN

		combined.bonds.count() shouldBe 1
		combined.bonds.isBonded(dstC, dstN) shouldBe true

		atomMap.assert(srcC, dstC)
		atomMap.assert(srcN, dstN)
	}

	test("two small") {

		val mol1 = Molecule("small1")
		val srcC1 = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol1.atoms.add(it) }
		val srcN1 = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol1.atoms.add(it) }
		mol1.bonds.add(srcC1, srcN1)

		val mol2 = Molecule("small2")
		val srcC2 = Atom(Element.Carbon, "C", 11.0, 12.0, 13.0).also { mol2.atoms.add(it) }
		val srcN2 = Atom(Element.Nitrogen, "N", 14.0, 15.0, 16.0).also { mol2.atoms.add(it) }
		mol2.bonds.add(srcC2, srcN2)

		val (combined, atomMap) = listOf(mol1, mol2).combine("combined")

		combined.atoms.size shouldBe 4
		val dstC1 = combined.atoms.find { it.name == "C" && it.pos.x == 1.0 }!!
		val dstN1 = combined.atoms.find { it.name == "N" && it.pos.x == 4.0 }!!
		val dstC2 = combined.atoms.find { it.name == "C" && it.pos.x == 11.0 }!!
		val dstN2 = combined.atoms.find { it.name == "N" && it.pos.x == 14.0 }!!
		srcC1 shouldNotBeSameInstanceAs dstC1
		srcN1 shouldNotBeSameInstanceAs dstN1
		srcC2 shouldNotBeSameInstanceAs dstC2
		srcN2 shouldNotBeSameInstanceAs dstN2

		combined.bonds.count() shouldBe 2
		combined.bonds.isBonded(dstC1, dstN1) shouldBe true
		combined.bonds.isBonded(dstC2, dstN2) shouldBe true

		atomMap.assert(srcC1, dstC1)
		atomMap.assert(srcN1, dstN1)
		atomMap.assert(srcC2, dstC2)
		atomMap.assert(srcN2, dstN2)
	}

	test("one polymer") {

		val mol = Polymer("poly")
		val srcC = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol.atoms.add(it) }
		val srcN = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol.atoms.add(it) }
		mol.bonds.add(srcC, srcN)
		val srcChain = Polymer.Chain("A").also { mol.chains.add(it) }
		val srcRes = Polymer.Residue("1", "RES", listOf(srcC, srcN)).also { srcChain.residues.add(it) }

		val (combined, atomMap) = listOf(mol).combine("combined")

		combined as Polymer

		combined.atoms.size shouldBe 2
		val dstC = combined.atoms.find { it.name == "C" }!!
		val dstN = combined.atoms.find { it.name == "N" }!!
		srcC shouldNotBeSameInstanceAs dstC
		srcN shouldNotBeSameInstanceAs dstN

		combined.bonds.count() shouldBe 1
		combined.bonds.isBonded(dstC, dstN) shouldBe true

		combined.chains.size shouldBe 1
		val dstChain = combined.chains.find { it.id == "A" }!!
		val dstRes = dstChain.residues.find { it.id == "1" }!!
		srcChain shouldNotBeSameInstanceAs dstChain
		srcRes shouldNotBeSameInstanceAs dstRes
		dstRes.atoms shouldContainExactlyInAnyOrder  listOf(dstC, dstN)

		atomMap.assert(srcC, dstC)
		atomMap.assert(srcN, dstN)
	}

	test("two polymer") {

		val mol1 = Polymer("poly1")
		val srcC1 = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol1.atoms.add(it) }
		val srcN1 = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol1.atoms.add(it) }
		mol1.bonds.add(srcC1, srcN1)
		val srcChain1 = Polymer.Chain("A").also { mol1.chains.add(it) }
		val srcRes1 = Polymer.Residue("1", "RES", listOf(srcC1, srcN1)).also { srcChain1.residues.add(it) }

		val mol2 = Polymer("poly1")
		val srcC2 = Atom(Element.Carbon, "C", 11.0, 12.0, 13.0).also { mol2.atoms.add(it) }
		val srcN2 = Atom(Element.Nitrogen, "N", 14.0, 15.0, 16.0).also { mol2.atoms.add(it) }
		mol2.bonds.add(srcC2, srcN2)
		val srcChain2 = Polymer.Chain("B").also { mol2.chains.add(it) }
		val srcRes2 = Polymer.Residue("1", "RES", listOf(srcC2, srcN2)).also { srcChain2.residues.add(it) }

		val (combined, atomMap) = listOf(mol1, mol2).combine("combined")

		combined as Polymer

		combined.atoms.size shouldBe 4
		val dstC1 = combined.atoms.find { it.name == "C" && it.pos.x == 1.0 }!!
		val dstN1 = combined.atoms.find { it.name == "N" && it.pos.x == 4.0 }!!
		val dstC2 = combined.atoms.find { it.name == "C" && it.pos.x == 11.0 }!!
		val dstN2 = combined.atoms.find { it.name == "N" && it.pos.x == 14.0 }!!
		srcC1 shouldNotBeSameInstanceAs dstC1
		srcN1 shouldNotBeSameInstanceAs dstN1
		srcC2 shouldNotBeSameInstanceAs dstC2
		srcN2 shouldNotBeSameInstanceAs dstN2

		combined.bonds.count() shouldBe 2
		combined.bonds.isBonded(dstC1, dstN1) shouldBe true
		combined.bonds.isBonded(dstC2, dstN2) shouldBe true

		combined.chains.size shouldBe 2
		val dstChain1 = combined.chains.find { it.id == "A" }!!
		val dstRes1 = dstChain1.residues.find { it.id == "1" }!!
		val dstChain2 = combined.chains.find { it.id == "B" }!!
		val dstRes2 = dstChain2.residues.find { it.id == "1" }!!
		srcChain1 shouldNotBeSameInstanceAs dstChain1
		srcRes1 shouldNotBeSameInstanceAs dstRes1
		dstRes1.atoms shouldContainExactlyInAnyOrder listOf(dstC1, dstN1)
		srcChain2 shouldNotBeSameInstanceAs dstChain2
		srcRes2 shouldNotBeSameInstanceAs dstRes2
		dstRes2.atoms shouldContainExactlyInAnyOrder listOf(dstC2, dstN2)

		atomMap.assert(srcC1, dstC1)
		atomMap.assert(srcN1, dstN1)
		atomMap.assert(srcC2, dstC2)
		atomMap.assert(srcN2, dstN2)
	}

	test("two polymer, chain id collision") {

		val mol1 = Polymer("poly1")
		val srcC1 = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol1.atoms.add(it) }
		val srcN1 = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol1.atoms.add(it) }
		mol1.bonds.add(srcC1, srcN1)
		val srcChain1 = Polymer.Chain("A").also { mol1.chains.add(it) }
		Polymer.Residue("1", "RES", listOf(srcC1, srcN1)).also { srcChain1.residues.add(it) }

		val mol2 = Polymer("poly1")
		val srcC2 = Atom(Element.Carbon, "C", 11.0, 12.0, 13.0).also { mol2.atoms.add(it) }
		val srcN2 = Atom(Element.Nitrogen, "N", 14.0, 15.0, 16.0).also { mol2.atoms.add(it) }
		mol2.bonds.add(srcC2, srcN2)
		val srcChain2 = Polymer.Chain("A").also { mol2.chains.add(it) }
		Polymer.Residue("1", "RES", listOf(srcC2, srcN2)).also { srcChain2.residues.add(it) }

		shouldThrow<IllegalArgumentException> {
			listOf(mol1, mol2).combine("combined")
		}
	}

	test("two polymer, chain id resolve") {

		val mol1 = Polymer("poly1")
		val srcC1 = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol1.atoms.add(it) }
		val srcN1 = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol1.atoms.add(it) }
		mol1.bonds.add(srcC1, srcN1)
		val srcChain1 = Polymer.Chain("A").also { mol1.chains.add(it) }
		Polymer.Residue("1", "RES", listOf(srcC1, srcN1)).also { srcChain1.residues.add(it) }

		val mol2 = Polymer("poly1")
		val srcC2 = Atom(Element.Carbon, "C", 11.0, 12.0, 13.0).also { mol2.atoms.add(it) }
		val srcN2 = Atom(Element.Nitrogen, "N", 14.0, 15.0, 16.0).also { mol2.atoms.add(it) }
		mol2.bonds.add(srcC2, srcN2)
		val srcChain2 = Polymer.Chain("A").also { mol2.chains.add(it) }
		Polymer.Residue("1", "RES", listOf(srcC2, srcN2)).also { srcChain2.residues.add(it) }

		val generator = object : ChainIdGenerator {
			override fun setUsedIds(ids: Collection<String>) {}
			override fun generateId() = "B"
		}
		val (combined, _) = listOf(mol1, mol2).combine("combined", generator)

		(combined as Polymer).chains.map { it.id } shouldContainExactly listOf("A", "B")
	}

	test("one small, one polymer") {

		val mol1 = Molecule("small1")
		val srcC1 = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol1.atoms.add(it) }
		val srcN1 = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol1.atoms.add(it) }
		mol1.bonds.add(srcC1, srcN1)

		val mol2 = Polymer("poly1")
		val srcC2 = Atom(Element.Carbon, "C", 11.0, 12.0, 13.0).also { mol2.atoms.add(it) }
		val srcN2 = Atom(Element.Nitrogen, "N", 14.0, 15.0, 16.0).also { mol2.atoms.add(it) }
		mol2.bonds.add(srcC2, srcN2)
		val srcChain2 = Polymer.Chain("A").also { mol2.chains.add(it) }
		Polymer.Residue("1", "RES", listOf(srcC2, srcN2)).also { srcChain2.residues.add(it) }

		val (combined, _) = listOf(mol1, mol2).combine("combined")

		combined as Polymer
		combined.chains.map { it.id } shouldBe listOf("A")
		combined.chains.find { it.id == "A" }!!.apply {
            residues.size shouldBe 1
            residues.sumOf { it.atoms.size } shouldBe mol2.atoms.size
        }
		combined.atoms.size shouldBe mol1.atoms.size + mol2.atoms.size
	}

	test("one small, one polymer, generate chain") {

		val mol1 = Molecule("small1")
		val srcC1 = Atom(Element.Carbon, "C", 1.0, 2.0, 3.0).also { mol1.atoms.add(it) }
		val srcN1 = Atom(Element.Nitrogen, "N", 4.0, 5.0, 6.0).also { mol1.atoms.add(it) }
		mol1.bonds.add(srcC1, srcN1)

		val mol2 = Polymer("poly1")
		val srcC2 = Atom(Element.Carbon, "C", 11.0, 12.0, 13.0).also { mol2.atoms.add(it) }
		val srcN2 = Atom(Element.Nitrogen, "N", 14.0, 15.0, 16.0).also { mol2.atoms.add(it) }
		mol2.bonds.add(srcC2, srcN2)
		val srcChain2 = Polymer.Chain("A").also { mol2.chains.add(it) }
		Polymer.Residue("1", "RES", listOf(srcC2, srcN2)).also { srcChain2.residues.add(it) }

		val generator = object : ChainGenerator {
			override fun setUsedIds(ids: Collection<String>) {}
			override fun generateChain(nonPolymerMol: Molecule, polymerMol: Polymer, polymerAtoms: List<Atom>) =
				Polymer.Chain("B").apply {
					residues.add(Polymer.Residue(
						"1",
						"MOL",
						polymerAtoms
					))
				}
		}
		val (combined, _) = listOf(mol1, mol2).combine("combined", chainGenerator=generator)

		combined as Polymer
		combined.chains.map { it.id } shouldBe listOf("A", "B")
		combined.chains.find { it.id == "A" }!!.apply {
            residues.size shouldBe 1
            residues.sumOf { it.atoms.size } shouldBe mol2.atoms.size
        }
		combined.chains.find { it.id == "B" }!!.apply {
			residues.size shouldBe 1
			residues.first().apply {
				atoms.size shouldBe mol1.atoms.size
			}
		}
		combined.atoms.size shouldBe mol1.atoms.size + mol2.atoms.size
	}
})
