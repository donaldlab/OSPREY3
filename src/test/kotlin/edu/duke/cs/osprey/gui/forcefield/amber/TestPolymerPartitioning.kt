package edu.duke.cs.osprey.gui.forcefield.amber


import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.fromPDB
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.molscope.molecule.combine
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.maps.shouldNotContainKey
import io.kotest.matchers.shouldBe
import io.kotest.matchers.types.shouldBeTypeOf
import org.joml.Vector3d


class TestPolymerPartitioning : FunSpec({

	context("1CC8") {

		val mol = Molecule.fromPDB(OspreyGui.getResourceAsString("1cc8.pdb")) as Polymer

		test("partition") {

			val partition = mol
				.partition(combineSolvent = false)
				.groupBy(
					keySelector = { (type, _) -> type },
					valueTransform = { (_, mol) -> mol }
				)

			// check the molecule type counts
			partition[MoleculeType.Protein]!!.size shouldBe 1
			partition[MoleculeType.SmallMolecule]!!.size shouldBe 3
			partition[MoleculeType.Solvent]!!.size shouldBe 117
			partition shouldNotContainKey MoleculeType.DNA
			partition shouldNotContainKey MoleculeType.RNA
			partition shouldNotContainKey MoleculeType.AtomicIon
			partition shouldNotContainKey MoleculeType.Synthetic

			// check the protein
			partition[MoleculeType.Protein]!!.first().apply {

				this as Polymer

				chains.size shouldBe 1
				chains[0].apply {
					id shouldBe "A"
					residues.size shouldBe 72

					// TODO: check a few resides?
				}

				atoms.size shouldBe 567
			}

			partition[MoleculeType.SmallMolecule]!!.apply {

				// check the mercury atom
				find { it.type == "HG" }!!.apply {
					atoms.size shouldBe 1
					atoms[0].apply {
						name shouldBe "HG"
						element shouldBe Element.Mercury
					}
				}

				// check one of the benzamidines
				find { it.type == "BEN" }!!.apply {

					atoms.size shouldBe 9

					// spot check a couple atoms
					atoms.find { it.name == "C1" }!!.apply {
						element shouldBe Element.Carbon
						pos shouldBe Vector3d(6.778, 10.510, 20.665)
					}
					atoms.find { it.name == "N2" }!!.apply {
						element shouldBe Element.Nitrogen
						pos shouldBe Vector3d(4.965, 11.590, 21.821)
					}
				}
			}
		}

		test("combine solvent") {

			val partition = mol.partition(combineSolvent = true)
			val solvents = partition
				.filter { (moltype, _) -> moltype == MoleculeType.Solvent }

			solvents.size shouldBe 1

			val solvent = solvents.first().second

			// all the solvent molecules should be in a single "molecule"
			solvent.shouldBeTypeOf<Molecule>();
			solvent.atoms.size shouldBe 117
		}

		test("combine protein, Hg, benzamidines") {

			mol.partition()
				.filter { (type, _) -> type in setOf(MoleculeType.Protein, MoleculeType.SmallMolecule) }
				.map { (_, mol) -> mol }
				.combine("1CC8").first
				.apply {

					this as Polymer

					// the protein gets its own chain,
					// and all the small molecules don't get any chains
					chains.size shouldBe 1

					// check the protein
					chains.find { it.id == "A" }!!.apply {
						residues.size shouldBe 72
					}

					// check the mercury atom
					atoms.find { it.name == "HG" }!!.apply {
						element shouldBe Element.Mercury
					}

					// spot check a couple atoms in one of the benzamidines
					atoms.find { it.name == "C1" && it.pos.x == 6.778 }!!.apply {
						element shouldBe Element.Carbon
						pos shouldBe Vector3d(6.778, 10.510, 20.665)
					}
					atoms.find { it.name == "N2" && it.pos.x == 4.965 }!!.apply {
						element shouldBe Element.Nitrogen
						pos shouldBe Vector3d(4.965, 11.590, 21.821)
					}

					atoms.size shouldBe 567 + 1 + 9 + 9
				}
		}

		test("partition all, combine") {
			mol.partition(combineSolvent = false)
				.map { (_, mol) -> mol }
				.combine("combined").first
				.run {

					this as Polymer

					// should have just one chain, from the protein
					chains.size shouldBe 1
					chains.map { it.id }.toSet() shouldBe setOf("A")

					// but we should still have all of the atoms
					atoms.size shouldBe 567 + 117 + 1 + 9 + 9
				}
		}

		test("partition all, combined solvent, combine") {
			mol.partition(combineSolvent = true)
				.map { (_, mol) -> mol }
				.combine("combined").first
				.run {

					this as Polymer

					// should have just one chain, from the protein
					chains.size shouldBe 1
					chains.map { it.id }.toSet() shouldBe setOf("A")

					// but we should still have all of the atoms
					atoms.size shouldBe 567 + 117 + 1 + 9 + 9
				}
		}
	}
})
