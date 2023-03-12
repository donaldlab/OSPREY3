package edu.duke.cs.osprey.molscope.molecule

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.fromPDB
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldBeSameSizeAs
import io.kotest.matchers.ints.shouldBeExactly
import io.kotest.matchers.shouldBe
import io.kotest.matchers.types.shouldNotBeSameInstanceAs
import org.joml.Vector3d

class TestMoleculeTransform : FunSpec({

    fun flipOverZ(pos: Vector3d): Vector3d = Vector3d(pos.x, pos.y, -pos.z)

    fun List<Atom>.shouldBeFlippedOverZ(other: List<Atom>) {
        val src = this
        src.zip(other).forEach { (first, second) ->
            first.apply {
                name shouldBe second.name
                element shouldBe second.element
                pos.apply {
                    x shouldBe second.pos.x
                    y shouldBe second.pos.y
                    (z + second.pos.z) shouldBe 0.0
                }
            }
        }
    }

    fun List<Atom>.shouldBeSame(other: List<Atom>) {
        val src = this
        src.zip(other).forEach {(first, second) ->
            first.apply {
                name shouldBe second.name
                element shouldBe second.element
                pos.apply {
                    x shouldBe second.pos.x
                    y shouldBe second.pos.y
                    z shouldBe second.pos.z
                }
            }
        }
    }

    test("inverting chirality of polymer updates positions of all atoms") {

        val mol = Molecule.fromPDB(OspreyGui.getResourceAsString("1cc8.protein.pdb")) as Polymer
        val inverted = mol.transformCopy{ flipOverZ(it) }.first

        // basic checks on things that should be the same
        inverted shouldNotBeSameInstanceAs mol
        inverted.name shouldBe mol.name
        inverted.netCharge shouldBe mol.netCharge
        inverted.type shouldBe mol.type
        inverted.atoms shouldBeSameSizeAs  mol.atoms
        inverted.chains shouldBeSameSizeAs mol.chains
        inverted.chains.flatMap { it.residues } shouldBeSameSizeAs mol.chains.flatMap { it.residues }
        inverted.bonds.count() shouldBeExactly mol.bonds.count()

        // Check that the atom positions are different
        inverted.atoms.shouldBeFlippedOverZ(mol.atoms)

        // Check that round-tripping returns atoms to their original positions
        val roundTrip = inverted.transformCopy { flipOverZ(it) }.first
        roundTrip.atoms.shouldBeSame(mol.atoms)
    }
})