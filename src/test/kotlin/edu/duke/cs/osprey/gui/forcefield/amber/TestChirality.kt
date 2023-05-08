package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.fromPDB
import edu.duke.cs.osprey.gui.io.withService
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.booleans.shouldBeTrue
import io.kotest.matchers.collections.shouldContainExactly
import io.kotest.matchers.collections.shouldHaveSize
import io.kotest.matchers.comparables.shouldBeEqualComparingTo
import io.kotest.matchers.equality.shouldBeEqualToComparingFields
import io.kotest.matchers.equality.shouldNotBeEqualToComparingFields
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import io.kotest.matchers.types.shouldBeSameInstanceAs
import io.kotest.matchers.types.shouldNotBeSameInstanceAs
import org.joml.Vector3d

/* Tests of the preparation of molecules that are dextrorotatory. */
class TestChirality : FunSpec({

    fun invertOverZAxis(pos: Vector3d) = Vector3d(pos.x, pos.y, -pos.z)

    context("in-place versus copying inversions") {
        val lThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/thanatin-chain-b.pdb"))
        val dThanatin = lThanatin.invertedCopy()

        lThanatin shouldNotBeSameInstanceAs dThanatin

        dThanatin.atoms.zip(lThanatin.atoms).forEach { (cpyAtom, originalAtom) ->
            cpyAtom.pos shouldNotBeEqualToComparingFields originalAtom.pos
        }

        // flip it once and back
        dThanatin.atoms[0].invertInPlace()
        dThanatin.atoms[0].pos.z shouldBeEqualComparingTo lThanatin.atoms[0].pos.z
        dThanatin.atoms[0].invertInPlace()

        val inPlaceInversion = lThanatin.invertedInPlace()
        lThanatin shouldBeSameInstanceAs inPlaceInversion

        dThanatin.atoms.zip(inPlaceInversion.atoms).forEach { (cpyAtom, inPlaceAtom) ->
            cpyAtom.pos shouldBeEqualToComparingFields inPlaceAtom.pos
        }
    }

    context("missing atoms") {
        val lThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/thanatin-chain-b.pdb"))
        /* this structure is the above structure flipped over the z-axis */
        val dThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/d-thanatin-chain-b.pdb"))

        // The structure has all the heavy atoms. Remove some so that we can use LEaP to add them back
        val atx1ProteinMissingAtoms = Molecule.fromPDB(OspreyGui.getResourceAsString("1cc8.protein.pdb")).apply {
            val A = (this as Polymer).chains.find { it.id == "A" }!!
            val A23 = A.residues.find { it.id == "23" }!!
            val A23N = A23.atoms.find { it.name == "N" }!!
            this.atoms.remove(A23N)
            A23.atoms.remove(A23N)
            val A45 = A.residues.find { it.id == "45" }!!
            val A45CA = A45.atoms.find { it.name == "CA" }!!
            this.atoms.remove(A45CA)
            A45.atoms.remove(A45CA)
        }

        test("LEaP places D missing atoms not in mirror image of where it places L missing atoms") {
            // This demonstrates that to use the missing atoms function we need to flip into L space first

            withService {
                val invertedThanatin = lThanatin.transformCopy { invertOverZAxis(it) }.first
                val missingLAtoms = lThanatin.inferMissingAtomsAmber()
                val missingDAtoms = invertedThanatin.inferMissingAtomsAmber()

                // the list of missing atoms should be the same
                missingDAtoms.map { it.toString() } shouldContainExactly missingLAtoms.map { it.toString() }

                missingDAtoms.zip(missingLAtoms).forEach { (missingD, missingL) ->
                    missingD.atom.pos.x shouldNotBe missingL.atom.pos.x
                    missingD.atom.pos.y shouldNotBe missingL.atom.pos.y
                    missingD.atom.pos.z shouldNotBe missingL.atom.pos.z
                }

                // show this problem with 1cc8, too.
                val missingAtomsInAtx1 = atx1ProteinMissingAtoms.inferMissingAtomsAmber()
                missingAtomsInAtx1 shouldHaveSize 2
                val invertedAtx1 = atx1ProteinMissingAtoms.transformCopy { invertOverZAxis(it) }.first
                val invertedAtx1MissingAtoms = invertedAtx1.inferMissingAtomsAmber()

                invertedAtx1MissingAtoms.zip(missingAtomsInAtx1).forEach { (missingD, missingL) ->
                    missingD.atom.pos.x shouldNotBe missingL.atom.pos.x
                    missingD.atom.pos.y shouldNotBe missingL.atom.pos.y
                    missingD.atom.pos.z shouldNotBe missingL.atom.pos.z
                }
            }
        }

        test("Mirroring D-polymers into L space allows LEaP to correctly identify positions of missing atoms") {
            withService {
                val invertedDThanatin = dThanatin.transformCopy { invertOverZAxis(it) }.first // this is L-thanatin
                val missingInvertedDAtoms = invertedDThanatin.inferMissingAtomsAmber()
                val missingLAtoms = lThanatin.inferMissingAtomsAmber()

                missingInvertedDAtoms.zip(missingLAtoms).forEach{(invertedDAtom, lAtom) ->
                    invertedDAtom.atom.pos shouldBeEqualToComparingFields lAtom.atom.pos
                }
            }
        }
    }

    fun Molecule.Bonds.toContentSet() = toSet().map { it.toContent() }.toSet()

    context("bonds") {

        test("bonds can be inferred correctly by flipping D to L chirality") {
            val lThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/thanatin-chain-b.pdb"))
            /* this structure is the above structure flipped over the z-axis */
            val dThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/d-thanatin-chain-b.pdb"))

            withService {

                for (bond in lThanatin.inferBondsAmber()) {
                    lThanatin.bonds.add(bond)
                }

                val dInverted = dThanatin.transformCopy{ invertOverZAxis(it) }.first
                for (bond in dInverted.inferBondsAmber()) {
                    dInverted.bonds.add(bond)
                }

                lThanatin.bonds.toContentSet() shouldBe dInverted.bonds.toContentSet()
            }
        }

        test("bonds are also inferred correctly in D-polymers") {
            val lThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/thanatin-chain-b.pdb"))
            /* this structure is the above structure flipped over the z-axis */
            val dThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/d-thanatin-chain-b.pdb"))

            withService {
                val lBonds = lThanatin.inferBondsAmber()
                val dBonds = dThanatin.inferBondsAmber()

                dBonds shouldHaveSize lBonds.size
                dBonds.map { it.toString() } shouldContainExactly lBonds.map { it.toString() }
            }
        }
    }

    context("protonation") {

        fun List<ProtonatedAtom>.isSameAs(reference: List<ProtonatedAtom>) {
            for ((testAtom, referenceAtom) in zip(reference)) {
                testAtom.heavy shouldBeEqualToComparingFields referenceAtom.heavy
                testAtom.light shouldBeEqualToComparingFields referenceAtom.light
            }
        }

        fun List<ProtonatedAtom>.allAtomsAreMirroredFrom(reference: List<ProtonatedAtom>) {

            this.size shouldBe reference.size

            for ((testAtom, referenceAtom) in zip(reference)) {
                testAtom.location shouldBe referenceAtom.location // refers to location in chain/residue
                testAtom.heavy shouldBeEqualToComparingFields referenceAtom.heavy
                testAtom.light shouldBeEqualToComparingFields referenceAtom.light
            }
        }

        fun List<ProtonatedAtom>.atomPositionsAreNotMirroredFrom(reference: List<ProtonatedAtom>) {
            this.size shouldBe reference.size

            zip(reference)
                .any { (testAtom, referenceAtom) ->
                    val mirroredLight = testAtom.light.copy { invertOverZAxis(it) }
                    mirroredLight.pos != referenceAtom.light.pos
                }
                .shouldBeTrue()
        }

        test("protonating a molecule the second time returns the same results as the first time") {
            val lThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/thanatin-chain-b.pdb"))

            withService {
                val firstGo = lThanatin.inferProtonation()
                val secondGo = lThanatin.inferProtonation()

                secondGo.size shouldBe firstGo.size
                firstGo.isSameAs(secondGo)
            }
        }

        /* demonstrate that we need to flip into L space to correctly protonate D polymer */
        test("LEaP does not protonate D-proteins in the mirror image of their L-counterparts") {
            val lThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/thanatin-chain-b.pdb"))
            val dThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/d-thanatin-chain-b.pdb"))
            lThanatin.deprotonate()
            dThanatin.deprotonate()

            withService {
                val lProtonation = lThanatin.inferProtonation()
                val dProtonation = dThanatin.inferProtonation()

                dProtonation.atomPositionsAreNotMirroredFrom(lProtonation) // unfortunately...
            }
        }

        test("Flipping D-polymer into L-space, protonating, then flipping back correctly adds protons") {
            val lThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/thanatin-chain-b.pdb"))
            val dThanatin = Molecule.fromPDB(OspreyGui.getResourceAsString("dexDesign/d-thanatin-chain-b.pdb"))
            lThanatin.deprotonate()
            dThanatin.deprotonate()

            val invertedD = dThanatin.transformCopy { invertOverZAxis(it) }.first
            invertedD.inferProtonation().allAtomsAreMirroredFrom(lThanatin.inferProtonation())
        }
    }

    context("net charges") {
        // empty
    }

    context("minimization") {
        // empty
    }
})