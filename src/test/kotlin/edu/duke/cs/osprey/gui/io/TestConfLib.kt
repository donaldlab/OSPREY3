package edu.duke.cs.osprey.gui.io


import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.molscope.molecule.Element
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldContainExactly
import io.kotest.matchers.collections.shouldContainExactlyInAnyOrder
import io.kotest.matchers.collections.shouldExist
import io.kotest.matchers.maps.shouldContainExactly
import io.kotest.matchers.shouldBe
import io.kotest.matchers.shouldNotBe
import io.kotest.matchers.types.shouldBeTypeOf
import org.joml.Vector3d


class TestConfLib : FunSpec({

	test("read Lovell") {

		val conflib = ConfLib.from(OspreyGui.getResourceAsString("conflib/lovell.conflib"))

		// spot check a few bits
		conflib.name shouldBe "Amino Acids: The Penultimate Rotamer Library (with Hydroxyl rotamers)"
		conflib.description shouldNotBe null
		conflib.citation shouldNotBe null

		conflib.fragments.getValue("ALA").run {

			id shouldBe "ALA"
			name shouldBe "Alanine"
			type shouldBe "ALA"

			atoms.size shouldBe 6
			val ha = atoms.find { it.name == "HA" }!!.apply {
				element shouldBe Element.Hydrogen
			}
			val hb1 = atoms.find { it.name == "HB1" }!!.apply {
				element shouldBe Element.Hydrogen
			}
			val cb = atoms.find { it.name == "CB" }!!.apply {
				element shouldBe Element.Carbon
			}
			val h = atoms.find { it.name == "H" }!!.apply {
				element shouldBe Element.Hydrogen
			}

			bonds.size shouldBe 3
			bonds.shouldExist { it.a == hb1 && it.b == cb }

			anchors.size shouldBe 2
			val anchorCA = (anchors.find { it.id == 1 } as ConfLib.Anchor.Single).apply {
				bonds.shouldContainExactlyInAnyOrder(ha, cb)
			}
			val anchorN = (anchors.find { it.id == 2 } as ConfLib.Anchor.Single).apply {
				bonds.shouldContainExactlyInAnyOrder(h)
			}

			motions.size shouldBe 1
			motions[0].shouldBeTypeOf<ConfLib.ContinuousMotion.DihedralAngle> { motion ->
				motion.id shouldBe 1
				motion.a.shouldBeTypeOf<ConfLib.AnchorAtomPointer> { a ->
					a.anchor.id shouldBe 1
					a.index shouldBe 1
				}
				motion.b.shouldBeTypeOf<ConfLib.AnchorAtomPointer> { b ->
					b.anchor.id shouldBe 1
					b.index shouldBe 0
				}
				motion.c.shouldBeTypeOf<ConfLib.AtomInfo> { c ->
					c.name shouldBe "CB"
				}
				motion.d.shouldBeTypeOf<ConfLib.AtomInfo> { d ->
					d.name shouldBe "HB1"
				}
			}

			confs.size shouldBe 1
			confs.getValue("ALA").run {

				coords.size shouldBe 6
				coords[cb] shouldBe Vector3d(19.617, 5.198, 24.407)

				anchorCoords.size shouldBe 2
				anchorCoords.getValue(anchorCA).shouldBeTypeOf<ConfLib.AnchorCoords.Single> {
					it.b shouldBe Vector3d(20.267, 6.768, 22.647)
				}
				anchorCoords.getValue(anchorN).shouldBeTypeOf<ConfLib.AnchorCoords.Single> {
					it.c shouldBe Vector3d(10.668, 14.875, 20.325)
				}
			}
		}

		conflib.fragments.getValue("VAL").run {

			id shouldBe "VAL"
			name shouldBe "Valine"
			type shouldBe "VAL"

			atoms.size shouldBe 12
			bonds.size shouldBe 9
			anchors.size shouldBe 2
			motions.size shouldBe 3

			confs.size shouldBe 3
			confs.getValue("p").run {
				name shouldBe "p"
				coords.size shouldBe 12
				coords[atoms.find { it.name == "HA" }] shouldBe Vector3d(108.613, 17.002, -4.136)
				coords[atoms.find { it.name == "CG2" }] shouldBe Vector3d(106.309809, 15.583253, -4.049056)
			}
		}

		conflib.fragments.getValue("PRO").run {

			id shouldBe "PRO"
			name shouldBe "Proline"
			type shouldBe "PRO"

			atoms.size shouldBe 10
			bonds.size shouldBe 8

			anchors.size shouldBe 1
			val anchor = (anchors.find { it.id == 1 } as ConfLib.Anchor.Double).apply {
				bondsa.shouldContainExactlyInAnyOrder(
					atoms.find { it.name == "CB" },
					atoms.find { it.name == "HA" }
				)
				bondsb.shouldContainExactlyInAnyOrder(
					atoms.find { it.name == "CD" }
				)
			}

			motions.size shouldBe 0

			confs.size shouldBe 2
			confs.getValue("up").run {

				coords.size shouldBe 10
				coords[atoms.find { it.name == "HG3" }] shouldBe Vector3d(21.152899, 20.472575, 35.499174)

				anchorCoords.getValue(anchor).shouldBeTypeOf<ConfLib.AnchorCoords.Double> {
					it.d shouldBe Vector3d(21.102527, 24.472094, 37.053483)
				}
			}
		}

		// try an N-terminal residue
		conflib.fragments.getValue("VALn").run {

			id shouldBe "VALn"
			name shouldBe "Valine, N-terminal"
			type shouldBe "VAL"

			atoms.size shouldBe 11
			bonds.size shouldBe 9
			anchors.size shouldBe 1
			motions.size shouldBe 3

			confs.size shouldBe 3
			confs.getValue("p").run {
				name shouldBe "p"
				coords.size shouldBe 11
				coords[atoms.find { it.name == "HA" }] shouldBe Vector3d(108.613, 17.002, -4.136)
				coords[atoms.find { it.name == "CG2" }] shouldBe Vector3d(106.309809, 15.583253, -4.049056)
			}
		}
	}

	test("fragments list roundtrip") {

		// get a list of fragments
		val conflib = ConfLib.from(OspreyGui.getResourceAsString("conflib/lovell.conflib"))
		val expFrags = conflib.fragments
			.values
			.sortedBy { it.id }

		val (toml, idsByFrag) = expFrags.toToml(resolveIdCollisions = false)

		// check the ids
		for (frag in expFrags) {
			idsByFrag[frag] shouldBe frag.id
		}

		// do the roundtrip
		val obsFrags = ConfLib.fragmentsFrom(toml)
			.values
			.sortedBy { it.id }

		// see if we got the same frags
		obsFrags.size shouldBe expFrags.size
		for ((obsFrag, expFrag) in obsFrags.zip(expFrags)) {
			obsFrag shouldBeFrag expFrag
		}
	}

	fun emptyFrag(id: String) =
		ConfLib.Fragment(
			id = id,
			name = "Foo",
			type = "foo",
			atoms = emptyList(),
			bonds = emptyList(),
			anchors = emptyList(),
			confs = emptyMap(),
			motions = emptyList()
		)

	test("fragment name collision") {
		shouldThrow<IllegalArgumentException> {

			listOf(emptyFrag("foo"), emptyFrag("foo"))
				.toToml(resolveIdCollisions = false)
		}
	}

	test("fragment name collision resolved") {

		val frag1 = emptyFrag("foo")
		val frag2 = emptyFrag("foo")
		val frag3 = emptyFrag("foo")

		val (toml, idsByFrag) = listOf(frag1, frag2, frag3)
			.toToml(resolveIdCollisions = true)

		idsByFrag[frag1] shouldBe "foo"
		idsByFrag[frag2] shouldBe "foo2"
		idsByFrag[frag3] shouldBe "foo3"
	}
})

infix fun ConfLib.Fragment?.shouldBeFrag(exp: ConfLib.Fragment?) {

	// it's ok if they're both null
	if (this == null) {
		exp shouldBe null
	} else {

		// but otherwise they should both be not null
		exp shouldNotBe null
		val obsFrag = this
		val expFrag = exp!!

		obsFrag.id shouldBe expFrag.id
		obsFrag.name shouldBe expFrag.name
		obsFrag.type shouldBe expFrag.type

		obsFrag.atoms shouldContainExactly expFrag.atoms
		obsFrag.bonds shouldContainExactly expFrag.bonds
		obsFrag.anchors shouldContainExactly expFrag.anchors

		obsFrag.confs.keys shouldBe expFrag.confs.keys
		for ((obsConf, expConf) in obsFrag.confs.values.zip(expFrag.confs.values)) {
			obsConf shouldBeConf expConf
		}

		obsFrag.motions shouldContainExactly expFrag.motions
	}
}

infix fun ConfLib.Conf?.shouldBeConf(exp: ConfLib.Conf?) {

	// it's ok if they're both null
	if (this == null) {
		exp shouldBe null
	} else {

		// but otherwise they should both be not null
		exp shouldNotBe null
		val obsConf = this
		val expConf = exp!!

		obsConf.id shouldBe expConf.id
		obsConf.name shouldBe expConf.name
		obsConf.description shouldBe expConf.description

		// re-map these to replace the identity maps with regular maps
		// for value comparison instead of identity comparison
		obsConf.coords.toMap() shouldContainExactly expConf.coords.toMap()
		obsConf.anchorCoords.toMap() shouldContainExactly expConf.anchorCoords.toMap()
	}
}
