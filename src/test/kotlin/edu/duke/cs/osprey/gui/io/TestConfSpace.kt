package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.motions.ConfMotion
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.motions.MolMotion
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.prep.DesignPosition
import edu.duke.cs.osprey.gui.prep.Proteins
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.molscope.tools.identityHashMapOf
import io.kotest.assertions.throwables.shouldThrow
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldNotBeEmpty
import io.kotest.matchers.shouldBe
import io.kotest.matchers.types.shouldNotBeSameInstanceAs
import org.joml.Vector3d


/** This conf space has a little bit of everything! */
fun makeTestConfSpace(): ConfSpace {

	// load some molecules
	val protein = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0] as Polymer
	val smallmol = Molecule.fromOMOL(OspreyGui.getResourceAsString("benzamidine.omol"))[0]

	// load some confs
	val conflib = ConfLib.from(OspreyGui.getResourceAsString("conflib/lovell.conflib"))

	// make the conf space
	return ConfSpace(listOf(
		MoleculeType.Protein to protein,
		MoleculeType.SmallMolecule to smallmol
	)).apply {

		name = "The Awesomest Conformation Space Evarrrr!!!"

		conflibs.add(conflib)

		// add some positions to the protein
		val leu26 = protein.findChainOrThrow("A").findResidueOrThrow("26")
		val thr27 = protein.findChainOrThrow("A").findResidueOrThrow("27")
		val pos1 = Proteins.makeDesignPosition(protein, leu26, "Pos1")
		val pos2 = Proteins.makeDesignPosition(protein, thr27, "Pos2")
		designPositionsByMol[protein] = mutableListOf(pos1, pos2)

		// set the pos conf spaces
		positionConfSpaces.getOrMake(pos1).apply {

			val leu = conflib.fragments.getValue("LEU")
			val ala = conflib.fragments.getValue("ALA")
			val pro = conflib.fragments.getValue("PRO")

			// add the wt frag
			val wtFrag = pos1.makeFragment("wt-${pos1.name.toTomlKey()}", "WildType", motions = leu.motions)
			wildTypeFragment = wtFrag

			// add some mutations
			mutations.add(ala.type)
			mutations.add(pro.type)

			// add some confs
			confs.addAll(wtFrag)
			confs.addAll(leu, "pp", "tp", "tt")
			confs.addAll(ala)
			confs.addAll(pro)

			// add some continuous flexibility
			for (frag in listOf(wtFrag, leu)) {
				val settings = DihedralAngle.LibrarySettings(
					radiusDegrees = 9.0,
					includeHydroxyls = true,
					includeNonHydroxylHGroups = false
				)
				for (space in confs.getByFragment(frag)) {
					space.motions.addAll(DihedralAngle.ConfDescription.makeFromLibrary(pos1, space.frag, space.conf, settings))
				}
			}
			for (space in confs.getByFragment(ala)) {
				val settings = DihedralAngle.LibrarySettings(
					radiusDegrees = 30.0,
					includeHydroxyls = true,
					includeNonHydroxylHGroups = true
				)
				space.motions.addAll(DihedralAngle.ConfDescription.makeFromLibrary(pos1, space.frag, space.conf, settings))
			}
		}
		positionConfSpaces.getOrMake(pos2).apply {

			val gly = conflib.fragments.getValue("GLY")

			// set only one mutation
			mutations.add(gly.type)

			// add the conf
			confs.addAll(gly, "GLY")

			// no continuous flexibility, it's glycine ...
		}

		// add a position to the small mol
		val pos3 = DesignPosition(
			name = "Pos3",
			type = "FOO",
			mol = smallmol
		).apply {

			// add a single anchor
			anchorGroups.add(mutableListOf(
				anchorSingle(
					smallmol.atoms.findOrThrow("C2"),
					smallmol.atoms.findOrThrow("C1"),
					smallmol.atoms.findOrThrow("C3")
				)
			))

			// add the source atom
			sourceAtoms.add(smallmol.atoms.findOrThrow("H10"))
		}
		designPositionsByMol[smallmol] = mutableListOf(pos3)

		// set the pos conf spaces
		positionConfSpaces.getOrMake(pos3).apply {

			// add the wt frag
			val wtFrag = pos3.makeFragment("wt-${pos3.name.toTomlKey()}", "WildType")
			wildTypeFragment = wtFrag

			// add a mutation
			mutations.add("BAR")

			// make a fragment from a hydroxyl group
			val ob = ConfLib.AtomInfo(1, "OB", Element.Oxygen)
			val hb = ConfLib.AtomInfo(2, "HB", Element.Hydrogen)
			val anchor = ConfLib.Anchor.Single(
				id = 1,
				bonds = listOf(ob)
			)
			val bar = ConfLib.Fragment(
				id = "bar",
				name = "Bar",
				type = "BAR",
				atoms = listOf(ob, hb),
				bonds = listOf(ConfLib.Bond(ob, hb)),
				anchors = listOf(anchor),
				confs = mapOf(
					"bar1" to ConfLib.Conf(
						id = "bar1",
						name = "Bar 1",
						description = null,
						// NOTE: coords borrowed from tyrosine hydroxyl group
						coords = identityHashMapOf(
							// none of the actual coords matter at all
							ob to Vector3d(-0.131, -6.026, 1.092),
							hb to Vector3d(0.132, -6.557, 1.849)
						),
						anchorCoords = identityHashMapOf(
							anchor to ConfLib.AnchorCoords.Single(
								// NOTE: don't make anchor atoms co-linear
								a = Vector3d(0.378, -4.766, 1.126), // CZ
								b = Vector3d(1.192, -4.367, 2.193), // CE1
								c = Vector3d(0.086, -3.866, 0.093) // CE2
							)
						)
					),
					"bar2" to ConfLib.Conf(
						id = "bar2",
						name = "Bar 2",
						description = "Bar 2 description, yup",
						// NOTE: coords borrowed from tyrosine hydroxyl group,
						// but HB is perturbed slightly
						coords = identityHashMapOf(
							// none of the actual coords matter at all
							ob to Vector3d(-0.131, -6.026, 1.092),
							hb to Vector3d(0.232, -6.457, 1.949)
						),
						anchorCoords = identityHashMapOf(
							anchor to ConfLib.AnchorCoords.Single(
								// NOTE: don't make anchor atoms co-linear
								a = Vector3d(0.378, -4.766, 1.126), // CZ
								b = Vector3d(1.192, -4.367, 2.193), // CE1
								c = Vector3d(0.086, -3.866, 0.093) // CE2
							)
						)
					)
				),
				motions = emptyList()
			)
			conflibs.add(ConfLib(
				id = "barlib",
				name = "BAR",
				fragments = mapOf(bar.id to bar)
			))

			// add some confs
			confs.addAll(wtFrag)
			confs.addAll(bar)

			// no motion settings, ie no continuous flexibility
		}

		// add a dihedral angle to the small molecule
		molMotions
			.getOrPut(smallmol) { ArrayList() }
			.add(DihedralAngle.MolDescription(
				smallmol,
				smallmol.atoms.findOrThrow("C6"),
				smallmol.atoms.findOrThrow("C1"),
				smallmol.atoms.findOrThrow("C"),
				smallmol.atoms.findOrThrow("N1"),
				0.0,
				60.0
			)
		)
	}
}

class TestConfSpace : FunSpec({

	infix fun ConfMotion.Description.shouldBeMotion(other: ConfMotion.Description) {

		// just make sure the types are correct for now
		this::class shouldBe other::class
	}

	infix fun MolMotion.Description.shouldBeMotion(other: MolMotion.Description) {

		// just make sure the types are correct for now
		this::class shouldBe other::class
	}

	infix fun ConfSpace.shouldBeConfSpace(exp: ConfSpace) {

		val obs = this

		obs.name shouldBe exp.name

		// check the conformation libraries
		obs.conflibs.map { it.id } shouldBe exp.conflibs.map { it.id }

		// check the molecules
		for ((obsPair, expPair) in obs.mols.zip(exp.mols)) {
			val (obsType, obsMol) = obsPair
			val (expType, expMol) = expPair

			obsType shouldBe expType
			obsMol.name shouldBe expMol.name
			obsMol.type shouldBe expMol.type

			// check the molecule motions
			val obsMotions = obs.molMotions.getOrDefault(obsMol, mutableListOf())
			val expMotions = exp.molMotions.getOrDefault(expMol, mutableListOf())
			obsMotions.size shouldBe expMotions.size
			for ((obsMotion, expMotion) in obsMotions.zip(expMotions)) {
				obsMotion shouldBeMotion expMotion
			}
		}

		// check the design positions
		obs.positions().size shouldBe exp.positions().size
		for ((obsPos, expPos) in obs.positions().zip(exp.positions())) {

			obsPos.name shouldBe expPos.name
			obsPos.type shouldBe expPos.type
			obsPos.mol.name shouldBe expPos.mol.name

			// check the anchor groups
			obsPos.anchorGroups.size shouldBe expPos.anchorGroups.size
			for ((obsGroup, expGroup) in obsPos.anchorGroups.zip(expPos.anchorGroups)) {

				obsGroup.size shouldBe expGroup.size
				for ((obsAnchor, expAnchor) in obsGroup.zip(expGroup)) {

					obsAnchor::class shouldBe expAnchor::class
					obsAnchor.anchorAtoms shouldBe expAnchor.anchorAtoms
				}
			}

			// check the atoms
			// re-do the sets to drop the indentity comparisons in favor of value comparisons
			obsPos.sourceAtoms.toSet() shouldBe expPos.sourceAtoms.toSet()

			// check the position conf spaces
			val obsPosConfSpace = obs.positionConfSpaces[obsPos]!!
			val expPosConfSpace = exp.positionConfSpaces[expPos]!!

			obsPosConfSpace.wildTypeFragment shouldBeFrag expPosConfSpace.wildTypeFragment
			if (expPosConfSpace.wildTypeFragment != null && obsPosConfSpace.wildTypeFragment != null) {
				obsPosConfSpace.wildTypeFragment shouldNotBeSameInstanceAs expPosConfSpace.wildTypeFragment
			}
			obsPosConfSpace.mutations shouldBe expPosConfSpace.mutations

			// check the confs
			obsPosConfSpace.confs.size shouldBe expPosConfSpace.confs.size
			for ((obsSpace, expSpace) in obsPosConfSpace.confs.zip(expPosConfSpace.confs)) {

				obsSpace.frag shouldBeFrag expSpace.frag
				obsSpace.conf shouldBeConf expSpace.conf

				val isWildtype = expSpace.frag === expPosConfSpace.wildTypeFragment
				if (isWildtype) {
					obsSpace.frag shouldNotBeSameInstanceAs expSpace.frag
				}

				// check the motions
				obsSpace.motions.size shouldBe expSpace.motions.size
				for ((obsMotion, expMotion) in obsSpace.motions.zip(expSpace.motions)) {
					obsMotion shouldBeMotion expMotion
				}
			}
		}
	}

	test("roundtrip") {

		// do the roundtrip
		val expConfSpace = makeTestConfSpace()
		val toml = expConfSpace.toToml()
		val obsConfSpace = ConfSpace.fromToml(toml)

		// make sure we got the same conf space back
		obsConfSpace shouldBeConfSpace expConfSpace
	}

	test("copy") {

		// do the roundtrip
		val expConfSpace = makeTestConfSpace()
		val obsConfSpace = expConfSpace.copy()

		// make sure we got the same conf space back
		obsConfSpace shouldBeConfSpace expConfSpace
	}

	test("copy protein") {

		val expConfSpace = makeTestConfSpace()
		val mol = expConfSpace.mols
			.mapNotNull { (type, mol) -> mol.takeIf { type == MoleculeType.Protein } }
			.first()
		val obsConfSpace = expConfSpace.copy(listOf(mol))

		// make sure we got just the protein
		obsConfSpace.mols shouldBe listOf(MoleculeType.Protein to mol)

		// and just the protein design positions
		obsConfSpace.positions().run {
			size shouldBe 2
			this[0].run {
				this.mol shouldBe mol
				name shouldBe expConfSpace.designPositionsByMol.getValue(mol)[0].name
			}
			this[1].run {
				this.mol shouldBe mol
				name shouldBe expConfSpace.designPositionsByMol.getValue(mol)[1].name
			}
		}
	}

	test("copy small mol") {

		val expConfSpace = makeTestConfSpace()
		val mol = expConfSpace.mols
			.mapNotNull { (type, mol) -> mol.takeIf { type == MoleculeType.SmallMolecule } }
			.first()
		val obsConfSpace = expConfSpace.copy(listOf(mol))

		// make sure we got just the small molecule
		obsConfSpace.mols shouldBe listOf(MoleculeType.SmallMolecule to mol)

		// and just the small molecule design positions
		obsConfSpace.positions().run {
			size shouldBe 1
			this[0].run {
				this.mol shouldBe mol
				name shouldBe expConfSpace.designPositionsByMol.getValue(mol)[0].name
			}
		}
	}

	test("can add per-molecule conflibs in a confspace") {
		val emptyConfLib = ConfLib("empty-conflib", "empty-conflib", mapOf())
		val lConflib = ConfLib.from(OspreyGui.getResourceAsString("conflib/lovell.conflib"))
		val dConflib = lConflib.invertChirality("D-${lConflib.id}", "D-${lConflib.name}")
		val protein = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0] as Polymer
		val leu26 = protein.findChainOrThrow("A").findResidueOrThrow("26")
		val pos1 = Proteins.makeDesignPosition(protein, leu26, "Pos1")

		ConfSpace(listOf(MoleculeType.Protein to protein)).apply {
			name = "confspace with invalid conflib"
			conflibs.add(emptyConfLib)
			designPositionsByMol[protein] = mutableListOf(pos1)
			positionConfSpaces.getOrMake(pos1).apply {
				shouldThrow<IllegalArgumentException> { // this should fail when there's no GLY fragments in the conflibs
					addConformationsFromLibraries(pos1, "GLY")
				}
			}
		}

		ConfSpace(listOf(MoleculeType.Protein to protein)).apply {
			name = "invalid conflib for protein specified so should fail"
			conflibs.addAll(listOf(emptyConfLib, lConflib))
			designPositionsByMol[protein] = mutableListOf(pos1)
			addConflibByMol(protein, emptyConfLib)
			positionConfSpaces.getOrMake(pos1).apply {
				shouldThrow<IllegalArgumentException> { // this should fail when there's no GLY fragments in the conflibs
					addConformationsFromLibraries(pos1, "GLY")
				}
			}
		}

		ConfSpace(listOf(MoleculeType.Protein to protein)).apply {
			name = "Conflibs with fragments with same ids but per-molecule conflib specified"
			conflibs.addAll(listOf(dConflib, lConflib))
			designPositionsByMol[protein] = mutableListOf(pos1)
			addConflibByMol(protein, lConflib)
			positionConfSpaces.getOrMake(pos1).apply {
				addConformationsFromLibraries(pos1, "ALA")
				confs.shouldNotBeEmpty()
				confs.size shouldBe lConflib.fragments["ALA"]?.confs?.size
			}
		}

		ConfSpace(listOf(MoleculeType.Protein to protein)).apply {
			name = "Conflibs with fragments with same types should be picked up (even if they're from different conflibs) unless conflibsByMol specified."
			conflibs.addAll(listOf(dConflib, lConflib))
			designPositionsByMol[protein] = mutableListOf(pos1)
			positionConfSpaces.getOrMake(pos1).apply {
				addConformationsFromLibraries(pos1, "VAL")
				val expectedConfs = (lConflib.fragments["VAL"]?.confs?.size ?: -1) + (lConflib.fragments["VAL"]?.confs?.size ?: -1)
				confs.size shouldBe expectedConfs
			}
		}
	}
})
