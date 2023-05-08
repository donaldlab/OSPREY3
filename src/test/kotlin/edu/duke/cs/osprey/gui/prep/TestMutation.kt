package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.io.fromOMOL
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.show
import edu.duke.cs.osprey.molscope.molecule.*
import io.kotest.core.spec.style.FunSpec
import io.kotest.matchers.collections.shouldBeSameSizeAs
import io.kotest.matchers.collections.shouldContain
import io.kotest.matchers.comparables.beLessThanOrEqualTo
import io.kotest.matchers.should
import io.kotest.matchers.shouldBe
import io.kotest.matchers.types.shouldBeSameInstanceAs
import io.kotest.matchers.types.shouldBeTypeOf
import org.joml.Vector3d


class TestMutation : FunSpec({

	fun Vector3d.shouldBeNear(p: Vector3d) =
		distance(p) should beLessThanOrEqualTo(1e-3)

	fun Vector3d.shouldBeNear(x: Double, y: Double, z: Double) =
		shouldBeNear(Vector3d(x, y, z))

	fun Polymer.Residue.shouldHaveAtomNear(atomName: String, p: Vector3d) =
		findAtomOrThrow(atomName).pos.shouldBeNear(p)

	fun Polymer.Residue.shouldHaveAtomNear(atomName: String, x: Double, y: Double, z: Double) =
		shouldHaveAtomNear(atomName, Vector3d(x, y, z))

	fun Molecule.shouldBeConsistent() {

		val atomsLookup = atoms.toIdentitySet()
		val residueAtoms = if (this is Polymer) {
			chains.flatMap { chain -> chain.residues }
				.flatMap { residue -> residue.atoms }
				.toMutableSet()
		} else {
			mutableSetOf()
		}.toIdentitySet()

		// check that all residue atoms are also in the main atoms list
		(residueAtoms - atomsLookup) shouldBeSameSizeAs emptySet<Atom>()

		// check that all the bonds are to atoms also in the molecule
		val bondedAtoms = atoms.flatMap { bonds.bondedAtoms(it) }.toIdentitySet()
		(bondedAtoms - atomsLookup) shouldBeSameSizeAs emptySet<Atom>()
	}

	/**
	 * Use this to visually inspect the mutation, to see if the atoms and bonds look correct.
	 */
	fun DesignPosition.show() {
		mol.show(focusAtom=anchorGroups.first().first().anchorAtoms.first(), wait=true)
	}

	/**
	 * Then use this to dump the atom coords to build the regression test.
	 */
	fun Polymer.Residue.dump() {
		for (atom in atoms) {
			println("res.shouldHaveAtomNear(\"%s\", %.6f, %.6f, %.6f)"
				.format(atom.name, atom.pos.x, atom.pos.y, atom.pos.z))
		}
	}

	fun Polymer.Residue.capturePositions() =
		atoms.associate { it.name to Vector3d(it.pos) }


	context("protein") {

		val conflib = ConfLib.from(OspreyGui.getResourceAsString("conflib/lovell.conflib"))
		val protein1cc8 = Molecule.fromOMOL(OspreyGui.getResourceAsString("1cc8.protein.omol"))[0] as Polymer

		data class Instance(
			val res: Polymer.Residue,
			val pos: DesignPosition,
			val confSpace: ConfSpace
		) {
			constructor(res: Polymer.Residue, pos: DesignPosition) : this(
				res,
				pos,
				ConfSpace(listOf(MoleculeType.Protein to pos.mol)).apply {
					designPositionsByMol[pos.mol] = mutableListOf(pos)
				}
			)
		}
		fun ala2(): Instance {

			// copy the molecule, so we don't destroy the original
			val mol = protein1cc8.copy()
			val res = mol.findChainOrThrow("A").findResidueOrThrow("2")

			// make the design position
			val pos = Proteins.makeDesignPosition(mol, res, "ala2")

			// the current atoms should be glycine sidechain atoms
			pos.sourceAtoms.map { it.name }.toSet() shouldBe setOf("HA", "CB", "HB1", "HB2", "HB3")

			return Instance(res, pos)
		}
		fun gly17(): Instance {

			// copy the molecule, so we don't destroy the original
			val mol = protein1cc8.copy()
			val res = mol.findChainOrThrow("A").findResidueOrThrow("17")

			// make the design position
			val pos = Proteins.makeDesignPosition(mol, res, "gly17")
			
			// the current atoms should be glycine sidechain atoms
			pos.sourceAtoms.map { it.name }.toSet() shouldBe setOf("HA2", "HA3", "H")
			
			return Instance(res, pos)
		}
		fun pro52(): Instance {

			// copy the molecule, so we don't destroy the original
			val mol = protein1cc8.copy()
			val res = mol.findChainOrThrow("A").findResidueOrThrow("52")

			// make the design position
			val pos = Proteins.makeDesignPosition(mol, res, "pro52")

			// the current atoms should be proline sidechain atoms
			pos.sourceAtoms.map { it.name }.toSet() shouldBe setOf("CD", "HD2", "HD3", "CG", "HG2", "HG3", "CB", "HB2", "HB3", "HA")

			return Instance(res, pos)
		}

		fun DesignPosition.AnchorMatch.shouldHaveAnchors_CA_N(atoms: Set<Atom>, namesCA: Set<String>, namesN: Set<String>) {

			// should have gotten the 2x single anchors
			pairs.size shouldBe 2

			val (posAnchorCA, fragAnchorCA) = pairs[0]

			posAnchorCA.shouldBeTypeOf<Anchor.Single>()
			posAnchorCA.a.name shouldBe "CA"
			posAnchorCA.getConnectedAtoms(atoms).map { it.name }.toSet() shouldBe namesCA

			fragAnchorCA.shouldBeTypeOf<ConfLib.Anchor.Single>()
			fragAnchorCA.id shouldBe 1
			val (posAnchorN, fragAnchorN) = pairs[1]
			posAnchorN.shouldBeTypeOf<Anchor.Single>()
			posAnchorN.a.name shouldBe "N"
			posAnchorN.getConnectedAtoms(atoms).map { it.name }.toSet() shouldBe namesN
			fragAnchorN.shouldBeTypeOf<ConfLib.Anchor.Single>()
			fragAnchorN.id shouldBe 2
		}

		fun DesignPosition.AnchorMatch.shouldHaveAnchors_CA(atoms: Set<Atom>, namesCA: Set<String>) {

			// should have gotten the 1x single anchor
			pairs.size shouldBe 1

			val (posAnchorCA, fragAnchorCA) = pairs[0]
			posAnchorCA.shouldBeTypeOf<Anchor.Single>()
			posAnchorCA.a.name shouldBe "CA"
			posAnchorCA.getConnectedAtoms(atoms).map { it.name }.toSet() shouldBe namesCA
			fragAnchorCA.shouldBeTypeOf<ConfLib.Anchor.Single>()
			fragAnchorCA.id shouldBe 1
		}

		fun DesignPosition.AnchorMatch.shouldHaveAnchors_CAN(atoms: Set<Atom>, names: Set<String>) {

			// should have gotten the 1x double anchor
			pairs.size shouldBe 1

			val (posAnchor, fragAnchor) = pairs[0]
			posAnchor.shouldBeTypeOf<Anchor.Double>()
			posAnchor.a.name shouldBe "CA"
			posAnchor.b.name shouldBe "N"
			posAnchor.getConnectedAtoms(atoms).map { it.name }.toSet() shouldBe names
			fragAnchor.shouldBeTypeOf<ConfLib.Anchor.Double>()
			fragAnchor.id shouldBe 1
		}

		fun Assignments.AssignmentInfo.mapRes(res: Polymer.Residue) =
			(molInfo.maps as PolymerMaps).residues.getBOrThrow(res)

		// just spot-check a few amino acids

		test("glycine->glycine") {

			val (res, pos, confSpace) = gly17()

			// find the glycine conformation in the library
			val frag = conflib.fragments.getValue("GLY")
			val conf = frag.confs.getValue("GLY")

			// check the anchors
			pos.findAnchorMatch(frag)!!.shouldHaveAnchors_CA_N(
				atoms = pos.sourceAtoms,
				namesCA = setOf("HA2", "HA3"),
				namesN = setOf("H")
			)

			// do eeeet!
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// does this look like glycine?
			info.confSwitcher.type shouldBe "GLY"
			info.mapRes(res).apply {
				type shouldBe "GLY"
				atoms.size shouldBe 7
				shouldHaveAtomNear("N", 4.400000, -0.515000, 7.533000)
				shouldHaveAtomNear("CA", 5.791000, -0.751000, 7.871000)
				shouldHaveAtomNear("C", 6.672000, 0.451000, 7.612000)
				shouldHaveAtomNear("O", 7.716000, 0.674000, 8.236000)
				shouldHaveAtomNear("HA2", 5.865571, -1.028786, 8.932637)
				shouldHaveAtomNear("HA3", 6.167383, -1.604053, 7.287993)
				shouldHaveAtomNear("H", 3.958053, -0.848456, 6.691013)
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("glycine->valine") {

			val (res, pos, confSpace) = gly17()

			val frag = conflib.fragments.getValue("VAL")
			val conf = frag.confs.getValue("t")

			// check the anchors
			pos.findAnchorMatch(frag)!!.shouldHaveAnchors_CA_N(
				atoms = pos.sourceAtoms,
				namesCA = setOf("HA2", "HA3"),
				namesN = setOf("H")
			)

			// do eeeet!
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// does this look like valine?
			info.confSwitcher.type shouldBe "VAL"
			info.mapRes(res).apply {
				type shouldBe "VAL"
				atoms.size shouldBe 16
				shouldHaveAtomNear("N", 4.400000, -0.515000, 7.533000)
				shouldHaveAtomNear("CA", 5.791000, -0.751000, 7.871000)
				shouldHaveAtomNear("C", 6.672000, 0.451000, 7.612000)
				shouldHaveAtomNear("O", 7.716000, 0.674000, 8.236000)
				shouldHaveAtomNear("HA", 5.843315, -0.997520, 8.942108)
				shouldHaveAtomNear("CB", 6.366537, -1.934244, 7.054276)
				shouldHaveAtomNear("HB", 6.229395, -1.728927, 5.982360)
				shouldHaveAtomNear("CG1", 7.853488, -2.084830, 7.332262)
				shouldHaveAtomNear("CG2", 5.625587, -3.217240, 7.393457)
				shouldHaveAtomNear("HG11", 8.253188, -2.927268, 6.748425)
				shouldHaveAtomNear("HG12", 8.375506, -1.160221, 7.046118)
				shouldHaveAtomNear("HG13", 8.009019, -2.276337, 8.404558)
				shouldHaveAtomNear("HG21", 6.043071, -4.048968, 6.807168)
				shouldHaveAtomNear("HG22", 5.737316, -3.432624, 8.466278)
				shouldHaveAtomNear("HG23", 4.558468, -3.099393, 7.152487)
				shouldHaveAtomNear("H", 3.958053, -0.848456, 6.691013)
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("glycine->tryptophan") {

			val (res, pos, confSpace) = gly17()

			val frag = conflib.fragments.getValue("TRP")
			val conf = frag.confs.getValue("t90")

			// check the anchors
			pos.findAnchorMatch(frag)!!.shouldHaveAnchors_CA_N(
				atoms = pos.sourceAtoms,
				namesCA = setOf("HA2", "HA3"),
				namesN = setOf("H")
			)

			// do eeeet!
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// does this look like tryptophan?
			info.confSwitcher.type shouldBe "TRP"
			info.mapRes(res).apply {
				type shouldBe "TRP"
				atoms.size shouldBe 24
				shouldHaveAtomNear("N", 4.400000, -0.515000, 7.533000)
				shouldHaveAtomNear("CA", 5.791000, -0.751000, 7.871000)
				shouldHaveAtomNear("C", 6.672000, 0.451000, 7.612000)
				shouldHaveAtomNear("O", 7.716000, 0.674000, 8.236000)
				shouldHaveAtomNear("HA", 5.843864, -0.979317, 8.945833)
				shouldHaveAtomNear("CB", 6.346565, -1.946694, 7.093316)
				shouldHaveAtomNear("HB2", 5.705770, -2.821187, 7.279278)
				shouldHaveAtomNear("HB3", 6.288212, -1.728491, 6.016449)
				shouldHaveAtomNear("CG", 7.764233, -2.300739, 7.440858)
				shouldHaveAtomNear("CD1", 8.178192, -3.168377, 8.410530)
				shouldHaveAtomNear("CD2", 8.954448, -1.792557, 6.822732)
				shouldHaveAtomNear("HD1", 7.509462, -3.732152, 9.076669)
				shouldHaveAtomNear("NE1", 9.550764, -3.235475, 8.434387)
				shouldHaveAtomNear("CE2", 10.053614, -2.401242, 7.470743)
				shouldHaveAtomNear("CE3", 9.200747, -0.883001, 5.784103)
				shouldHaveAtomNear("HE1", 10.093753, -3.801440, 9.054872)
				shouldHaveAtomNear("CZ2", 11.378703, -2.130800, 7.114313)
				shouldHaveAtomNear("CZ3", 10.519567, -0.612321, 5.431754)
				shouldHaveAtomNear("HE3", 8.367570, -0.393433, 5.258512)
				shouldHaveAtomNear("HZ2", 12.219517, -2.616621, 7.631324)
				shouldHaveAtomNear("CH2", 11.590892, -1.237031, 6.095765)
				shouldHaveAtomNear("HZ3", 10.727176, 0.101402, 4.621522)
				shouldHaveAtomNear("HH2", 12.622241, -1.003771, 5.791199)
				shouldHaveAtomNear("H", 3.958053, -0.848456, 6.691013)
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("glycine->proline") {

			val (res, pos, confSpace) = gly17()

			val frag = conflib.fragments.getValue("PRO")
			val conf = frag.confs.getValue("up")

			// check the anchors
			pos.findAnchorMatch(frag)!!.shouldHaveAnchors_CAN(
				atoms = pos.sourceAtoms,
				names = setOf("HA2", "HA3", "H")
			)

			// do eeeet!
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// does this look like proline?
			// this part of the backbone has TOTALLY the wrong phi/psi conformation to support proline,
			// but otherwise, the mutation looks correct
			info.confSwitcher.type shouldBe "PRO"
			info.mapRes(res).apply {
				type shouldBe "PRO"
				atoms.size shouldBe 14
				shouldHaveAtomNear("N", 4.400000, -0.515000, 7.533000)
				shouldHaveAtomNear("CA", 5.791000, -0.751000, 7.871000)
				shouldHaveAtomNear("C", 6.672000, 0.451000, 7.612000)
				shouldHaveAtomNear("O", 7.716000, 0.674000, 8.236000)
				shouldHaveAtomNear("HG2", 5.006665, -3.076945, 5.560371)
				shouldHaveAtomNear("CB", 6.224608, -1.804685, 6.845759)
				shouldHaveAtomNear("HB3", 6.965117, -2.499234, 7.269964)
				shouldHaveAtomNear("HD2", 2.908340, -1.706927, 6.556454)
				shouldHaveAtomNear("HA", 5.837578, -1.102552, 8.915911)
				shouldHaveAtomNear("CG", 4.931207, -2.528918, 6.511270)
				shouldHaveAtomNear("HG3", 4.647473, -3.239651, 7.301579)
				shouldHaveAtomNear("CD", 3.948287, -1.376688, 6.416172)
				shouldHaveAtomNear("HB2", 6.653983, -1.333461, 5.948938)
				shouldHaveAtomNear("HD3", 4.019604, -0.850028, 5.452965)
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("glycine->proline->back") {

			val (res, pos, confSpace) = gly17()

			// record all the atom positions
			val wtPositions = res.capturePositions()

			// make the wildtype fragment
			val wtFrag = pos.makeFragment("wt", "WT")
			val wtConf = wtFrag.confs.values.first()

			// check the fragment anchors
			wtFrag.anchors.run {
				size shouldBe 2
				this[0].shouldBeTypeOf<ConfLib.Anchor.Single>()
				this[1].shouldBeTypeOf<ConfLib.Anchor.Single>()
			}

			// mutate to proline
			val frag = conflib.fragments.getValue("PRO")
			val conf = frag.confs.getValue("up")
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// mutate back to wildtype
			info.confSwitcher.setConf(wtFrag, wtConf)

			// does this look like the original residue?
			info.confSwitcher.type shouldBe "GLY"
			info.mapRes(res).apply {
				type shouldBe "GLY"
				atoms.size shouldBe wtPositions.size
				for ((atomName, atomPos) in wtPositions) {
					shouldHaveAtomNear(atomName, atomPos)
				}
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("glycine->valine->serine") {

			val (res, pos, confSpace) = gly17()

			// mutate to valine
			var frag = conflib.fragments.getValue("VAL")
			var conf = frag.confs.getValue("t")
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// mutate to serine
			frag = conflib.fragments.getValue("SER")
			conf = frag.confs.getValue("t_60")
			info.confSwitcher.setConf(frag, conf)

			// does this look like serine?
			info.confSwitcher.type shouldBe "SER"
			info.mapRes(res).apply {
				type shouldBe "SER"
				atoms.size shouldBe 11
				shouldHaveAtomNear("N", 4.400000, -0.515000, 7.533000)
				shouldHaveAtomNear("CA", 5.791000, -0.751000, 7.871000)
				shouldHaveAtomNear("C", 6.672000, 0.451000, 7.612000)
				shouldHaveAtomNear("O", 7.716000, 0.674000, 8.236000)
				shouldHaveAtomNear("HA", 5.839607, -0.973788, 8.946599)
				shouldHaveAtomNear("CB", 6.332726, -1.937088, 7.077001)
				shouldHaveAtomNear("HB2", 5.701465, -2.820324, 7.253298)
				shouldHaveAtomNear("HB3", 6.283402, -1.716092, 6.000542)
				shouldHaveAtomNear("OG", 7.667992, -2.225024, 7.443905)
				shouldHaveAtomNear("HG", 7.706109, -2.447422, 8.418113)
				shouldHaveAtomNear("H", 3.958053, -0.848456, 6.691013)
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("glycine->valine->back") {

			val (res, pos, confSpace) = gly17()

			// record all the atom positions
			val wtPositions = res.capturePositions()

			// make the wildtype fragment
			val wtFrag = pos.makeFragment("wt", "WT")
			val wtConf = wtFrag.confs.values.first()

			// check the fragment anchors
			wtFrag.anchors.run {
				size shouldBe 2
				this[0].shouldBeTypeOf<ConfLib.Anchor.Single>()
				this[1].shouldBeTypeOf<ConfLib.Anchor.Single>()
			}

			// mutate to valine
			val frag = conflib.fragments.getValue("VAL")
			val conf = frag.confs.getValue("t")
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// mutate back to wildtype
			info.confSwitcher.setConf(wtFrag, wtConf)

			// does this look like the original residue?
			info.confSwitcher.type shouldBe "GLY"
			info.mapRes(res).apply {
				type shouldBe "GLY"
				atoms.size shouldBe wtPositions.size
				for ((atomName, atomPos) in wtPositions) {
					shouldHaveAtomNear(atomName, atomPos)
				}
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("proline->glycine") {

			val (res, pos, confSpace) = pro52()

			// find the glycine conformation in the library
			val frag = conflib.fragments.getValue("GLY")
			val conf = frag.confs.getValue("GLY")

			// check the anchors
			pos.findAnchorMatch(frag)!!.shouldHaveAnchors_CA_N(
				atoms = pos.sourceAtoms,
				namesCA = setOf("HB2", "CD", "CB", "CG", "HG3", "HA", "HD2", "HG2", "HB3", "HD3"),
				namesN = setOf("HB2", "CD", "CB", "CG", "HG3", "HD2", "HG2", "HB3", "HD3")
			)

			// do eeeet!
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// does this look like glycine?
			info.confSwitcher.type shouldBe "GLY"
			info.mapRes(res).apply {
				type shouldBe "GLY"
				atoms.size shouldBe 7
				shouldHaveAtomNear("N", 17.783000, 20.016000, 15.734000)
				shouldHaveAtomNear("CA", 17.915000, 20.746000, 14.478000)
				shouldHaveAtomNear("C", 17.289000, 20.039000, 13.277000)
				shouldHaveAtomNear("O", 17.273000, 18.803000, 13.190000)
				shouldHaveAtomNear("HA2", 17.447047, 21.736009, 14.581452)
				shouldHaveAtomNear("HA3", 18.980680, 20.914459, 14.265382)
				shouldHaveAtomNear("H", 18.497148, 19.433047, 16.140951)
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("proline->glycine->back") {

			val (res, pos, confSpace) = pro52()

			// record all the atom positions
			val wtPositions = res.capturePositions()

			// make the wildtype fragment
			val wtFrag = pos.makeFragment("wt", "WT")
			val wtConf = wtFrag.confs.values.first()

			// check the fragment anchors
			wtFrag.anchors.run {
				size shouldBe 1
				val confAnchor = this[0]
				confAnchor.shouldBeTypeOf<ConfLib.Anchor.Double>()
				val confCoords = wtConf.anchorCoords[confAnchor]
				confCoords.shouldBeTypeOf<ConfLib.AnchorCoords.Double>()
				res.shouldHaveAtomNear("CA", confCoords.a)
				res.shouldHaveAtomNear("N", confCoords.b)
				protein1cc8.findChainOrThrow("A").findResidueOrThrow("51") // C in previous residue
					.shouldHaveAtomNear("C", confCoords.c)
				res.shouldHaveAtomNear("C", confCoords.d)
			}

			// mutate to glycine
			val frag = conflib.fragments.getValue("GLY")
			val conf = frag.confs.getValue("GLY")
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// mutate back to wildtype
			info.confSwitcher.setConf(wtFrag, wtConf)

			// does this look like the original residue?
			info.confSwitcher.type shouldBe "PRO"
			info.mapRes(res).apply {
				type shouldBe "PRO"
				atoms.size shouldBe wtPositions.size
				for ((atomName, atomPos) in wtPositions) {
					shouldHaveAtomNear(atomName, atomPos)
				}
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("Nala->glycine") {

			val (res, pos, confSpace) = ala2()

			// find the glycine conformation in the library
			val frag = conflib.fragments.getValue("GLYn")
			val conf = frag.confs.getValue("GLY")

			// check the anchors
			pos.findAnchorMatch(frag)!!.shouldHaveAnchors_CA(
				atoms = pos.sourceAtoms,
				namesCA = setOf("HA", "CB", "HB1", "HB2", "HB3")
			)

			// do eeeet!
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// does this look like glycine?
			info.confSwitcher.type shouldBe "GLY"
			info.mapRes(res).apply {
				type shouldBe "GLY"
				atoms.size shouldBe 9
				shouldHaveAtomNear("N", 14.789000, 27.073000, 24.130000)
				shouldHaveAtomNear("CA", 13.936000, 25.892000, 24.216000)
				shouldHaveAtomNear("C", 12.753000, 25.845000, 23.241000)
				shouldHaveAtomNear("O", 11.656000, 25.370000, 23.562000)
				shouldHaveAtomNear("H1", 14.423345, 27.671948, 23.403599)
				shouldHaveAtomNear("H2", 15.693612, 26.728407, 23.841847)
				shouldHaveAtomNear("H3", 14.686749, 27.541349, 25.018984)
				shouldHaveAtomNear("HA2", 14.542628, 24.991303, 24.041222)
				shouldHaveAtomNear("HA3", 13.526968, 25.811973, 25.233619)
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("Nala->glycine->back") {

			val (res, pos, confSpace) = ala2()

			// record all the atom positions
			val wtPositions = res.capturePositions()

			// make the wildtype fragment
			val wtFrag = pos.makeFragment("wt", "WT")
			val wtConf = wtFrag.confs.values.first()

			// check the fragment anchors
			wtFrag.anchors.run {
				size shouldBe 1
				this[0].shouldBeTypeOf<ConfLib.Anchor.Single>()
			}

			// mutate to glycine
			val frag = conflib.fragments.getValue("GLYn")
			val conf = frag.confs.getValue("GLY")
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// mutate back to wildtype
			info.confSwitcher.setConf(wtFrag, wtConf)

			// does this look like the original residue?
			info.confSwitcher.type shouldBe "ALA"
			info.mapRes(res).apply {
				type shouldBe "ALA"
				atoms.size shouldBe wtPositions.size
				for ((atomName, atomPos) in wtPositions) {
					shouldHaveAtomNear(atomName, atomPos)
				}
			}
			info.molInfo.assignedMol.shouldBeConsistent()
		}

		test("glycine->valine chi1") {

			val (_, pos, confSpace) = gly17()

			// mutate to valine
			val frag = conflib.fragments.getValue("VAL")
			val conf = frag.confs.getValue("t")
			val assignment = PosAssignment(pos, frag, conf)
			val info = confSpace.assign(assignment).assignmentInfos.getValue(assignment)

			// build the dihedral angle
			val chi1 = frag.motions[0] as ConfLib.ContinuousMotion.DihedralAngle
			DihedralAngle.ConfDescription(pos, chi1, conf, radiusDegrees = 9.0).run {
				make(info.molInfo.assignedMol, info.confSwitcher.atomResolverOrThrow).run {
					mol shouldBeSameInstanceAs info.molInfo.assignedMol
					a.name shouldBe "N"
					b.name shouldBe "CA"
					c.name shouldBe "CB"
					d.name shouldBe "CG1"
				}
			}
		}
	}

	// TODO: small molecules?
})
