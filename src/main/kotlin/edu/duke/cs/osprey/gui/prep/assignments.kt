package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.MoleculeMaps
import edu.duke.cs.osprey.molscope.tools.associateIdentity
import edu.duke.cs.osprey.gui.io.ConfLib


class PosAssignment(
	val pos: DesignPosition,
	val frag: ConfLib.Fragment,
	val conf: ConfLib.Conf
) {
	override fun toString() =
		"PosAssignment[pos=${pos.name}, frag=${frag.id}, conf=${conf.id}]"
}


/**
 * Copies the molecules from the conformation space,
 * and set the conformations at each specified design position.
 */
class Assignments(val confSpace: ConfSpace, val assignments: List<PosAssignment>) {

	class MolInfo(
		/** a copy of the molecule from the conf space, with the new conformations applied */
		val assignedMol: Molecule,
		/** the original molecule from the conf space */
		val confSpaceMol: Molecule,
		/** molecule maps from the conf space molecule to the molecule copy in the assignments */
		val maps: MoleculeMaps
	) {

		fun getConfSpaceAtom(assignedAtom: Atom) =
			maps.atoms.getA(assignedAtom)

		fun getConfSpaceAtomOrThrow(assignedAtom: Atom) =
			getConfSpaceAtom(assignedAtom)
				?: throw NoSuchElementException("assigned atom $assignedAtom not found")

		fun getAssignedAtom(confSpaceAtom: Atom) =
			maps.atoms.getB(confSpaceAtom)

		fun getAssignedAtomOrThrow(confSpaceAtom: Atom) =
			getAssignedAtom(confSpaceAtom)
				?: throw NoSuchElementException("conf space atom $confSpaceAtom not found")
	}

	// copy the molecules
	val molInfos = confSpace.mols.map { (_, confSpaceMol) ->
		val (assignedMol, maps) = confSpaceMol.copyWithMaps()
		MolInfo(assignedMol, confSpaceMol, maps)
	}

	fun molInfoByConfSpaceMol(confSpaceMol: Molecule): MolInfo =
		molInfos.find { it.confSpaceMol === confSpaceMol }
			?: throw NoSuchElementException("no molecule info found for conf space molecule: $confSpaceMol")

	class AssignmentInfo(
		val molInfo: MolInfo,
		val confSwitcher: ConfSwitcher
	)

	// apply the assignments to the copied molecules
	val assignmentInfos: Map<PosAssignment,AssignmentInfo> = assignments
		.associateIdentity { assignment ->
			val molInfo = molInfoByConfSpaceMol(assignment.pos.mol)
			val confSwitcher = ConfSwitcher(assignment.pos, molInfo).apply {
				setConf(assignment.frag, assignment.conf)
			}
			assignment to AssignmentInfo(molInfo, confSwitcher)
		}

	fun assignmentInfo(pos: DesignPosition) =
		assignments.find { it.pos === pos }
}

fun ConfSpace.assign(assignments: List<PosAssignment>) = Assignments(this, assignments)
fun ConfSpace.assign(vararg assignments: PosAssignment) = Assignments(this, assignments.toList())
