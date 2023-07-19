package edu.duke.cs.osprey.gui.compiler

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.tools.HashCalculator
import edu.duke.cs.osprey.gui.forcefield.AtomIndex
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.prep.*
import edu.duke.cs.osprey.gui.tools.ArrayMap
import kotlin.collections.ArrayList


/**
 * Establishes an authritative order for all the positions and conformations in the space.
 * Ignores any design positions that have no position conformation space.
 */
class ConfSpaceIndex(val confSpace: ConfSpace) {

	inner class ConfInfo(
		val posInfo: PosInfo,
		val fragInfo: FragInfo,
		val confConfSpace: ConfSpace.ConfConfSpace,
		val index: Int
	) {
		val conf get() = confConfSpace.conf

		val id = "${fragInfo.frag.id}:${conf.id}"
	}

	inner class FragInfo(
		val posInfo: PosInfo,
		val frag: ConfLib.Fragment,
		val index: Int
	) {

		/**
		 * Append the conf atoms to the dynamic fixed atoms for this design position,
		 * in a consistent order.
		 *
		 * Only return atoms from the assigned molecule.
		 */
		fun orderAtoms(fixedAtoms: FixedAtoms, assignmentInfo: Assignments.AssignmentInfo): List<Atom> {
			return ArrayList<Atom>().apply {

				// add the dynamic fixed atoms in the existing order
				for (atomInfo in fixedAtoms[posInfo].dynamics) {
					add(assignmentInfo.molInfo.getAssignedAtomOrThrow(atomInfo.atom))
				}

				// add the conformation atoms in the existing order
				for (atomInfo in frag.atoms) {
					add(assignmentInfo.confSwitcher.atomResolverOrThrow.resolveOrThrow(atomInfo))
				}
			}
		}
	}

	inner class PosInfo(
		val pos: DesignPosition,
		val moli: Int,
		val posConfSpace: ConfSpace.PositionConfSpace,
		val index: Int
	) {

		// index the fragments
		val fragments =
			posConfSpace.confs.fragments()
				.mapIndexed { i, frag -> FragInfo(this, frag, i) }

		// index the conformations
		val confs =
			fragments
				.flatMap { fragInfo ->
					posConfSpace.confs.getByFragment(fragInfo.frag)
						.sortedBy { space -> space.conf.id }
						.map { space -> fragInfo to space }
				}
				.mapIndexed { i, (fragInfo, space) ->
					ConfInfo(this, fragInfo, space, i)
				}

		/**
		 * Iterates over the conformations in the position conf space,
		 * returning assignments for each conformation in turn
		 */
		inline fun forEachConf(block: (Assignments, Assignments.AssignmentInfo, ConfInfo) -> Unit) {

			for (confInfo in confs) {

				// make the assignments
				val assignment = PosAssignment(pos, confInfo.fragInfo.frag, confInfo.conf)
				val assignments = confSpace.assign(assignment)
				val assignmentInfo = assignments.assignmentInfos.getValue(assignment)

				block(assignments, assignmentInfo, confInfo)
			}
		}

		/**
		 * Iterates over the fragments in the position conf space,
		 * returning assignments for an arbitrary conformation of each fragment in turn.
		 */
		inline fun forEachFrag(block: (Assignments, Assignments.AssignmentInfo, ConfInfo) -> Unit) {

			for (fragInfo in fragments) {

				// choose an arbitrary conformation from the fragment
				val confInfo = confs.first { it.fragInfo === fragInfo }

				// make the assignment
				val assignment = PosAssignment(pos, fragInfo.frag, confInfo.conf)
				val assignments = confSpace.assign(assignment)
				val assignmentInfo = assignments.assignmentInfos.getValue(assignment)

				block(assignments, assignmentInfo, confInfo)
			}
		}
	}

	// choose an order for the molecules
	val mols =
		confSpace.mols
			.map { (_, mol) -> mol }
			.sortedBy { it.name }

	fun findMoli(mol: Molecule): Int? =
		mols.indexOfFirst { it === mol }

	fun findMoliOrThrow(mol: Molecule): Int =
		findMoli(mol)
			?: throw NoSuchElementException("molecule $mol not indexed")

	// choose an order for the design positions and assign indices
	// alas, we need a named tuple here, since we need to associate three things together rather than just two
	private class PosTuple(val pos: DesignPosition, val moli: Int, val posConfSpace: ConfSpace.PositionConfSpace)

	val positions =
		mols.withIndex()
			.mapNotNull { (moli, mol) ->
				confSpace.designPositionsByMol[mol]
					?.mapNotNull { pos ->
						confSpace.positionConfSpaces[pos]?.let { PosTuple(pos, moli, it) }
					}
			}
			.flatten()
			.mapIndexed { i, t -> PosInfo(t.pos, t.moli, t.posConfSpace, i) }

	fun findPosi(pos: DesignPosition): Int? =
		positions
			.indexOfFirst { it.pos === pos }
			.takeIf { it >= 0 }

	fun getPosi(pos: DesignPosition): Int =
		findPosi(pos)
			?: throw NoSuchElementException("design position not found in this conf space")

	/**
	 * Lists of fixed atoms, indexed by molecule index
	 */
	val fixedAtoms: List<List<Atom>> = mols.map { confSpace.fixedAtoms(it) }


	// Assign an index to every atom in every conformation/fragment/position too,
	// so it's easy to associate forcefield parameters with them

	interface AtomKey

	class PositionAtomKey(
		val posi: Int,
		val frag: ConfLib.Fragment?,
		val atomi: Int
	) : AtomKey {

		// make sure to hash and compare by the *identity* of the fragment, not its value

		override fun hashCode() =
			HashCalculator.combineHashes(posi, System.identityHashCode(frag), atomi)

		override fun equals(other: Any?) =
			other is PositionAtomKey
				&& this.posi == other.posi
				&& this.frag === other.frag
				&& this.atomi == other.atomi
	}

	data class WildTypeAtomKey(
		val moli: Int,
		val atomi: Int
	) : AtomKey

	private val atomIndexWildTypeByMol = ArrayMap<AtomIndex>()
	private val atomIndexByKey = HashMap<AtomKey,Int>()

	private fun addKey(key: AtomKey): Int {
		val index = atomIndexByKey.size
		atomIndexByKey[key] = index
		return index
	}

	init {

		// wild-type atoms first
		for ((moli, mol) in mols.withIndex()) {

			val atomIndex = AtomIndex()
			atomIndexWildTypeByMol[moli] = atomIndex

			for ((atomi, atom) in mol.atoms.withIndex()) {
				atomIndex[atom] = addKey(WildTypeAtomKey(moli, atomi))
			}
		}

		// then designed atoms
		for ((posi, posInfo) in positions.withIndex()) {
			for (fragInfo in posInfo.fragments) {
				for (atomInfo in fragInfo.frag.atoms) {
					addKey(PositionAtomKey(posi, fragInfo.frag, atomInfo.id))
				}
			}
		}
	}

	/**
	 * Get the atom index for the wild type atoms in the given molecule.
	 */
	fun atomIndexWildType(moli: Int): AtomIndex =
		atomIndexWildTypeByMol[moli] ?: throw NoSuchElementException("no atom index exists for molecule $moli")

	/**
	 * Build an atom index for assigned molecules.
	 */
	fun matchAtoms(assignments: Assignments) : List<AtomIndex> {

		val atomIndices = mols.map { AtomIndex() }

		// start with the wild-type atoms
		for (moli in mols.indices) {
			val molInfo = assignments.molInfoByConfSpaceMol(mols[moli])
			for (atom in molInfo.assignedMol.atoms) {
				val csAtom = molInfo.getConfSpaceAtom(atom) ?: continue
				atomIndices[moli][atom] = atomIndexWildType(moli).getOrThrow(csAtom)
			}
		}

		// then add the designed atoms
		for ((assignment, info) in assignments.assignmentInfos) {
			val posi = getPosi(assignment.pos)
			val moli = positions[posi].moli
			for (atomInfo in assignment.frag.atoms) {
				val atom = info.confSwitcher.atomResolverOrThrow.resolveOrThrow(atomInfo)
				atomIndices[moli][atom] = atomIndexByKey.getValue(PositionAtomKey(posi, assignment.frag, atomInfo.id))
			}
		}

		return atomIndices
	}
}


inline fun ConfSpace.forEachFragIn(posInfo1: ConfSpaceIndex.PosInfo, posInfo2: ConfSpaceIndex.PosInfo, block: (Assignments, Assignments.AssignmentInfo, ConfSpaceIndex.ConfInfo, Assignments.AssignmentInfo, ConfSpaceIndex.ConfInfo) -> Unit) {

	for (fragInfo1 in posInfo1.fragments) {

		// choose an arbitrary conformation from the fragment
		val confInfo1 = posInfo1.confs.first { it.fragInfo === fragInfo1 }

		for (fragInfo2 in posInfo2.fragments) {

			// choose an arbitrary conformation from the fragment
			val confInfo2 = posInfo2.confs.first { it.fragInfo === fragInfo2 }

			// make the assignments
			val assignment1 = PosAssignment(posInfo1.pos, fragInfo1.frag, confInfo1.conf)
			val assignment2 = PosAssignment(posInfo2.pos, fragInfo2.frag, confInfo2.conf)
			val assignments = assign(assignment1, assignment2)
			val assignmentInfo1 = assignments.assignmentInfos.getValue(assignment1)
			val assignmentInfo2 = assignments.assignmentInfos.getValue(assignment2)

			block(assignments, assignmentInfo1, confInfo1, assignmentInfo2, confInfo2)
		}
	}
}
