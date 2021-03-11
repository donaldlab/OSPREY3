package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.tools.UnsupportedClassException
import org.joml.Vector3d
import java.util.*


/**
 * Sets conformations on a molecule at a design position.
 */
class ConfSwitcher(
	val pos: DesignPosition,
	val molInfo: Assignments.MolInfo
) {

	// initialize the current atoms to the source atoms at the design position
	val currentAtoms =
		pos.sourceAtoms
			.map { molInfo.getAssignedAtomOrThrow(it) }
			.toIdentitySet()

	var type: String = pos.type
		private set

	var frag: ConfLib.Fragment? = null
		private set

	var conf: ConfLib.Conf? = null
		private set

	var anchorMatch: DesignPosition.AnchorMatch? = null
		private set

	var atomResolver: ConfLib.AtomPointer.Resolver? = null
		private set

	val atomResolverOrThrow get() =
		atomResolver ?: throw NoSuchElementException("design position currently has no atom resolver")

	/**
	 * Removes the current atoms from the molecule, including associated bonds.
	 * Then adds the mutated atoms to the molecule, adding bonds and aligning to the anchor as necessary.
	 *
	 * If the molecule is a polymer, the new atoms are added to the residues of their anchor atoms.
	 *
	 * Also set the atom pointer resolver.
	 */
	fun setConf(frag: ConfLib.Fragment, conf: ConfLib.Conf) {

		this.frag = frag
		this.conf = conf

		// find the anchor match and copy it to our molecule
		val anchorMatch = pos.findAnchorMatch(frag)
			?.copyToMol(molInfo.assignedMol, molInfo.maps.atoms)
			?: throw DesignPosition.IncompatibleAnchorsException(pos, frag)
		this.anchorMatch = anchorMatch

		// remove the existing atoms
		for (atom in currentAtoms) {
			molInfo.assignedMol.atoms.remove(atom)
			if (molInfo.assignedMol is Polymer) {
				molInfo.assignedMol.findResidue(atom)?.atoms?.remove(atom)
			}

			// update the atom map too
			molInfo.maps.atoms.removeB(atom)
		}
		currentAtoms.clear()

		// update the type
		type = frag.type

		// copy the atoms from the conf and add them to the molecule
		// update the current atoms with the newly added atoms
		val atomsByInfo = IdentityHashMap<ConfLib.AtomInfo,Atom>()
		for ((posAnchor, fragAnchor) in anchorMatch.pairs) {

			// find the residue, if any
			val res = posAnchor.findResidue()

			// update the residue type, if needed
			res?.type = frag.type

			for (atomInfo in frag.getAtomsFor(fragAnchor)) {

				// make the atom
				val atom = Atom(atomInfo.element, atomInfo.name, Vector3d(conf.coords[atomInfo]))
				atomsByInfo[atomInfo] = atom

				// add it to the mol
				molInfo.assignedMol.atoms.add(atom)
				res?.atoms?.add(atom)

				// update the current atoms
				currentAtoms.add(atom)
			}
		}

		// add the bonds
		for (bond in frag.bonds) {
			val atoma = atomsByInfo.getValue(bond.a)
			val atomb = atomsByInfo.getValue(bond.b)
			molInfo.assignedMol.bonds.add(atoma, atomb)
		}

		try {
			for ((posAnchor, fragAnchor) in anchorMatch.pairs) {

				// add the anchor bonds
				posAnchor.bondToAnchors(fragAnchor) { atomInfos, anchorAtom ->
					for (bondedInfo in atomInfos) {
						val atom = atomsByInfo.getValue(bondedInfo)
						molInfo.assignedMol.bonds.add(anchorAtom, atom)
					}
				}

				// align the conf coords
				posAnchor.align(currentAtoms, conf.anchorCoords.getValue(fragAnchor))
			}
		} catch (ex: IllegalAlignmentException) {

			// add extra information to the exception before passing it along
			throw IllegalConformationException(this, frag, conf, "Can't set conformation, bad alignment", ex)
		}

		// set the atom pointer resolver
		atomResolver = object : ConfLib.AtomPointer.Resolver {
			override fun resolve(p: ConfLib.AtomPointer) =
				when (p) {
					is ConfLib.AtomInfo -> {
						atomsByInfo[p]
					}
					is ConfLib.AnchorAtomPointer -> {
						anchorMatch
							.findPosAnchor(p.anchor)
							?.anchorAtoms
							?.get(p.index)
					}
					else -> throw UnsupportedClassException("unrecognized atom pointer type", p)
				}
		}
	}

	class IllegalConformationException(
		val switcher: ConfSwitcher,
		val frag: ConfLib.Fragment,
		val conf: ConfLib.Conf,
		msg: String,
		cause: Throwable? = null
	) : IllegalArgumentException("Design position ${switcher.pos.name} can't set conformation:\n\tfragment = ${frag.id}\n\tconformation = ${conf.id}\n$msg", cause)
}