package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.tools.UnsupportedClassException
import org.joml.Vector3d
import java.util.*
import kotlin.NoSuchElementException


/**
 * A position in a molecule that allows removing the existing atoms and
 * aligning and attaching other conformations via an anchor atom system.
 */
class DesignPosition(
	var name: String,
	var type: String,
	val mol: Molecule
) {

	/**
	 * Anchor atoms are used to bond and align conformations to the molecule.
	 */
	val anchorGroups: MutableList<MutableList<Anchor>> = ArrayList()

	/**
	 * The atoms that will be replaced or re-located by the next conformation.
	 * Needs to be initialized when creating a new design position.
	 */
	val sourceAtoms: MutableSet<Atom> = Atom.identitySet()

	/**
	 * Returns true iff all the fragment's anchors are compatible with this design position.
	 */
	fun isFragmentCompatible(frag: ConfLib.Fragment) =
		findAnchorMatch(frag) != null

	fun compatibleFragments(conflib: ConfLib) =
		conflib.fragments
			.values
			.filter { isFragmentCompatible(it) }

	data class AnchorMatch(
		val posAnchors: List<Anchor>,
		val fragAnchors: List<ConfLib.Anchor>
	) {
		val pairs = posAnchors.zip(fragAnchors)

		fun findPosAnchor(fragAnchor: ConfLib.Anchor) =
			pairs
				.find { it.second === fragAnchor }
				?.first

		// NOTE: don't try to add a resolveCoords() function here
		// it won't work the way you want, since library conformation atoms are
		// in a different coordinate space than design position atoms
		// instead, use AtomPointer.resolveCoords(Conf)

		/** get the name from an atom pointer */
		fun resolveName(p: ConfLib.AtomPointer): String? =
			when (p) {
				is ConfLib.AtomInfo -> {
					// easy peasy
					p.name
				}
				is ConfLib.AnchorAtomPointer -> {
					// a bit more work ...
					val posAnchor = findPosAnchor(p.anchor)
						?: throw RuntimeException("no matched anchor")
					posAnchor.anchorAtoms.getOrNull(p.index)?.name
				}
				else -> throw UnsupportedClassException("unrecognized atom pointer type", this)
			}

		fun resolveNameOrThrow(p: ConfLib.AtomPointer) =
			resolveName(p) ?: throw NoSuchElementException("no atom found for pointer $p")

		/** get the element from an atom pointer */
		fun resolveElement(p: ConfLib.AtomPointer): Element? =
			when (p) {
				is ConfLib.AtomInfo -> {
					// easy peasy
					p.element
				}
				is ConfLib.AnchorAtomPointer -> {
					// a bit more work ...
					val posAnchor = findPosAnchor(p.anchor)
						?: throw RuntimeException("no matched anchor")
					posAnchor.anchorAtoms.getOrNull(p.index)?.element
				}
				else -> throw UnsupportedClassException("unrecognized atom pointer type", this)
			}

		fun resolveElementOrThrow(p: ConfLib.AtomPointer) =
			resolveElement(p) ?: throw NoSuchElementException("no atom found for pointer $p")

		fun copyToMol(mol: Molecule, atomMap: AtomMap) =
			AnchorMatch(
				posAnchors.map { it.copyToMol(mol, atomMap) },
				fragAnchors
			)
	}

	fun findAnchorMatch(frag: ConfLib.Fragment): AnchorMatch? =
		// check anchor types and order
		anchorGroups
			.firstOrNull { posAnchors ->

				// we should have the same number of anchors
				posAnchors.size == frag.anchors.size

				// check anchor types and order
				&& posAnchors
					.zip(frag.anchors)
					.all { (posAnchor, fragAnchor) ->
						posAnchor.isAnchorCompatible(fragAnchor)
					}

			}
			?.let { posAnchors ->
				AnchorMatch(posAnchors, frag.anchors)
			}

	class IncompatibleAnchorsException(val pos: DesignPosition, val frag: ConfLib.Fragment) :
		RuntimeException("design position ${pos.name} anchors are not compatible with fragment ${frag.id} anchors")

	/**
	 * Looks at the current atoms and the bonding pattern to see which anchor group is currently in use.
	 */
	fun getCurrentAnchorGroups() =
		anchorGroups
			.filter { anchors ->

				// get the connected atoms for each anchor
				val atoms = anchors.map { it.getConnectedAtoms(sourceAtoms) }

				for (i in 0 until atoms.size) {

					// connected atoms should be a single connected component
					if (!anchors[i].connectedAtomsIsSingleComponent(sourceAtoms)) {
						return@filter false
					}

					// no two pairs of connected atoms should overlap
					for (j in 0 until i) {
						if (atoms[i].intersection(atoms[j]).isNotEmpty()) {
							return@filter false
						}
					}
				}

				// all the atoms should exactly cover the source atoms
				return@filter atoms.union() == sourceAtoms
			}

	class IllegalAnchorsException(msg: String, val pos: DesignPosition) : IllegalStateException(msg)

	/**
	 * Makes a fragment from the source atoms and coords.
	 */
	fun makeFragment(
		fragId: String,
		fragName: String,
		confId: String = fragId,
		confName: String = fragName,
		motions: List<ConfLib.ContinuousMotion> = emptyList()
	): ConfLib.Fragment {

		val posAnchors = getCurrentAnchorGroups()
			.let { groups ->
				when {
					groups.isEmpty() -> throw IllegalAnchorsException("can't make fragment for this design position, no active anchors", this)
					groups.size > 1 -> throw IllegalAnchorsException("can't make fragment for this design position, multiple active anchors", this)
					else -> groups[0]
				}
			}

		// sort the atoms so we get a deterministic order from the set
		val sortedAtoms = sourceAtoms.sortedBy { it.name }

		// make the atom infos
		val atomInfos = Atom.identityMap<ConfLib.AtomInfo>().apply {
			for ((i, atom) in sortedAtoms.withIndex()) {
				put(atom, ConfLib.AtomInfo(i + 1, atom.name, atom.element))
			}
		}

		fun Atom.info() =
			atomInfos[this]
				?: throw NoSuchElementException("no atom info for $name")

		// translate the anchors
		var nextAnchorId = 1
		val anchorsFragByPos = posAnchors.associateWith { posAnchor ->
			posAnchor.makeLibraryAnchor(sourceAtoms, nextAnchorId++) { atom -> atom.info() }
		}

		val atoms = atomInfos.values
			.sortedBy { it.id }
			.toList()
		val anchors = anchorsFragByPos.values
			.sortedBy { it.id }

		return ConfLib.Fragment(
			id = fragId,
			name = fragName,
			type = type,
			atoms = atoms,
			bonds = sortedAtoms
				.flatMap { atom ->
					val atomInfo = atom.info()
					mol.bonds.bondedAtomsSorted(atom)
						.filter { it in sourceAtoms }
						.mapNotNull { bondedAtom ->
							val bondedInfo = bondedAtom.info()
							// keep only bonds in one direction
							if (atomInfo.id < bondedInfo.id) {
								ConfLib.Bond(atomInfo, bondedInfo)
							} else {
								null
							}
						}
				},
			anchors = anchors,
			confs = mapOf(
				confId to ConfLib.Conf(
					id = confId,
					name = confName,
					description = null,
					coords = sourceAtoms
						.associate { atom ->
							atom.info() to Vector3d(atom.pos)
						},
					anchorCoords = anchorsFragByPos
						.map { (posAnchor, fragAnchor) ->
							fragAnchor to posAnchor.makeLibraryCoords()
						}
						.associate { it }
				)
			),
			motions = motions
				.map { it.copyTo(atoms, anchors, it.id) }
		)
	}

	// convience shortcuts to emulate inner classes
	fun anchorSingle(a: Atom, b: Atom, c: Atom) = Anchor.Single(mol, a, b, c)
	fun anchorDouble(a: Atom, b: Atom, c: Atom, d: Atom) = Anchor.Double(mol, a, b, c, d)
}
