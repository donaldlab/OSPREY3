package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.toIdentitySet
import edu.duke.cs.osprey.molscope.tools.identityHashSet
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.motions.ConfMotion
import edu.duke.cs.osprey.gui.motions.MolMotion
import java.math.BigInteger
import java.util.*


class ConfSpace(val mols: List<Pair<MoleculeType,Molecule>>) {

	var name = "Conformation Space"

	fun findMol(name: String): Molecule? =
		mols.map { (_, mol) -> mol }
			.find { it.name == name }

	val designPositionsByMol: MutableMap<Molecule,MutableList<DesignPosition>> = IdentityHashMap()

	/**
	 * Collect all the positions for all the molecules in a stable order.
	 */
	fun positions(): List<DesignPosition> =
		designPositionsByMol
			.toList()
			.sortedBy { (mol, _) -> mol.name }
			.flatMap { (_, positions) -> positions }

	class DuplicateConfLibException(val conflib: ConfLib) : RuntimeException("Conformation library already loaded: ${conflib.name}")

	class ConfLibs : Iterable<ConfLib> {

		private val conflibs = HashMap<String,ConfLib>()

		override fun iterator() = conflibs.values
			.sortedBy { it.name }
			.iterator()

		fun contains(conflib: ConfLib) =
			conflibs.containsKey(conflib.id)

		fun get(id: String) =
			conflibs[id]

		fun getOrThrow(id: String) =
			get(id) ?: throw NoSuchElementException("no conformation library with id $id")

		fun add(conflib: ConfLib) {

			// don't load the same library more than once
			if (contains(conflib)) {
				throw DuplicateConfLibException(conflib)
			}

			conflibs[conflib.id] = conflib
		}

		fun addAll(conflibs: Iterable<ConfLib>) {
			for (conflib in conflibs) {
				add(conflib)
			}
		}
	}
	val conflibs = ConfLibs()

	/**
	 * Collect all the unique library (ie non-wild-type) fragments used by the conf space.
	 */
	fun libraryFragments(): List<ConfLib.Fragment> =
		designPositionsByMol
			.values
			.flatten()
			.mapNotNull { positionConfSpaces[it] }
			.flatMap { posConfSpace ->
				posConfSpace.confs.fragments()
					.filter { it !== posConfSpace.wildTypeFragment }
			}
			.toCollection(identityHashSet())
			.sortedBy { it.id }

	/**
	 * Collect all the wild-type fragments from the conf space, whether used or not.
	 */
	fun wildTypeFragments(): List<ConfLib.Fragment> =
		designPositionsByMol
			.values
			.flatten()
			.mapNotNull { positionConfSpaces[it]?.wildTypeFragment }
			.sortedBy { it.id }

	class ConfConfSpace(val frag: ConfLib.Fragment, val conf: ConfLib.Conf) {

		val motions: MutableList<ConfMotion.Description> = ArrayList()
	}

	class PositionConfSpace {

		var wildTypeFragment: ConfLib.Fragment? = null
		val mutations: MutableSet<String> = HashSet()

		/**
		 * Returns true iff the position allows a sequence type other than the wildtype.
		 */
		fun isMutable() =
			mutations.any { it != wildTypeFragment?.type }

		inner class Confs : Iterable<ConfConfSpace> {

			private val byFragConf = IdentityHashMap<ConfLib.Fragment,MutableMap<ConfLib.Conf,ConfConfSpace>>()

			fun fragments() =
				byFragConf
					.keys
					.sortedBy { frag -> frag.id }

			override fun iterator() =
				fragments()
					.mapNotNull { frag -> byFragConf[frag] }
					.flatMap { spaces -> spaces.values.sortedBy { it.conf.id } }
					.iterator()

			val size get() =
				byFragConf
					.values
					.sumBy { it.size }

			fun get(frag: ConfLib.Fragment, conf: ConfLib.Conf) =
				byFragConf[frag]?.get(conf)

			fun getByFragment(frag: ConfLib.Fragment): List<ConfConfSpace> =
				byFragConf[frag]?.values?.toList() ?: emptyList()

			fun contains(frag: ConfLib.Fragment, conf: ConfLib.Conf) =
				get(frag, conf) != null

			fun getOrAdd(frag: ConfLib.Fragment, conf: ConfLib.Conf) =
				get(frag, conf) ?: add(frag, conf)

			fun add(frag: ConfLib.Fragment, conf: ConfLib.Conf): ConfConfSpace {

				// don't add duplicates
				if (get(frag, conf) != null) {
					throw IllegalStateException("position already has conformation ${frag.id}, ${conf.id}")
				}

				// add it
				val confConfSpace = ConfConfSpace(frag, conf)
				byFragConf
					.getOrPut(frag) { IdentityHashMap() }
					.put(conf, confConfSpace)
				return confConfSpace
			}

			fun addAll(frag: ConfLib.Fragment, confs: Iterable<ConfLib.Conf>) {
				for (conf in confs) {
					add(frag, conf)
				}
			}

			fun addAll(frag: ConfLib.Fragment, vararg confIds: String) =
				addAll(frag, frag.getConfs(*confIds))

			fun addAll(frag: ConfLib.Fragment) =
				addAll(frag, frag.confs.values)

			fun remove(frag: ConfLib.Fragment, conf: ConfLib.Conf) =
				byFragConf[frag]?.remove(conf)

			fun removeByFragmentType(type: String) =
				byFragConf.values.forEach { spaces ->
					spaces.values.removeIf { space -> space.frag.type == type }
				}
		}
		val confs = Confs()
	}
	class PositionConfSpaces {

		private val confSpaces: MutableMap<DesignPosition,PositionConfSpace> = IdentityHashMap()

		operator fun get(pos: DesignPosition) = confSpaces[pos]
		fun getOrMake(pos: DesignPosition) = confSpaces.getOrPut(pos) { PositionConfSpace() }
		fun remove(pos: DesignPosition) = confSpaces.remove(pos)

		fun sequenceSpaceSize(): BigInteger =
			confSpaces.values
				.takeIf { it.isNotEmpty() }
				?.map { it.mutations.size.toBigInteger() }
				?.reduce { a, b -> a.multiply(b) }
				?: BigInteger.ZERO

		fun confSpaceSize(): BigInteger =
			confSpaces.values
				.takeIf { it.isNotEmpty() }
				?.map { it.confs.size.toBigInteger() }
				?.reduce { a, b -> a.multiply(b) }
				?: BigInteger.ZERO
	}
	val positionConfSpaces = PositionConfSpaces()

	val molMotions = IdentityHashMap<Molecule,MutableList<MolMotion.Description>>()

	/**
	 * Get all the atoms that aren't part of design positions.
	 */
	fun fixedAtoms(mol: Molecule): List<Atom> {
		val posAtoms = (designPositionsByMol[mol]
			?.flatMap { it.sourceAtoms }
			?: emptyList())
			.toIdentitySet()
		return mol.atoms
			.filter { it !in posAtoms }
	}

	/**
	 * Makes a copy of the selected molecules along with all their
	 * design positions into a new conformation space
	 */
	fun copy(selMols: List<Molecule> = mols.map { (_, mol) -> mol }): ConfSpace {

		val old = this

		// match the selected molecules to the conf space
		val oldTypedMols = selMols.map { selMol ->
			old.mols
				.find { (_, mol) -> mol === selMol }
				?: throw NoSuchElementException("this conformation space didn't contain the selected molecule: $selMol")
		}
		val oldMols = oldTypedMols.map { (_, mol) -> mol }

		// make the new conf space from copies of the old molecules
		val new = ConfSpace(oldTypedMols.map { (type, mol) ->
			type to mol.copy()
		})
		val newMols = new.mols.map { (_, mol) -> mol }

		// copy the name
		new.name = old.name

		// copy over the conformation libraries
		new.conflibs.addAll(old.conflibs)

		for ((oldMol, newMol) in oldMols.zip(newMols)) {

			// make an atom map across the two molecules
			val oldToNew = oldMol.mapAtomsByValueTo(newMol)

			// copy the design positions, if needed
			old.designPositionsByMol[oldMol]?.let { oldPositions ->

				val newPositions = ArrayList<DesignPosition>()
				new.designPositionsByMol[newMol] = newPositions

				for (oldPos in oldPositions) {

					val newPos = DesignPosition(
						oldPos.name,
						oldPos.type,
						newMol
					)

					// copy the anchors
					newPos.anchorGroups.addAll(oldPos.anchorGroups.map { oldAnchors ->
						oldAnchors
							.map { oldAnchor ->
								oldAnchor.copyToMol(newMol, oldToNew)
							}
							.toMutableList()
					})

					// copy the current atoms
					newPos.sourceAtoms.addAll(oldPos.sourceAtoms.map { atom ->
						oldToNew.getBOrThrow(atom)
					})

					newPositions.add(newPos)

					// copy the position conf space
					val oldPosConfSpace = old.positionConfSpaces[oldPos] ?: continue
					val newPosConfSpace = new.positionConfSpaces.getOrMake(newPos)

					// copy the wild-type fragment, if needed
					oldPosConfSpace.wildTypeFragment?.let { oldWTFrag ->
						val oldWTConf = oldWTFrag.confs.values.first()
						newPosConfSpace.wildTypeFragment = newPos.makeFragment(
							fragId = oldWTFrag.id,
							fragName = oldWTFrag.name,
							confId = oldWTConf.id,
							confName = oldWTConf.name,
							motions = oldWTFrag.motions
						)

						// TODO: reset the atom ids to match the old atom ids?
					}

					// copy the mutations
					newPosConfSpace.mutations.addAll(oldPosConfSpace.mutations)

					for (oldSpace in oldPosConfSpace.confs) {

						// copy the fragment/conformation
						val isWildtype = oldSpace.frag === oldPosConfSpace.wildTypeFragment
						val newSpace = if (isWildtype) {
							val newWTFrag = newPosConfSpace.wildTypeFragment!!
							newPosConfSpace.confs.add(newWTFrag, newWTFrag.confs.values.first())
						} else {
							newPosConfSpace.confs.add(oldSpace.frag, oldSpace.conf)
						}

						// copy the motions
						for (oldMotion in oldSpace.motions) {
							newSpace.motions.add(oldMotion.copyTo(newPos))
						}
					}
				}
			}

			// copy the molecule motion settings
			val oldMotions = old.molMotions[oldMol]
			if (oldMotions != null) {
				new.molMotions[newMol] = oldMotions
					.map { it.copyTo(newMol) }
					.toMutableList()
			}
		}

		return new
	}

	companion object {
		// blank for now, but defined so it can be extended
	}
}
