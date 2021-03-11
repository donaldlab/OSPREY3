package edu.duke.cs.osprey.gui.forcefield

import edu.duke.cs.osprey.molscope.molecule.Atom


/**
 * An bijection between a set of atoms and a set of integers.
 * The integers need not be contiguous.
 *
 * Useful to assign identites to atoms in molecules.
 */
class AtomIndex() {

	private val atomsToInts = Atom.identityMap<Int>()
	private val intsToAtoms = HashMap<Int,Atom>()

	constructor(atoms: List<Atom>) : this() {
		for ((atomi, atom) in atoms.withIndex()) {
			add(atomi, atom)
		}
	}

	val size: Int get() = atomsToInts.size

	fun add(atomi: Int, atom: Atom) {
		atomsToInts[atom] = atomi
		intsToAtoms[atomi] = atom
	}

	operator fun get(atom: Atom?) = atomsToInts[atom]

	fun getOrThrow(atom: Atom) =
		get(atom) ?: throw NoSuchElementException("doesn't have atom: $atom")

	operator fun get(atomi: Int?) = intsToAtoms[atomi]

	fun getOrThrow(atomi: Int) =
		get(atomi) ?: throw NoSuchElementException("doesn't have atom index: $atomi")

	operator fun set(atomi: Int, atom: Atom) = add(atomi, atom)
	operator fun set(atom: Atom, atomi: Int) = add(atomi, atom)

	operator fun contains(atom: Atom) = atomsToInts.containsKey(atom)
	operator fun contains(atomi: Int) = intsToAtoms.containsKey(atomi)

	override fun toString() =
		intsToAtoms.toString()
}
