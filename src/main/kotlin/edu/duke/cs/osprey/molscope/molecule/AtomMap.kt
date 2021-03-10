package edu.duke.cs.osprey.molscope.molecule

import edu.duke.cs.osprey.molscope.tools.Bijection


/**
 * A bijection between two sets of atoms.
 */
class AtomMap : Bijection<Atom>() {

	companion object {

		fun identity(atoms: Collection<Atom>) =
			AtomMap().apply {
				for (atom in atoms) {
					add(atom, atom)
				}
			}
	}
}

class MoleculeMap : Bijection<Molecule>()

open class MoleculeMaps(
	val mols: MoleculeMap,
	val atoms: AtomMap
) {

	constructor (other: MoleculeMaps) : this(other.mols, other.atoms)
}
