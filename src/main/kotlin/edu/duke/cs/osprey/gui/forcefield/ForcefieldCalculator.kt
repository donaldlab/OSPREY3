package edu.duke.cs.osprey.gui.forcefield

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule


object ForcefieldCalculator {

	class MolInfo(
		moli: Int,
		mol: Molecule,
		atoms: Iterable<Atom>,
		atomIndex: AtomIndex,
		val atomsParams: ForcefieldParams.AtomsParams
	) : AtomPairer.MolInfo(moli, mol, atoms, atomIndex)

	/**
	 * Calculate the forcefield energy of the selected atoms.
	 *
	 * Called by the conf space compiler, so it doesn't need to be ultra fast, but it also needs to not be ultra slow.
	 */
	fun calc(atomPairsParams: ForcefieldParams.AtomPairsParams, infos: List<MolInfo>, unconnectedDistance: Int? = null): Double {

		var energy = 0.0

		// add the internal energies
		for (info in infos) {
			for (atom in info.atoms) {
				val atomi = info.atomIndex.getOrThrow(atom)
				energy += info.atomsParams[atomi]?.internalEnergy() ?: continue
			}
		}

		// add the atom pair energies
		for (molPair in AtomPairer.molPairs(infos)) {
			molPair.forEach(unconnectedDistance) { info1, atomi1, info2, atomi2, distance ->
				atomPairsParams[info1.moli, atomi1, info2.moli, atomi2, distance]?.let { params ->
					val atom1 = info1.atomIndex.getOrThrow(atomi1)
					val atom2 = info2.atomIndex.getOrThrow(atomi2)
					energy += params.calcEnergy(atom1.pos.distance(atom2.pos))
				}
			}
		}

		return energy
	}
}
