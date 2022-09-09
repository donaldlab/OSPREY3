package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer


/**
 * Detects atom duplications in a molecule, and provides tools for removing one of the duplicated atoms.
 *
 * PDB files from the Protein Data Bank are usually very high quality files,
 * but rarely a third-party tool will create a malformed PDB file with duplicated atoms.
 * This tool can help correct those errors.
 */
class DuplicateAtoms private constructor(
	val mol: Molecule,
	private val groups: MutableList<AtomGroup>
) : List<DuplicateAtoms.AtomGroup> by groups {

	constructor(mol: Molecule) : this(mol, ArrayList())

	inner class AtomGroup(
		val name: String,
		val chain: Polymer.Chain?,
		val res: Polymer.Residue?,
		val atoms: List<Atom>
	) {
		fun included(atomi: Int) =
			mol.atoms.any { it === atoms[atomi] }

		fun add(atomi: Int) {

			// don't add it more than once
			if (included(atomi)) {
				return
			}

			mol.atoms.add(atoms[atomi])
			res?.atoms?.add(atoms[atomi])
		}

		fun remove(atomi: Int) {
			mol.atoms.remove(atoms[atomi])
			res?.atoms?.remove(atoms[atomi])
		}

		/** a friendly description for the location of the atom */
		val location: String? get() =
			if (res != null && chain != null) {
				"${chain.id}${res.id}"
			} else {
				null
			}

		override fun toString() =
			location
				?.let { "$name @ $it" }
				?: name
	}

	init {
		// look for duplicated atoms
		when (mol) {

			// check atoms for each residue
			is Polymer -> {
				for (chain in mol.chains) {
					for (res in chain.residues) {
						detectDupeAtoms(chain, res, res.atoms)
					}
				}
			}

			// check all the atoms at once
			else -> {
				detectDupeAtoms(null, null, mol.atoms)
			}
		}
	}

	private fun detectDupeAtoms(chain: Polymer.Chain?, res: Polymer.Residue?, atoms: List<Atom>) {

		// look for duplicates by atom name
		val dupesByName = HashMap<String,MutableList<Atom>>()
		for (atom in atoms) {
			dupesByName.getOrPut(atom.name) { ArrayList() }.add(atom)
		}

		// find the actual duplicates and flag them for user review
		for ((name, dupes) in dupesByName) {
			if (dupes.size > 1) {
				groups.add(AtomGroup(name, chain, res, dupes))
			}
		}
	}
}