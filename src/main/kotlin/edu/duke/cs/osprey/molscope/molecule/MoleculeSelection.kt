package edu.duke.cs.osprey.molscope.molecule


/** A function that returns a subset of atoms in a molecule */
typealias MoleculeSelector = (Molecule) -> List<Atom>

object MoleculeSelectors {

	/** Selects no atoms */
	val none: MoleculeSelector = { _ -> emptyList() }

	/** Selects all the atoms */
	val all: MoleculeSelector = { mol -> mol.atoms }

	/**
	 * If the molecule is a polymer, the selection returns the mainchain atoms.
	 * Otherwise, the selections returns no atoms.
	 */
	val mainchain: MoleculeSelector = { mol ->
		if (mol is Polymer) {
			// TODO: try to find the backbone using the bond network
			throw Error("not implemented yet")
		} else {
			emptyList()
		}
	}

	fun atomsByName(name: String): MoleculeSelector = { mol ->
		fun String.normalize() = uppercase()
		mol.atoms.filter { it.name.normalize() == name.normalize() }
	}
}


/**
 * Applies the selector and returns a new molecule containing only the selected atoms.
 * Polymer structure is preserved where possible
 */
fun MoleculeSelector.filter(mol: Molecule): Molecule {

	if (mol is Polymer) {

		val atomsLookup = this(mol).toIdentitySet()

		val filteredMol = Polymer(mol.name)
		for (chain in mol.chains) {

			val filteredChain = Polymer.Chain(chain.id)
			filteredMol.chains.add(filteredChain)

			for (res in chain.residues) {

				val filteredRes = Polymer.Residue(res.id, res.type)
				filteredChain.residues.add(filteredRes)

				val filteredAtoms = res.atoms.filter { it in atomsLookup }
				filteredRes.atoms.addAll(filteredAtoms)
				filteredMol.atoms.addAll(filteredAtoms)
			}
		}
		return filteredMol

	} else {

		val filteredMol = Molecule(mol.name, mol.type)
		filteredMol.atoms.addAll(this(mol))
		return filteredMol
	}
}