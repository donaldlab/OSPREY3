package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType


/**
 * Generates chain ids in the range A-Z.
 */
open class ChainIdGeneratorAZ : ChainIdGenerator {

	private val usedIds = HashSet<String>()
	private var nextChainId = 'A'

	override fun setUsedIds(ids: Collection<String>) {
		usedIds.clear()
		usedIds.addAll(ids)
	}

	private fun getNextId(): String {
		if (nextChainId > 'Z') {
			throw IllegalStateException("out of unique chain ids in A-Z")
		}
		return "${nextChainId++}"
	}

	override fun generateId(): String {
		var chainId = getNextId()
		while (chainId in usedIds) {
			chainId = getNextId()
		}
		usedIds.add(chainId)
		return chainId
	}
}


/**
 * Generates single-residue chains for non-polymer molecules.
 * Useful for small molecules.
 */
class ChainGeneratorSingleResidue(val idGenerator: ChainIdGenerator) : ChainGenerator {

	override fun setUsedIds(ids: Collection<String>) =
		idGenerator.setUsedIds(ids)

	override fun generateChain(nonPolymerMol: Molecule, polymerMol: Polymer, polymerAtoms: List<Atom>) =
		Polymer.Chain(idGenerator.generateId()).apply {
			residues.add(Polymer.Residue(
				"1",
				nonPolymerMol.type ?: "SML",
				polymerAtoms
			))
		}
}


/**
 * Generates a chain with one residue for each group of bonded atoms in the molecule.
 * Useful for small molecules and solvents, as long as we have bonds.
 */
class ChainGeneratorBondedGroup(val idGenerator: ChainIdGenerator) : ChainGenerator {

	override fun setUsedIds(ids: Collection<String>) =
		idGenerator.setUsedIds(ids)

	override fun generateChain(nonPolymerMol: Molecule, polymerMol: Polymer, polymerAtoms: List<Atom>) =
		Polymer.Chain(idGenerator.generateId()).apply {

			var resNum = 1

			// get the connected components in the bond graph
			for (component in polymerMol.bonds.connectedComponents(polymerAtoms)) {

				residues.add(Polymer.Residue(
					(resNum++).toString(),
					nonPolymerMol.type ?: "SML",
					component.toList()
				))
			}
		}
}


class ChainGeneratorByMolType(val idGenerator: ChainIdGenerator) : ChainGenerator {

	override fun setUsedIds(ids: Collection<String>) =
		idGenerator.setUsedIds(ids)

	override fun generateChain(nonPolymerMol: Molecule, polymerMol: Polymer, polymerAtoms: List<Atom>): Polymer.Chain {

		val gen = when (nonPolymerMol.type?.let { MoleculeType[it] }) {
			MoleculeType.Solvent -> ChainGeneratorBondedGroup(idGenerator)
			else -> ChainGeneratorSingleResidue(idGenerator)
		}

		return gen.generateChain(nonPolymerMol, polymerMol, polymerAtoms)
	}
}


/**
 * A wrapper for Molecule.combine() that tries to intelligently
 * generate polymer chains for non-polymer molecules to satisfy
 * the limitations of the PDB file format.
 */
fun Collection<Molecule>.combineForPDB(name: String): Pair<Molecule,AtomMap> {
	val chainIdGen = ChainIdGeneratorAZ()
	val chainGen = ChainGeneratorByMolType(chainIdGen)
	return combine(name, chainIdGen, chainGen)
}
