package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.fromMol2
import edu.duke.cs.osprey.gui.io.toPDB
import edu.duke.cs.osprey.molscope.molecule.AtomPair
import edu.duke.cs.osprey.service.services.BondsRequest
import kotlinx.coroutines.runBlocking


fun Molecule.inferBondsAmberBlocking() =
	runBlocking { inferBondsAmber() }

/**
 * Uses Amber forecfields to infer atom connectivity,
 * but not bond order.
 */
suspend fun Molecule.inferBondsAmber(): List<AtomPair> {

	val dst = this
	val dstBonds = ArrayList<AtomPair>()

	// treat each molecule in the partition with the appropriate forcefield and ambertools
	val (partition, atomMap) = dst.partitionAndAtomMap(combineSolvent = true)
	partition@for ((type, src) in partition) {

		// TODO: allow user to pick the forcefields?
		// get the bonded mol from osprey service
		val request = when (type) {

			// treat molecules with either leap or antechamber
			MoleculeType.Protein,
			MoleculeType.DNA,
			MoleculeType.RNA,
			MoleculeType.Solvent -> BondsRequest(src.toPDB(), type.defaultForcefieldNameOrThrow.name)

			MoleculeType.SmallMolecule -> BondsRequest(src.toPDB(), null)

			// atomic ions don't have bonds
			MoleculeType.AtomicIon -> continue@partition

			// synthetics aren't real molecules, just ignore them
			MoleculeType.Synthetic -> continue@partition
		}
		val bondedMol = Molecule.fromMol2(OspreyService.bonds(request).mol2)

		// translate the bonds to the input mol
		val srcBonds = bondedMol.translateBonds(src)
		for ((srcA1, srcA2) in srcBonds) {
			dstBonds.add(AtomPair(
				atomMap.getAOrThrow(srcA1),
				atomMap.getAOrThrow(srcA2)
			))
		}
	}

	return dstBonds
}

private fun Molecule.translateBonds(dst: Molecule): List<Pair<Atom,Atom>> {

	val src: Molecule = this

	val mapper = src.mapTo(dst)

	// map the bonds from src to dst
	return src.bonds
		.toSet()
		.mapNotNull { (srcA1, srcA2) ->
			val dstA1 = mapper.mapAtom(srcA1) ?: return@mapNotNull null
			val dstA2 = mapper.mapAtom(srcA2) ?: return@mapNotNull null
			dstA1 to dstA2
		}
}
