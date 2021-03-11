package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.fromMol2
import edu.duke.cs.osprey.gui.io.toPDB
import edu.duke.cs.osprey.service.services.MissingAtomsRequest


/**
 * Uses Amber forecfields to infer missing heavy atoms and their positions.
 */
suspend fun Molecule.inferMissingAtomsAmber(): List<Pair<Atom,Polymer.Residue?>> {

	val dst = this
	val dstAtoms = ArrayList<Pair<Atom,Polymer.Residue?>>()

	// treat each molecule in the partition with the appropriate forcefield and ambertools
	partition@for ((type, src) in partition(combineSolvent = true)) {

		// TODO: allow user to pick the forcefields?
		val srcAtoms = when (type) {

			// treat regular molecules with leap
			MoleculeType.Protein,
			MoleculeType.DNA,
			MoleculeType.RNA -> runLeap(src, type.defaultForcefieldNameOrThrow)

			// ignore everything else
			else -> continue@partition
		}

		fun Polymer.Residue.translate(): Polymer.Residue {
			val srcRes = this
			val srcChain = (src as Polymer).chains.find { srcRes in it.residues }!!

			return (dst as Polymer).chains.find { it.id == srcChain.id }!!
				.residues.find { it.id == srcRes.id }!!
		}

		for ((atom, srcRes) in srcAtoms) {
			val dstRes = srcRes?.translate()
			dstAtoms.add(atom to dstRes)
		}
	}

	return dstAtoms
}


private suspend fun runLeap(mol: Molecule, ffname: ForcefieldName): List<Pair<Atom,Polymer.Residue?>> {

	// run LEaP to infer all the missing atoms
	val response = OspreyService.missingAtoms(MissingAtomsRequest(mol.toPDB(), ffname.name))
	val src = Molecule.fromMol2(response.mol2)

	val dst = mol
	val mapper = src.mapTo(dst)

	// find all the added heavy atoms
	val dstAtoms = ArrayList<Pair<Atom,Polymer.Residue?>>()
	if (src is Polymer) {

		for (srcChain in src.chains) {
			for (srcRes in srcChain.residues) {
				val dstRes = mapper.mapResidue(srcRes)

				srcRes.atoms
					.filter { atom -> atom.element != Element.Hydrogen }
					.filter { atom -> dstRes.atoms.none { it.name == atom.name } }
					.forEach { atom ->
						dstAtoms.add(atom to dstRes)
					}
			}
		}

	} else {
		throw UnsupportedOperationException("we only need this for polymers, right?")
	}

	return dstAtoms
}
