package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer


/**
 * Maps from a Antechamber/LEaP-generated mol2 molecule to a Molscope molecule.
 */
class MoleculeMapper(val src: Molecule, val dst: Molecule) {

	// LEaP doesn't preserve chain ids or residue ids in the mol2 files,
	// so we need to translate based on residue order rather than id
	val srcResidues = (src as? Polymer)?.chains?.flatMap { it.residues }
	val dstResidues = (dst as? Polymer)?.chains?.flatMap { it.residues }

	// TODO: could speed this up by caching mappings?

	fun mapResidue(srcRes: Polymer.Residue): Polymer.Residue {

		srcResidues!!
		dstResidues!!

		if (srcResidues.size != dstResidues.size) {
			throw Error("Amber didn't preserve residue list")
		}

		return dstResidues[srcResidues.indexOf(srcRes)]
	}

	fun mapAtom(srcAtom: Atom): Atom? {

		return if (src is Polymer) {
			val srcRes = src.findResidue(srcAtom) ?: return null
			val dstRes = mapResidue(srcRes)
			dstRes.atoms.find { it.name == srcAtom.name }
		} else {
			dst.atoms.find { it.name == srcAtom.name }
		}
	}

	fun mapAtomOrThrow(srcAtom: Atom) =
		mapAtom(srcAtom) ?: throw NoSuchElementException("can't match atom $srcAtom")
}

fun Molecule.mapTo(dst: Molecule) = MoleculeMapper(this, dst)
