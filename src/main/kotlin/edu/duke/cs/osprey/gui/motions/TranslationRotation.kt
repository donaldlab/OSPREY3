package edu.duke.cs.osprey.gui.motions

import edu.duke.cs.osprey.molscope.molecule.AtomMap
import edu.duke.cs.osprey.molscope.molecule.Molecule
import org.joml.Quaterniond
import org.joml.Vector3d


class TranslationRotation(val mol: Molecule): MolMotion {

	data class MolDescription(
		override val mol: Molecule,
		val maxTranslationDist: Double,
		val maxRotationDegrees: Double
	) : MolMotion.Description {

		override fun copyTo(mol: Molecule, atomMap: AtomMap) =
			MolDescription(mol, maxTranslationDist, maxRotationDegrees)

		override fun make() = TranslationRotation(mol)

		override fun getAffectedAtoms() = mol.atoms
	}

	// copy the coords
	val coords = mol.atoms
		.map { Vector3d(it.pos) }

	// compute the centroid
	val centroid = Vector3d().apply {
		for (atom in mol.atoms) {
			add(atom.pos)
		}
		div(mol.atoms.size.toDouble())
	}

	fun set(psi: Double, theta: Double, phi: Double, x: Double, y: Double, z: Double) {

		val qPsi = Quaterniond().rotationX(psi)
		val qTheta = Quaterniond().rotationY(theta)
		val qPhi = Quaterniond().rotationZ(phi)
		val t = Vector3d(x, y, z)

		// transform all the atoms, starting from the initial coords
		for (i in 0 until mol.atoms.size) {
			mol.atoms[i].pos.apply {
				set(coords[i])
				sub(centroid)
				rotate(qPsi)
				rotate(qTheta)
				rotate(qPhi)
				add(centroid)
				add(t)
			}
		}
	}

	fun reset() {
		for (i in 0 until mol.atoms.size) {
			mol.atoms[i].pos.apply {
				set(coords[i])
			}
		}
	}
}
