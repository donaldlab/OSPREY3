package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.AtomMap
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.tools.assert
import edu.duke.cs.osprey.gui.forcefield.AtomIndex
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.toVector3d
import edu.duke.cs.osprey.service.services.MinimizeRequest
import org.joml.Vector3d
import java.util.*
import kotlin.NoSuchElementException
import kotlin.collections.ArrayList


class MinimizerInfo(
	val mol: Molecule,
	val restrainedAtoms: List<Atom> = emptyList()
) {

	// sander minimizations are currently configured to use implicit solvent
	// so make sure to filter out all the solvent molecules before minimizing

	// NOTE: Don't call any heavy-weight computations (like AmberTools) in this constructor,
	// so this class is safe to construct on any UI thread.
	// This class is just to prepare for heavy-weight computations later,
	// which should probably shouldn't run in UI threads.

	val partition: List<Pair<MoleculeType,Molecule>>
	val atomMap: AtomMap
	/** in the original molecule, not the partition */
	val minimizableAtoms: List<Atom>
	init {
		// split up the molecule into the parts amber expects
		val (partition, atomMap) = mol.partitionAndAtomMap(combineSolvent = false)

		// remove the solvent
		this.partition = partition.filter { (moltype, _) -> moltype != MoleculeType.Solvent }

		this.atomMap = atomMap

		// get the atoms in the original molecule that correspond to the non-solvent molecules
		minimizableAtoms = this.partition
			.flatMap { (_, mol) ->
				mol.atoms.map { atomMap.getAOrThrow(it) }
			}
	}

	var unminimizedCoords: List<Vector3d>? = null

	fun captureCoords() {
		unminimizedCoords = minimizableAtoms.map { Vector3d(it.pos) }
	}

	var minimizedCoords: List<Vector3d>? = null

	fun setCoords(coords: List<Vector3d>) {
		minimizableAtoms.forEachIndexed { i, atom ->
			atom.pos.set(coords[i])
		}
	}
}

suspend fun List<MinimizerInfo>.minimize(numSteps: Int) {

	// capture the original coords
	for (info in this) {
		info.captureCoords()
	}

	val infosByAtomIndex = ArrayList<Pair<IntRange,MinimizerInfo>>()
	val atomIndex = AtomIndex()

	// get the amber params for the combined molecules
	val params = this
		.flatMap { info ->
			info.partition
				.map { (moltype, mol) ->

					// but keep track of the atom indices,
					// so we can match the minimized coordinates back later
					val indexRange = atomIndex.size until atomIndex.size + mol.atoms.size
					infosByAtomIndex.add(indexRange to info)
					for (atom in mol.atoms) {
						atomIndex[atom] = atomIndex.size
					}

					val types = mol.calcTypesAmber(moltype, atomIndex)
					val frcmods = mol.calcModsAmber(types, atomIndex)
						?.let { listOf(it) }
						?: emptyList()
					AmberMolParams(mol, atomIndex, types, frcmods)
				}
		}
		.calcParamsAmber()

	// just in case, check the atom indices by matching atom names
	assert {
		TopIO.read(params.top).atomNames == flatMap { info -> info.minimizableAtoms.map { it.name } }
	}

	// convert restrained atoms into a sander-style restraint mask
	val restraintIndices = flatMap { info ->
		info.restrainedAtoms.map { atomA ->

			// map to the partitioned molecule
			val atomB = info.atomMap.getBOrThrow(atomA)

			// atom indices in sander start with 1
			atomIndex.getOrThrow(atomB) + 1
		}
	}
	val restraintMask = if (restraintIndices.isNotEmpty()) {
		"@"	+ restraintIndices.joinToString(",")
	} else {
		null
	}

	// minimize it!
	val allCoords = OspreyService.minimize(MinimizeRequest(
		params.top,
		params.crd,
		numSteps,
		restraintMask
	)).coords

	// grab the minimized coords for each mol
	val coordsByInfo = IdentityHashMap<MinimizerInfo,MutableList<Vector3d>>()
	allCoords.forEachIndexed { i, pos ->

		val info = infosByAtomIndex
			.find { (range, _) -> i in range }
			?.second
			?: throw NoSuchElementException("no MolInfo for minimized atom index $i")

		coordsByInfo.getOrPut(info) { ArrayList() }.add(pos.toVector3d())
	}

	// update the info with the minimized coords
	for ((info, coords) in coordsByInfo) {

		// just in case...
		assert(info.minimizableAtoms.size == coords.size)

		info.minimizedCoords = coords
	}
}
