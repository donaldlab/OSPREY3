package edu.duke.cs.osprey.gui.forcefield.amber

import org.joml.Vector3d


object CrdIO {

	/**
	 * Read an AMBER coordinates file (eg .crd)
	 */
	fun read(content: String): List<Vector3d> {

		val lines = content.lines()

		// get the number of atoms
		val numAtoms = lines[1].trim().toInt()

		// read the coords, fmt = 6F12.7
		val coords = ArrayList<Vector3d>(numAtoms)
		for (i in 2 until lines.size) {

			val line = lines[i]
			if (line.isBlank()) {
				continue
			}

			fun IntRange.toCoord(): Vector3d = this
				.map { line.substring(it*12, (it+1)*12) }
				.let { (sx, sy, sz) ->
					val x = sx.toDoubleOrNull()
					val y = sy.toDoubleOrNull()
					val z = sz.toDoubleOrNull()
					if (x == null || y == null || z == null) {
						throw ParseException("unrecognized coords: $sx, $sy, $sz on line ${i + 1}")
					}
					return@let Vector3d(x, y, z)
				}

			(0 .. 2).toCoord().let { coords.add(it) }
			if (line.length > 12*3) {
				(3 .. 5).toCoord().let { coords.add(it) }
			}
		}
		return coords
	}
}
