package edu.duke.cs.osprey.molscope.tools

import org.joml.Vector3d
import java.util.ArrayList


/**
 * Sample points on a sphere quasi-uniformly by choosing the midpoionts of the faces of a regular icosohedron.
 * Subdivide faces of the icosohedron (and project back to the sphere) as needed to get the desired resolution.
 */
class SphereGrid(): Iterable<Vector3d> {

	class Face(a: Vector3d, b: Vector3d, c: Vector3d) {
		
		val vertices = arrayOf(a, b, c)

		val midpoint: Vector3d by lazy {
			Vector3d(vertices[0])
				.add(vertices[1])
				.add(vertices[2])
				.normalize()
		}
	}

	private var faces = ArrayList<Face>().apply {

		// compute the vertices of a regular icosahedron
		val x = 0.525731112119133606
		val z = 0.850650808352039932
		val vertices = listOf(
			Vector3d(-x, 0.0, z),
			Vector3d(x, 0.0, z),
			Vector3d(-x, 0.0, -z),
			Vector3d(x, 0.0, -z),
			Vector3d(0.0, z, x),
			Vector3d(0.0, z, -x),
			Vector3d(0.0, -z, x),
			Vector3d(0.0, -z, -x),
			Vector3d(z, x, 0.0),
			Vector3d(-z, x, 0.0),
			Vector3d(z, -x, 0.0),
			Vector3d(-z, -x, 0.0)
		)

		// these points are already pretty close to the unit sphere, but project them just to be sure
		vertices.forEach { it.normalize() }

		// create the faces
		add(Face(vertices[1], vertices[6], vertices[10]))
		add(Face(vertices[1], vertices[10], vertices[8]))
		add(Face(vertices[1], vertices[8], vertices[4]))
		add(Face(vertices[1], vertices[4], vertices[0]))
		add(Face(vertices[1], vertices[0], vertices[6]))

		add(Face(vertices[6], vertices[7], vertices[10]))
		add(Face(vertices[10], vertices[3], vertices[8]))
		add(Face(vertices[8], vertices[5], vertices[4]))
		add(Face(vertices[4], vertices[9], vertices[0]))
		add(Face(vertices[0], vertices[11], vertices[6]))

		add(Face(vertices[10], vertices[7], vertices[3]))
		add(Face(vertices[8], vertices[3], vertices[5]))
		add(Face(vertices[4], vertices[5], vertices[9]))
		add(Face(vertices[0], vertices[9], vertices[11]))
		add(Face(vertices[6], vertices[11], vertices[7]))

		add(Face(vertices[7], vertices[2], vertices[3]))
		add(Face(vertices[3], vertices[2], vertices[5]))
		add(Face(vertices[5], vertices[2], vertices[9]))
		add(Face(vertices[9], vertices[2], vertices[11]))
		add(Face(vertices[11], vertices[2], vertices[7]))
	}

	var subdivisions: Int = 0
		private set

	fun subdivide(times: Int) {
		for (i in 0 until times) {
			subdivide()
		}
	}

	fun subdivide() {

		val newFaces = ArrayList<Face>(faces.size*4)

		for (face in faces) {

			// get the vertices
			val a = face.vertices[0]
			val b = face.vertices[1]
			val c = face.vertices[2]

			// bisect the edges and project the midpoints to the unit sphere
			fun midpoint(s: Vector3d, t: Vector3d) =
				Vector3d(s)
					.add(t)
					.normalize()

			val ab = midpoint(a, b)
			val bc = midpoint(b, c)
			val ca = midpoint(c, a)

			// construct the new faces
			newFaces.add(Face(a, ab, ca))
			newFaces.add(Face(ab, bc, ca))
			newFaces.add(Face(ab, b, bc))
			newFaces.add(Face(ca, bc, c))
		}

		faces = newFaces
		subdivisions += 1
	}

	constructor (subdivisions: Int) : this() {
		subdivide(subdivisions)
	}

	override fun iterator(): Iterator<Vector3d> =
		faces.map { it.midpoint }.iterator()
}


// tool to generate sphere grids for use in e.g. the ambient occlusion shader
fun main() {

	// define some cardinal directions
	// NOTE: matches DIR_.. constants
	val dirs = listOf(
		Vector3d(-1.0, 0.0, 0.0), // 0
		Vector3d(+1.0, 0.0, 0.0), // 1
		Vector3d(0.0, -1.0, 0.0), // 2
		Vector3d(0.0, +1.0, 0.0), // 3
		Vector3d(0.0, 0.0, -1.0), // 4
		Vector3d(0.0, 0.0, +1.0)  // 5
	)

	val dirCounts = Array(dirs.size) { 0 }

	var numLines = 0
	for (v in SphereGrid(2)) {

		// which cardinal direction is closest?
		val iDir = dirs.indices.maxByOrNull { dirs[it].dot(v) }!!

		// filter out the points in "negative" directions
		if (iDir == 0 || iDir == 2 || iDir == 4) {
			continue
		}

		numLines += 1

		// apply the counts in both directions
		dirCounts[iDir] += 1
		dirCounts[iDir - 1] += 1

		// for LINES[]
		println("\t{ vec3(%8.5f,%8.5f,%8.5f), %1d },".format(v.x, v.y, v.z, iDir))
	}

	println("NUM_LINES = $numLines")

	// for NUM_RAYS_FOR_DIR[]
	for (i in dirs.indices) {
		println("\t${dirCounts[i]},")
	}
}
