package edu.duke.cs.osprey.gui.prep

import ags.utils.dataStructures.trees.thirdGenKD.KdTree
import ags.utils.dataStructures.trees.thirdGenKD.SquareEuclideanDistanceFunction
import cuchaz.kludge.tools.x
import cuchaz.kludge.tools.y
import cuchaz.kludge.tools.z
import edu.duke.cs.osprey.molscope.molecule.Atom
import org.joml.Vector3dc


/**
 * Builds a k-d tree that allows effecient nearest neighbor queries
 */
class AtomTree(val atoms: List<Atom>) {

	private val tree = KdTree<Atom>(3).apply {
		for (atom in atoms) {
			addPoint(doubleArrayOf(atom.pos.x, atom.pos.y, atom.pos.z), atom)
		}
	}

	private val distFunc = SquareEuclideanDistanceFunction()

	/**
	 * Iterates over atoms and their squared distances from the source position.
	 */
	fun nearest(pos: Vector3dc) = object : Iterator<Pair<Atom,Double>> {

		private val iter = tree.getNearestNeighborIterator(
			doubleArrayOf(pos.x, pos.y, pos.z),
			Int.MAX_VALUE,
			distFunc
		)

		override fun hasNext() = iter.hasNext()
		override fun next(): Pair<Atom,Double> {
			val atom = iter.next()
			val distSq = iter.distance()
			return atom to distSq
		}
	}

	fun nearest(atom: Atom) = nearest(atom.pos)
}

fun List<Atom>.toTree() = AtomTree(this)
