package edu.duke.cs.osprey.gui.compiler

import edu.duke.cs.osprey.gui.tools.pairs
import java.util.HashMap


class AtomPairs(index: ConfSpaceIndex) {

	data class AtomPair(
		val atomi1: Int,
		val atomi2: Int,
		val paramsi: Int
	)

	inner class Static(private val list: ArrayList<AtomPair> = ArrayList()) : List<AtomPair> by list {

		fun add(atomi1: Int, atomi2: Int, params: List<Double>) {
			list.add(AtomPair(atomi1, atomi2, paramsCache.index(params)))
		}
	}
	val static = Static()

	inner class PosSingles(index: ConfSpaceIndex) {

		// pre-allocate enough storage for every position and fragment
		private val list: List<List<MutableList<AtomPair>?>> =
			index.positions.map { posInfo ->
				posInfo.fragments
					.map {
						// IDEA is apparently lying about the warning here ...
						// looks like the compiler needs a little help figuring out this type
						@Suppress("RemoveExplicitTypeArguments")
						ArrayList<AtomPair>()
					}
			}

		fun add(posi1: Int, fragi1: Int, atomi1: Int, atomi2: Int, params: List<Double>) {
			list[posi1][fragi1]?.add(AtomPair(atomi1, atomi2, paramsCache.index(params)))
				?: throw NoSuchElementException("position:conformation $posi1:$fragi1 is not in this conf space")
		}

		operator fun get(posi1: Int, fragi1: Int): List<AtomPair> =
			list[posi1][fragi1]
				?: throw NoSuchElementException("position:conformation $posi1:$fragi1 is not in this conf space")
	}
	val pos = PosSingles(index)
	val posStatic = PosSingles(index)

	inner class PosPairs(index: ConfSpaceIndex) {

		// pre-allocate enough storage for every position pair and fragment pair
		private val list: List<List<List<MutableList<AtomPair>?>>> =
			index.positions.pairs()
				.map { (posInfo1, posInfo2) ->
					posInfo1.fragments.map {
						posInfo2.fragments.map {
							// IDEA is apparently lying about the warning here ...
							// looks like the compiler needs a little help figuring out this type
							@Suppress("RemoveExplicitTypeArguments")
							ArrayList<AtomPair>()
						}
					}
				}

		private fun posIndex(posi1: Int, posi2: Int): Int {
			if (posi2 < posi1) {
				return posi1*(posi1 - 1)/2 + posi2
			} else {
				throw IllegalArgumentException("posi2 $posi2 must be strictly less than posi1 $posi1")
			}
		}

		fun add(posi1: Int, fragi1: Int, posi2: Int, fragi2: Int, atomi1: Int, atomi2: Int, params: List<Double>) {
			list[posIndex(posi1, posi2)][fragi1][fragi2]?.add(AtomPair(atomi1, atomi2, paramsCache.index(params)))
				?: throw NoSuchElementException("position:conformation pair $posi1:$fragi2 - $posi2:$fragi2 is not in this conf space")
		}

		operator fun get(posi1: Int, fragi1: Int, posi2: Int, fragi2: Int): List<AtomPair> =
			list[posIndex(posi1, posi2)][fragi1][fragi2]
				?: throw NoSuchElementException("position:conformation pair $posi1:$fragi2 - $posi2:$fragi2 is not in this conf space")
	}
	val posPos = PosPairs(index)

	/**
	 * Different atom pairs will often yield the exact same forcefield parameters.
	 * This cache allows us to de-duplicate the params and save quite a bit of space/time.
	 */
	class ParamsCache(
		private val list: ArrayList<List<Double>> = ArrayList()
	) : List<List<Double>> by list {

		private val indices = HashMap<List<Double>,Int>()

		fun index(params: List<Double>): Int {

			// check the cache first
			indices[params]?.let { return it }

			// cache miss, add the params
			val index = list.size
			list.add(params)
			indices[params] = index
			return index
		}
	}
	val paramsCache = ParamsCache()
}
