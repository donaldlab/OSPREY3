package edu.duke.cs.osprey.gui.compiler

import edu.duke.cs.osprey.gui.forcefield.ForcefieldParams
import java.util.*


/**
 * A convenient place to store molecule forcefield parameterizations
 * for an indexed conformation space
 */
class ConfSpaceParams(val index: ConfSpaceIndex) {

	private inner class ForForcefield {

		val wildTypes = ArrayList<ForcefieldParams.AtomsParams?>().apply {
			for (moli in index.mols.indices) {
				add(null)
			}
		}

		// pre-allocate enough storage for every position and fragment
		val frags: List<MutableList<ForcefieldParams.AtomsParams?>> =
			index.positions.map { posInfo ->
				posInfo.fragments
					.map {
						// IDEA is lying about the useless cast warning ...
						// the compiler apparently needs a little help figuring out this type
						@Suppress("USELESS_CAST")
						null as ForcefieldParams.AtomsParams?
					}
					.toMutableList()
			}
	}

	private val paramsByFF = IdentityHashMap<ForcefieldParams,ForForcefield>()

	private operator fun get(ff: ForcefieldParams) =
		paramsByFF.getOrPut(ff) { ForForcefield() }

	operator fun set(ff: ForcefieldParams, moli: Int, params: ForcefieldParams.AtomsParams) {
		this[ff].wildTypes[moli] = params
	}

	operator fun get(ff: ForcefieldParams, moli: Int): ForcefieldParams.AtomsParams =
		this[ff].wildTypes[moli]
			?: throw NoSuchElementException("no wild-type params for molecule ${index.mols[moli]} in forcefield $ff")

	operator fun set(ff: ForcefieldParams, fragInfo: ConfSpaceIndex.FragInfo, params: ForcefieldParams.AtomsParams) {
		this[ff].frags[fragInfo.posInfo.index][fragInfo.index] = params
	}

	operator fun get(ff: ForcefieldParams, fragInfo: ConfSpaceIndex.FragInfo): ForcefieldParams.AtomsParams =
		this[ff].frags[fragInfo.posInfo.index][fragInfo.index]
			?: throw NoSuchElementException("no params for fragment ${fragInfo.posInfo.pos.name} = ${fragInfo.frag.id} in forcefield $ff")
}
