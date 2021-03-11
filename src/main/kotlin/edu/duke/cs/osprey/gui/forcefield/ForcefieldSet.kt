package edu.duke.cs.osprey.gui.forcefield


/** preserves add order though */
class ForcefieldSet private constructor(private val forcefields: ArrayList<ForcefieldParams>) : List<ForcefieldParams> by forcefields {
	
	constructor() : this(ArrayList())

	fun contains(ff: Forcefield): Boolean =
		forcefields.any { it.forcefield === ff }

	fun add(ff: Forcefield) =
		add(ff.parameterizer())

	fun add(ffparams: ForcefieldParams) {

		if (contains(ffparams)) {
			throw IllegalArgumentException("forcefield ${ffparams.forcefield} alread added")
		}

		forcefields.add(ffparams)
	}
}
