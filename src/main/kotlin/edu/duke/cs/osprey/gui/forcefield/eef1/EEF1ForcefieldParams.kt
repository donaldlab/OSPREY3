package edu.duke.cs.osprey.gui.forcefield.eef1

import cuchaz.kludge.tools.sqrt
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.forcefield.AtomIndex
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.forcefield.ForcefieldParams
import edu.duke.cs.osprey.gui.tools.ArrayMap
import edu.duke.cs.osprey.gui.tools.CombineCollisionException
import edu.duke.cs.osprey.gui.tools.combineMaps
import kotlin.math.PI
import kotlin.math.exp


class EEF1ForcefieldParams : ForcefieldParams {

	companion object {

		val trigConst = 2.0/(4.0*PI*PI.sqrt())

		/**
		 * Solvation interactions for atoms more than 9 A apart are already counted in dGref.
		 */
		const val cutoff = 9.0
	}

	override val forcefield = Forcefield.EEF1

	// EEF1 treats all 1-4 and farther as non-bonded
	override val unconnectedDistance = 4

	/**
	 * Scaling to apply to the solvent forcefield energy.
	 *
	 * The default value of 0.5 was determined empirically by Osprey developers
	 * to achieve good balance between the Amber96 and EEF1 forcefields.
	 */
	var scale = 0.5


	inner class AtomParams(
		val type: EEF1.AtomType
	) : ForcefieldParams.AtomParams {

		override fun internalEnergy(): Double? {
			return scale*type.dGref
		}

		override fun toString() = "$type"

		override fun hashCode() =
			type.hashCode()

		override fun equals(other: Any?) =
			other is AtomParams
				&& this.type == other.type
	}

	class AtomsParams : ForcefieldParams.AtomsParams {

		val atomsParams = ArrayMap<AtomParams>()

		override fun get(atomi: Int) = atomsParams[atomi]
	}

	override suspend fun parameterizeAtoms(mol: Molecule, atomIndex: AtomIndex, netCharge: Int?) =
		AtomsParams().apply {
			for (atom in mol.atoms) {
				val type = atom.atomTypeEEF1(mol) ?: continue
				val atomi = atomIndex.getOrThrow(atom)
				atomsParams[atomi] = AtomParams(type)
			}
		}


	inner class AtomPairParams(
		val vdwRadiusa: Double,
		val lambdaa: Double,
		val vdwRadiusb: Double,
		val lambdab: Double,
		val alpha1: Double,
		val alpha2: Double
	) : ForcefieldParams.AtomPairParams {

		constructor (a: EEF1.AtomType, b: EEF1.AtomType) : this(
			a.vdwRadius,
			a.lambda,
			b.vdwRadius,
			b.lambda,
			scale*trigConst*a.dGfree*b.volume/a.lambda,
			scale*trigConst*b.dGfree*a.volume/b.lambda
		)

		override val list = listOf(
			vdwRadiusa,
			lambdaa,
			vdwRadiusb,
			lambdab,
			alpha1,
			alpha2
		)

		override fun calcEnergy(r: Double): Double {
			return if (r <= cutoff) {
				val Xij = (r - vdwRadiusa)/lambdaa
				val Xji = (r - vdwRadiusb)/lambdab
				val r2 = r*r
				-(alpha1*exp(-Xij*Xij) + alpha2*exp(-Xji*Xji))/r2
			} else {
				0.0
			}
		}
	}

	inner class AtomPairsParams(
		val molsParams: Map<Int,AtomsParams>
	) : ForcefieldParams.AtomPairsParams {

		override fun get(moli1: Int, atomi1: Int, moli2: Int, atomi2: Int, dist: Int?): AtomPairParams? {

			// EEF1 only has interactions between atoms at least 1-4 bonded (ie >= 3 bonds away)
			// so skip the params if we're less than 3 bonds away
			if (dist != null && dist < 3) {
				return null
			}

			val params1 = molsParams[moli1]?.get(atomi1) ?: return null
			val params2 = molsParams[moli2]?.get(atomi2) ?: return null
			return AtomPairParams(params1.type, params2.type)
		}
	}

	override suspend fun parameterizeAtomPairs(infos: List<ForcefieldParams.MolInfo>) : AtomPairsParams {

		return AtomPairsParams(ArrayMap<AtomsParams>().apply {
			for (info in infos) {
				this[info.moli] = info.atomsParams as AtomsParams
			}
		})
	}

	override fun combineAtomsParams(info1: ForcefieldParams.AtomsInfo, info2: ForcefieldParams.AtomsInfo): ForcefieldParams.AtomsParams {

		info1.atomsParams as AtomsParams
		info2.atomsParams as AtomsParams

		val combined = AtomsParams()

		try {
			combineMaps(
				info1.atomsParams.atomsParams, info1.preferredAtomIndices, info1.ignoredAtomIndices,
				info2.atomsParams.atomsParams, info2.preferredAtomIndices, info2.ignoredAtomIndices,
				into = combined.atomsParams
			)
		} catch (ex: CombineCollisionException) {
			throw IllegalArgumentException("""
				|two EEF1 atoms params disagree on atom params for atom ${ex.key}:
				|	${ex.val1}
				|	${ex.val2}
			""".trimMargin())
		}

		return combined
	}
}
