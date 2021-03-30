package edu.duke.cs.osprey.gui.forcefield.amber

import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.forcefield.AtomIndex
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.forcefield.ForcefieldParams
import edu.duke.cs.osprey.gui.tools.ArrayMap
import edu.duke.cs.osprey.service.ServiceException
import java.util.*


class Amber96Params : AmberForcefieldParams(mapOf(
	MoleculeType.Protein to ForcefieldName.ff96
)) {
	override val forcefield = Forcefield.Amber96
}

class Amber14SBParams : AmberForcefieldParams(mapOf(
	MoleculeType.Protein to ForcefieldName.ff14SB
)) {
	override val forcefield = Forcefield.Amber14SB
}


object AmberDefaults {

	/**
	 * This appears to be a magic number (for now).
	 * Maybe 6 is reasonable for a protein interior?
	 *
	 * For reference, this value seems to have a simlar relative permittivity
	 * to neoperene (6.7), but is far less than water (~80).
	 * See: https://en.wikipedia.org/wiki/Relative_permittivity
	 */
	const val dielectric = 6.0

	const val distanceDependentDielectric = true

	/**
	 * The default value of 0.95 was determined empirically by early Osprey developers.
	 */
	const val vdwScale = 0.95

	/**
	 * The default AM1BCC method is currently recommended by Amber for most purposes.
	 */
	val chargeMethod = AmberChargeMethod.AM1BCC

	/**
	 * The default value of 0 assumes the input structure should not be perturbed via minimization.
	 * ie, 0 minimization steps results in a single point energy calculation.
	 */
	const val sqmMinimizationSteps = 0
}

abstract class AmberForcefieldParams(val ffnameOverrides: Map<MoleculeType,ForcefieldName>) : ForcefieldParams {

	// Amber has special non-bonded parameters for 1-4 bonds
	// anything farther is treated as generally non-bonded
	override val unconnectedDistance = 5

	/**
	 * The dielectric constant of the environment (aka its relative permittivity),
	 * which influences electrostatic calculations.
	 */
	var dielectric = AmberDefaults.dielectric

	/**
	 * If true, multiply the dielectric contstant by the atom pair distance (r)
	 * for electrostatic interactions.
	 */
	var distanceDependentDielectric = AmberDefaults.distanceDependentDielectric

	/**
	 * Scaling to apply to the van der Waals calculations.
	 */
	var vdwScale = AmberDefaults.vdwScale

	/**
	 * Method to generate partial charges for small molecules.
	 */
	var chargeMethod = AmberDefaults.chargeMethod

	/**
	 * Number of steps of `xmin` minimization to perform during partial charge calculation
	 * for small molecules.
	 *
	 * Minimization with `xmin` is slow, so adding minimization steps can add large amounts of time
	 * to your conformation space compilations.
	 */
	var sqmMinimizationSteps = 0


	override fun settings() = LinkedHashMap<String,Any>().apply {
		this["distanceDependentDielectric"] = distanceDependentDielectric
	}

	// no internal energies for Amber forcefields, so just re-use the same instance every time we need an AtomParams
	class AtomParams(val type: String, val charge: String) : ForcefieldParams.AtomParams {

		/** Amber forcefields don't have internal energies */
		override fun internalEnergy() = null

		override fun toString() = "$type:$charge"

		override fun hashCode() =
			type.hashCode()

		override fun equals(other: Any?) =
			other is AtomParams
				&& this.type == other.type
				&& this.charge == other.charge
	}

	class AtomsParams(
		val types: AmberTypes,
		val frcmods: List<String>
	) : ForcefieldParams.AtomsParams {

		override fun get(atomi: Int): AtomParams? {
			val atomType = types.atomTypes[atomi]
			val atomCharge = types.atomCharges[atomi]
			if (atomType == null || atomCharge == null) {
				return null
			}
			return AtomParams(atomType, atomCharge)
		}
	}

	override suspend fun parameterizeAtoms(mol: Molecule, atomIndex: AtomIndex, netCharge: Int?): AtomsParams {

		// find the forcefield name
		val molType = mol.findTypeOrThrow()
		val ffname = ffnameOverrides[molType] ?: molType.defaultForcefieldNameOrThrow

		// get charge generation settings if small molecule
		val generateCharges = when (molType) {
			MoleculeType.SmallMolecule -> AmberChargeGeneration(
				chargeMethod,
				netCharge ?: throw IllegalArgumentException("net charge needed to generate partial charges for small molecule"),
				sqmMinimizationSteps
			)
			else -> null
		}

		// get the amber types for the molecule
		val types = try {

			mol.calcTypesAmber(
				mol.findTypeOrThrow(),
				atomIndex,
				ffname,
				generateCharges
			)

		} catch (ex: ServiceException) {

			// TEMP: just re-throw the exception
			throw ex

			/* TODO: try to diagnose problems and give helpful hints to the user

			fun List<Pair<SQM.ErrorType,String>>.format(): String {
				val buf = StringBuffer()
				for ((errType, errMsg) in this) {

					buf.append("\n")

					// show the error message
					errMsg
						.split("\n")
						.forEachIndexed { i, line ->
							if (i == 0) {
								buf.append(" * $line\n")
							} else {
								buf.append("   $line\n")
							}
						}

					// suggest potential fixes
					val suggestedFix = when (errType) {
						SQM.ErrorType.NoConvergence ->
							"Maybe try fixing issues with molecule chmesitry or structure?"
						SQM.ErrorType.BadNetCharge ->
							"Try changing the net charge for this conformation?"
					}
					buf.append("TO FIX: $suggestedFix\n")
				}
				return buf.toString()
			}

			// parameterizing small molecules is especially error-prone
			// since we need to call SQM to compute the partial charges
			// if SQM failed, try to give a friendly(er) error message to the user
			ex.results.sqm
				?.takeIf { it.errors.isNotEmpty() }
				?.let {
					throw RuntimeException("""
							|Can't generate partial charges for $mol
							|SQM failed with errors:
							|${it.errors.format()}
						""".trimMargin(), ex)
				}
				// if we got here, we didn't parse any errors from sqm, so just re-throw the original exception
				?: throw ex
			*/
		}

		// calculate the frcmod, if needed
		val frcmods = mol.calcModsAmber(types, atomIndex)
			?.let { listOf(it) }
			?: emptyList()

		return AtomsParams(types, frcmods)
	}

	inner class AtomPairParams(
		val esQ: Double,
		val vdwA: Double,
		val vdwB: Double
	) : ForcefieldParams.AtomPairParams {

		// TODO: should we include torsion parameters?!

		override val list = listOf(
			esQ,
			vdwA,
			vdwB
		)

		/**
		 * See Amber manual, Eqn 14.1
		 * But we've modified the electrostatic calculations to use a
		 * distance-dependent dielectric (ie 1/r2 instead of 1/r) when needed.
		 */
		override fun calcEnergy(r: Double): Double {

			// calculate the electrostatics energy
			val r2 = r*r
			val es = if (distanceDependentDielectric) {
				esQ/r2
			} else {
				esQ/r
			}

			// calculate the van der Waals energy
			val r6 = r2*r2*r2
			val r12 = r6*r6
			val vdw = vdwA/r12 - vdwB/r6

			return es + vdw
		}
	}

	inner class AtomPairsParams(
		val top: AmberTopology.Mapped,
		val indexToTop: Map<Int,Map<Int,Int>>
	) : ForcefieldParams.AtomPairsParams {

		override fun get(moli1: Int, atomi1: Int, moli2: Int, atomi2: Int, dist: Int?): AtomPairParams? {

			// Amber forcefields only have "non-bonded" interactions between atoms at least 1-4 bonded (ie >= 3 bonds away)
			// so skip the params if we're less than 3 bonds away
			if (dist != null && dist < 3) {
				return null
			}

			// map the mol,atom indces to the topology info
			val ia = indexToTop[moli1]?.get(atomi1) ?: return null
			val ib = indexToTop[moli2]?.get(atomi2) ?: return null

			// get the electrostatic params
			// NOTE: the coulomb factor (~322.05) is already pre-multiplied into the charges
			// see: Amber manual 14.1.7, Partial Charges
			var esQ = (top.charge(ia) * top.charge(ib))/dielectric

			// get the van der Waals params
			// NOTE: Lorentz/Berthelot mixing rules already pre-applied to A and B
			// see: Amber manual 14.1.7, Van der Waals Parameters
			var (vdwA, vdwB) = top.vdw(ia, ib)

			/* Now, to apply Osprey's special brand of van der Waals scaling:

				We actually want to scale the vdW radii, rather than the vdW potential value,
				so let's breakdown the potential in terms of radii and other parameters:

				The following equation references are from the Amber manual v19.

				The 6-12 van der Waals potential has the general form given by Eqn 14.10:

				Vij = Aij/rij^12 - Bij/rij^6

				where rij is the radius (ie distance) between atoms i and j,
				and Aij and Bij are defined in Eqn 14.12:

				Aij = eij*Rmin^12
				Bij = 2*eij*Rmin^6

				where eij is defined in Eqn 14.14:

				eij = sqrt(eii*ejj)
				assuming eii = ei and ejj = ej, we can rewrite:
				eij = sqrt(ei*ej)

				and Rmin is defined in the text after Eqn 14.8:

				Rmin = Ri + Rj

				Here, ei and Ri are parameters defined for the Amber forcefields for
				atom i and we can assume they're given and constant.

				NOTE: Eqn 14.13 suggests a different mixing equation for Rmin:
				Rminij = 0.5*(Rmini + Rminj)
				assuming, Rminij = Rmin, Rmini = Ri, and Rminj = Rj, we can re-write:
				Rmin = 0.5*(Ri + Rj)

				After analying the outputs of LEaP, it appears that LEaP does NOT use Eqn 14.13 at all,
				not even when the two atom types are different.
				So the 0.5 scaling factor does not appear to be applied during LEaP's calculations.
				Dropping the 0.5 scaling factor allows us to produce values for Aij and Bij
				that are identical to LEaP's values.

				Ok, back to the derivation...

				To apply Osprey's version of vdW scaling, we need to multiply
				Ri and Rj by a factor of s:
				Let's massage Eqn 14.12 a bit to expose Ri and Rj:

				Aij = eij*Rmin^12
					= eij*( Ri + Rj )^12

				Now we apply the scaling factor s to Ri and Rj:

				   eij*( Ri*s + Rj*s )^12
				 = eij*[ s*(Ri + Rj) ]^12
				 = eij*(Ri + Rj)^12*s^12
				 = eij*Rmin^12*s^12
				 = Aij*s^12

				Therefore, Aij gets scaled by s^12,
				and similarly, Bij gets scaled by s^6.
			*/

			// apply the parametric vdW scaling
			val s2 = vdwScale*vdwScale
			val s6 = s2*s2*s2
			val s12 = s6*s6
			vdwA *= s12
			vdwB *= s6

			// 1-4 bonded atoms have scaled interactions
			// see: Amber manual 14.1.6, 1-4 Non-Bonded Interaction Scaling
			if (dist == 3) {
				esQ /= top.chargeDivisor14(ia, ib)
				top.vdwDivisor14(ia, ib).let {
					vdwA /= it
					vdwB /= it
				}
			}

			return AtomPairParams(esQ, vdwA, vdwB)
		}
	}

	override suspend fun parameterizeAtomPairs(infos: List<ForcefieldParams.MolInfo>) : AtomPairsParams {

		// make sure there are no duplicated molecules
		val molis = infos.map { it.moli }
		if (molis.toSet().size != molis.size) {
			throw IllegalArgumentException("duplicate molecules submitted")
		}

		// run amber to get the params
		val params = infos
			.map {
				it.atomsParams as AtomsParams
				AmberMolParams(
					it.mol,
					it.atomIndex,
					it.atomsParams.types,
					it.atomsParams.frcmods
				)
			}
			.calcParamsAmber()

		// NOTE: calcParamsAmber() eventually calls LEaP in a separate process
		// profiling shows this is by far the bottleneck in compiling atom pairs
		// maybe there's something we can do to speed this up so compiling goes faster?

		// read the amber forcefield params from the topology file
		val top = TopIO.read(params.top).mapTo(infos.map { it.mol })

		// build a map from the given indices to the topology indices
		val indexMaps = ArrayMap<ArrayMap<Int>>()
		for (i in 0 until top.numAtoms) {
			val atom = top.atom(i) ?: continue
			val info = infos
				.filter { atom in it.atomIndex }
				.also { if (it.size > 1) throw IllegalArgumentException("atom $atom is in multiple atom indices, can't determine molecule") }
				.firstOrNull()
				?: continue
			val atomi = info.atomIndex.getOrThrow(atom)
			indexMaps.getOrPut(info.moli) { ArrayMap() }[atomi] = i
		}

		return AtomPairsParams(top, indexMaps)
	}

	override fun combineAtomsParams(info1: ForcefieldParams.AtomsInfo, info2: ForcefieldParams.AtomsInfo): ForcefieldParams.AtomsParams {

		info1.atomsParams as AtomsParams
		info2.atomsParams as AtomsParams

		return AtomsParams(
			AmberTypes.combine(
				info1.atomsParams.types, info1.preferredAtomIndices, info1.ignoredAtomIndices,
				info2.atomsParams.types, info2.preferredAtomIndices, info2.ignoredAtomIndices
			),
			info1.atomsParams.frcmods + info2.atomsParams.frcmods
		)
	}
}
