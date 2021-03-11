package edu.duke.cs.osprey.gui.prep

import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.compiler.NetCharges
import edu.duke.cs.osprey.gui.forcefield.*
import kotlinx.coroutines.runBlocking


class ForcefieldAnalyzer(val mols: List<Molecule>) {

	val forcefields = ForcefieldSet()
	val netCharges = NetCharges()

	interface Parameterized {
		val forcefields: List<ForcefieldParams>
		fun calcEnergy(ff: Forcefield): Double
		// TODO: add functions to return every atom pair?
	}

	fun parameterize(): Parameterized {

		val ffinfos = HashMap<Forcefield,ForcefieldInfo>()
		val atomIndices = mols.map { AtomIndex(it.atoms) }

		runBlocking {

			for (ff in forcefields) {

				// get the atoms params
				val atomsParams = mols.mapIndexed { moli, mol ->
					ff.parameterizeAtoms(
						mol,
						atomIndices[moli],
						netCharges[mol]?.netChargeOrThrow
					)
				}

				val molInfos = mols.mapIndexed { moli, mol ->
					ForcefieldParams.MolInfo(moli, mol, atomsParams[moli], atomIndices[moli])
				}

				// get the atom pairs params
				val atomPairsParams = ff.parameterizeAtomPairs(molInfos)

				ffinfos[ff.forcefield] = ForcefieldInfo(ff, atomsParams, molInfos, atomPairsParams)
			}
		}

		return ParameterizedImpl(ffinfos)
	}

	private class ForcefieldInfo(
		val ff: ForcefieldParams,
		val atomsParams: List<ForcefieldParams.AtomsParams>,
		val molInfos: List<ForcefieldParams.MolInfo>,
		val atomPairsParams: ForcefieldParams.AtomPairsParams
	)

	private inner class ParameterizedImpl(val ffinfos: Map<Forcefield,ForcefieldInfo>) : Parameterized {

		override val forcefields = this@ForcefieldAnalyzer.forcefields

		private fun ffinfoOrThrow(ff: Forcefield): ForcefieldInfo =
			ffinfos[ff] ?: throw NoSuchElementException()

		override fun calcEnergy(ff: Forcefield): Double {
			val ffinfo = ffinfoOrThrow(ff)
			return ForcefieldCalculator.calc(
				ffinfo.atomPairsParams,
				ffinfo.molInfos.map {
					ForcefieldCalculator.MolInfo(it.moli, it.mol, it.mol.atoms, it.atomIndex, ffinfo.atomsParams[it.moli])
				},
				ffinfo.ff.unconnectedDistance
			)
		}
	}
}
