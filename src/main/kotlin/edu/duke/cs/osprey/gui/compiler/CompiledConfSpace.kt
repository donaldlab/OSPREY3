package edu.duke.cs.osprey.gui.compiler

import org.joml.Vector3d


class CompiledConfSpace(
	val name: String,
	val forcefields: List<ForcefieldInfo>,
	val molInfos: List<MolInfo>,
	val resInfos: List<ResInfo>,
	val staticAtoms: List<AtomInfo>,
	/** in the same order as the forcefields */
	val staticEnergies: List<Double>,
	val positions: List<PosInfo>,
	/** in the same order as the forcefields */
	val atomPairs: List<AtomPairs>
) {

	data class ForcefieldInfo(
		val name: String,
		val ospreyImplementation: String,
		val settings: List<Pair<String,Any>>
	)

	data class MolInfo(
		val name: String,
		val type: String?
	) {
		val motions = ArrayList<MotionInfo>()
	}

	data class ResInfo(
		val chainId: String,
		val id: String,
		val type: String,
		val indexInChain: Int
	)

	data class PosInfo(
		val name: String,
		val wildType: String,
		val fragments: List<FragInfo>,
		val confs: List<ConfInfo>
	)

	data class FragInfo(
		val name: String,
		val atomNames: List<String>
	)

	data class ConfInfo(
		val id: String,
		val type: String,
		val atoms: List<AtomInfo>,
		val motions: List<MotionInfo>,
		/** to index into the atom pairs */
		val fragIndex: Int,
		/** in the same order as the forcefields */
		val internalEnergies: List<Double>
	)

	data class AtomInfo(
		val name: String,
		val pos: Vector3d,
		val molIndex: Int,
		val resIndex: Int
	)

	sealed class MotionInfo {

		data class DihedralAngle(
			val minDegrees: Double,
			val maxDegrees: Double,
			val abcd: List<Int>,
			val rotated: List<Int>
		) : MotionInfo()

		data class TranslationRotation(
			val maxTranslationDistance: Double,
			val maxRotationRadians: Double,
			val centroid: Vector3d
		) : MotionInfo()
	}
}
