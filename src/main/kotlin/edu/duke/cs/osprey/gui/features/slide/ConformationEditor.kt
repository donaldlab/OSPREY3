package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.features.components.ConfPosEditor
import edu.duke.cs.osprey.gui.features.components.MolMotionEditor
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.motions.MolMotion
import edu.duke.cs.osprey.gui.motions.TranslationRotation
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.prep.DesignPosition
import java.math.BigInteger
import java.text.NumberFormat
import kotlin.collections.ArrayList


class ConformationEditor(val confSpace: ConfSpace) : SlideFeature {

	override val id = FeatureId("edit.conformations")

	private var posEditor: ConfPosEditor? = null
	private var motionEditor: MolMotionEditor? = null
	private val winState = WindowState()

	class MotionInfo(val index: Int, val desc: MolMotion.Description) {

		private fun describe() =
			when (desc) {
				is DihedralAngle.MolDescription -> {
					val atomNames = listOf(desc.a, desc.b, desc.c, desc.d)
						.joinToString(", ") { it.name }
					"Dihedral Angle @ ${desc.mol.name}: $atomNames"
				}
				is TranslationRotation.MolDescription -> {
					"Translation/Rotation @ ${desc.mol.name}"
				}
				else -> "(Unrecognized Motion)"
			}

		val label = "${describe()}##$index"
	}

	inner class MolInfo(val mol: Molecule, val molType: MoleculeType) {

		val label = mol.toString()
		val posInfos = ArrayList<PosInfo>()
		var numConfs: BigInteger = BigInteger.ZERO

		val motions get() =
			confSpace.molMotions
				.getOrPut(mol) { ArrayList() }

		val motionInfos = motions
			.mapIndexed { i, motion -> MotionInfo(i, motion) }
			.toMutableList()

		fun updateCounts() {
			numConfs = BigInteger.ONE
			for (posInfo in posInfos) {
				numConfs *= posInfo.posConfSpace.confs.size.toBigInteger()
			}
		}

		fun makeNewPosition(): PosInfo {

			val positions = confSpace.designPositionsByMol
				.getOrPut(mol) { ArrayList() }

			// choose a default but unique name
			val prefix = "Pos "
			val maxNum = positions
				.mapNotNull {
					it.name
						.takeIf { it.startsWith(prefix) }
						?.substring(prefix.length)
						?.toIntOrNull()
				}
				.max()
				?: 0
			val num = maxNum + 1

			// create the position and add it
			val pos = DesignPosition("$prefix$num", "none", mol)
			positions.add(pos)

			val posInfo = PosInfo(pos, posInfos.size)
			posInfos.add(posInfo)
			return posInfo
		}

		fun removePosition(posInfo: PosInfo) {

			posInfos.remove(posInfo)
			confSpace.designPositionsByMol[mol]?.remove(posInfo.pos)
			confSpace.positionConfSpaces.remove(posInfo.pos)
		}
	}

	inner class PosInfo(val pos: DesignPosition, val index: Int) {

		val label = "${pos.name}###pos$index"
		val posConfSpace get() = confSpace.positionConfSpaces.getOrMake(pos)
		val isMutable get() = posConfSpace.isMutable()
		val isFlexible get() = !isMutable
	}

	private val molInfos = ArrayList<MolInfo>()

	private fun resetInfos() {

		molInfos.clear()

		for ((moltype, mol) in confSpace.mols) {
			molInfos.add(MolInfo(mol, moltype).apply {
				confSpace.designPositionsByMol[mol]?.forEachIndexed { index, pos ->
					posInfos.add(PosInfo(pos, index))
				}
			})
		}
	}

	private var numConfs: BigInteger = BigInteger.ZERO
	private val confFormatter = NumberFormat.getIntegerInstance()
		.apply {
			isGroupingUsed = true
		}

	private var selectedPosInfo = null as PosInfo?
	private var selectedMotionInfo = null as MotionInfo?

	private fun updateCounts() {
		if (molInfos.isEmpty()) {
			numConfs = BigInteger.ZERO
		} else {
			numConfs = BigInteger.ONE
			for (molInfo in molInfos) {
				molInfo.updateCounts()
				numConfs *= molInfo.numConfs
			}
		}
	}

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Flexibility")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		winState.render(
			onOpen = {
				resetInfos()
				updateCounts()
			},
			whenOpen = {

				// draw the window
				window("Flexibility Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					var resetInfos = false

					tabBar("tabs") {

						for (molInfo in molInfos) {
							tabItem(molInfo.label) {

								withId("mutable") {

									text("Mutable Positions:")
									spacing()
									columns(2)
									indent(20f)
									val mutInfos = molInfo.posInfos.filter { it.isMutable }
									for (posInfo in mutInfos) {

										if (selectable(posInfo.pos.name, selectedPosInfo === posInfo)) {
											selectedPosInfo = posInfo
										}
										nextColumn()
										text("${posInfo.posConfSpace.confs.size} confs")
										nextColumn()
									}
									columns(1)
									if (mutInfos.isEmpty()) {
										text("(no mutable positions)")
									}
									unindent(20f)
								}

								spacing()

								withId("flexible") {

									text("Flexible Positions:")
									child("flexpos", 360f, 200f, true) {
										columns(2)
										for (posInfo in molInfo.posInfos.filter { it.isFlexible }) {

											if (selectable(posInfo.label, selectedPosInfo === posInfo)) {
												selectedPosInfo = posInfo
											}
											nextColumn()
											text("${posInfo.posConfSpace.confs.size} confs")
											nextColumn()
										}
										columns(1)
									}

									if (button("Add")) {

										// start the conformation editor
										posEditor = ConfPosEditor(
											confSpace,
											molInfo,
											molInfo.makeNewPosition(),
											onClose = {
												resetInfos()
												updateCounts()
												posEditor = null
											}
										)
									}

									sameLine()

									val selectedPosInfo = selectedPosInfo

									enabledIf(selectedPosInfo != null) {
										if (button("Edit")) {

											// sadly the compiler isn't quite smart enough to figure out this can't be null
											// so put in a runtime check to make flow typing work correctly
											selectedPosInfo!!

											// start the conformation editor
											posEditor = ConfPosEditor(
												confSpace,
												molInfo,
												selectedPosInfo,
												onClose = {
													resetInfos()
													updateCounts()
													posEditor = null
												}
											)
										}
									}

									sameLine()

									// allow removing flexible positions only
									enabledIf(selectedPosInfo != null && selectedPosInfo.isFlexible) {
										if (button("Remove")) {

											// sadly the compiler isn't quite smart enough to figure out this can't be null
											// so put in a runtime check just to make the compiler happy
											selectedPosInfo!!

											// remove the design position
											molInfo.removePosition(selectedPosInfo)
											updateCounts()
										}
									}
								}

								spacing()

								text("Conformations: ${confFormatter.format(molInfo.numConfs)}")

								spacing()
								spacing()
								spacing()

								withId("molMotions") {

									text("Molecule Motions:")
									child("motions", 360f, 200f, true) {
										for (motionInfo in molInfo.motionInfos) {

											if (selectable(motionInfo.label, selectedMotionInfo === motionInfo)) {
												selectedMotionInfo = motionInfo
											}
										}
									}

									// add button
									if (button("Add")) {

										// start the motion editor
										motionEditor = MolMotionEditor(
											molInfo,
											null,
											onClose = {
												motionEditor = null
												resetInfos()
												updateCounts()
											}
										)
									}

									sameLine()

									val selectedMotionInfo = selectedMotionInfo

									// edit button
									enabledIf(selectedMotionInfo != null) {
										if (button("Edit")) {

											// sadly the compiler isn't quite smart enough to figure out this can't be null
											// so put in a runtime check to make flow typing work correctly
											selectedMotionInfo!!

											// start the motion editor
											motionEditor = MolMotionEditor(
												molInfo,
												selectedMotionInfo.desc,
												onClose = {
													motionEditor = null
													resetInfos()
													updateCounts()
												}
											)
										}
									}

									sameLine()

									// remove button
									enabledIf(selectedMotionInfo != null) {
										if (button("Remove")) {

											// sadly the compiler isn't quite smart enough to figure out this can't be null
											// so put in a runtime check to make flow typing work correctly
											selectedMotionInfo!!

											molInfo.motions.remove(selectedMotionInfo.desc)
											this@ConformationEditor.selectedMotionInfo = null

											resetInfos = true
										}
									}
								}
							}
						}
					}

					// for multiple molecules, show the combined conformation count
					if (molInfos.size > 1) {
						spacing()
						text("Combined Conformations: ${confFormatter.format(numConfs)}")
					}

					// render the editors, when active
					posEditor?.gui(imgui, slide, slidewin)
					motionEditor?.gui(imgui, slide, slidewin)

					if (resetInfos) {
						resetInfos()
						updateCounts()
					}
				}
			},
			onClose = {

				// cleanup infos
				molInfos.clear()
			}
		)
	}
}
