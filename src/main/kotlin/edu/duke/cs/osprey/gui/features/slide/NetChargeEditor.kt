package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.gui.forcefield.amber.MoleculeType
import edu.duke.cs.osprey.gui.prep.MoleculePrep
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.columns
import edu.duke.cs.osprey.molscope.gui.enabledIf
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView


class NetChargeEditor(val prep: MoleculePrep) : SlideFeature {

	override val id = FeatureId("edit.netCharges")

	private val winState = WindowState()

	inner class MolInfo(val mol: Molecule, val type: MoleculeType) {

		val active = Ref.of(mol.netCharge != null)
		val value: Ref<Int> = Ref.of(
			getter = { mol.netCharge ?: 0 },
			setter = { mol.netCharge = it }
		)

		 fun reset() {
			 if (active.value) {
				 mol.netCharge = 0
			 } else {
				 mol.netCharge = null
			 }
		 }
	}

	private val molInfos = ArrayList<MolInfo>()

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Net Charges")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		winState.render(
			onOpen = {

				// init the molecule infos
				molInfos.clear()
				for ((type, mol) in prep.getIncludedTypedMols()) {
					molInfos.add(MolInfo(mol, type))
				}
			},
			whenOpen = {

				// draw the window
				setNextWindowSize(400f, 0f)
				window("Net Charge Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					columns(2)

					for ((infoi, info) in molInfos.withIndex()) {

						spacing()
						spacing()
						spacing()

						text("${info.mol}")
						nextColumn()

						spacing()
						spacing()
						spacing()

						when (info.type) {

							MoleculeType.SmallMolecule -> {
								if (checkbox("##active_$infoi", info.active)) {
									info.reset()
								}
								sameLine()
								enabledIf(info.active.value) {
									inputInt("##value_$infoi", info.value)
								}
							}

							else -> {
								text("(net charge not needed)")
							}
						}
						nextColumn()
					}

					columns(1)
				}
			},
			onClose = {

				// cleanup the molecule infos
				molInfos.clear()
			}
		)
	}
}
