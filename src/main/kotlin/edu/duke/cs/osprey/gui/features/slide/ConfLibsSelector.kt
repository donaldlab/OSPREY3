package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import edu.duke.cs.osprey.gui.features.components.ConfLibPicker
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState


class ConfLibsSelector(val confSpace: ConfSpace) : SlideFeature {

	override val id = FeatureId("edit.confLibs")

	private val winState = WindowState()

	private val conflibPicker = ConfLibPicker(confSpace)

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Conformation Libraries")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		winState.render(
			whenOpen = {

				// draw the window
				window("Conformation Libraries##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					tabBar("tabs") {
						for (mol in confSpace.mols) {
							tabItem(mol.second.name, flags = IntFlags.of(Commands.TabItemFlags.None)) {
								conflibPicker.render(imgui, mol.second)
							}
						}
					}
				}
			}
		)
	}
}
