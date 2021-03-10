package edu.duke.cs.osprey.molscope.gui.features.win

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.WindowFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId


class DevImguiDemo : WindowFeature {

	override val id = FeatureId("dev.imguidemo")

	val pOpen = Ref.of(false)

	override fun menu(imgui: Commands, win: WindowCommands) = imgui.run {
		if (menuItem("ImGUI Demo")) {
			pOpen.value = true
		}
	}

	override fun gui(imgui: Commands, win: WindowCommands) = imgui.run {
		if (pOpen.value) {
			showDemoWindow(pOpen)
		}
	}
}
