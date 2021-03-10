package edu.duke.cs.osprey.molscope.gui.features.win

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.imgui.Imgui
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.WindowFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.view.ColorsMode


class MenuColorsMode : WindowFeature {

	override val id = FeatureId("colorsmode")

	override fun menu(imgui: Commands, win: WindowCommands) = imgui.run {
		if (beginMenu("Colors")) {
			for (mode in ColorsMode.values()) {
				if (selectable(mode.name, ColorsMode.current == mode)) {
					ColorsMode.current = mode
					Imgui.styleColors = when (mode) {
						ColorsMode.Dark -> Imgui.StyleColors.Dark
						ColorsMode.Light -> Imgui.StyleColors.Light
					}
				}

			}
			endMenu()
		}
	}
}
