package edu.duke.cs.osprey.gui.features.win

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.window.FileDialog
import cuchaz.kludge.window.FilterList
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.WindowFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.gui.io.*
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.gui.prep.ConfSpacePrep
import java.nio.file.Path


class OpenConfSpace : WindowFeature {

	override val id = FeatureId("open.confspace")

	val filterList = FilterList(listOf("confspace"))

	override fun menu(imgui: Commands, win: WindowCommands) = imgui.run {
		if (menuItem("Open Conformation Space")) {
			FileDialog.openFile(filterList, UserSettings.openSaveDir)?.let { path ->
				UserSettings.openSaveDir = path.parent
				open(win, path)
			}
		}
	}

	private fun open(win: WindowCommands, path: Path) = win.showExceptions {

		// read the source file
		val toml = path.read()

		// resume a previous prep
		ConfSpacePrep(win, ConfSpace.fromToml(toml))
	}
}
