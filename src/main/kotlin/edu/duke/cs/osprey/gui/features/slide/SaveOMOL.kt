package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.window.FileDialog
import cuchaz.kludge.window.FilterList
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.io.UserSettings
import edu.duke.cs.osprey.gui.io.toOMOL
import edu.duke.cs.osprey.gui.io.write
import java.nio.file.Path


class SaveOMOL(val getter: () -> List<Molecule>) : SlideFeature {

	override val id = FeatureId("save.omol")

	companion object {
		const val extension = "omol"
	}

	val filterList = FilterList(listOf(extension))

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Save OMOL")) {
			FileDialog.saveFile(filterList, UserSettings.openSaveDir)?.let { path ->
				UserSettings.openSaveDir = path.parent
				save(slide, slidewin, path)
			}
		}
	}

	private fun save(slide: Slide.Locked, slidewin: SlideCommands, path: Path) = slidewin.showExceptions {

		// append the file extension if needed
		var filename = path.fileName.toString()
		if (!filename.endsWith(".$extension")) {
			filename += ".$extension"
		}
		val pathWithExt = path.parent.resolve(filename)

		// get the molecules and save them to the file
		getter()
			.toOMOL()
			.write(pathWithExt)

		// TODO: feedback to the user that the save worked?
	}
}
