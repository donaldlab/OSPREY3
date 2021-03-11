package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.window.FileDialog
import cuchaz.kludge.window.FilterList
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.combine
import edu.duke.cs.osprey.gui.io.*
import java.nio.file.Path


class ExportMol2(val getter: () -> Pair<List<Molecule>,String>) : SlideFeature {

	override val id = FeatureId("export.mol2")

	companion object {
		const val extension = "mol2"
	}

	val filterList = FilterList(listOf(extension))

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Export Mol2")) {
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

		// combine the mols and save the file
		val chainIdGenerator = ChainIdGeneratorAZ()
		val chainGenerator = ChainGeneratorSingleResidue(chainIdGenerator)
		getter()
			.let { (mols, name) ->
				mols.combine(name, chainIdGenerator, chainGenerator).first
			}
			.toMol2()
			.write(pathWithExt)

		// TODO: feedback to the user that the save worked?
	}
}
