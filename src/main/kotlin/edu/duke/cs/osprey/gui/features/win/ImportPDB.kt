package edu.duke.cs.osprey.gui.features.win

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.window.FileDialog
import cuchaz.kludge.window.FilterList
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.WindowFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.gui.io.UserSettings
import edu.duke.cs.osprey.gui.io.fromPDB
import edu.duke.cs.osprey.gui.io.read
import edu.duke.cs.osprey.gui.prep.MoleculePrep
import java.nio.file.Path


class ImportPDB : WindowFeature {

	override val id = FeatureId("import.pdb")

	val filterList = FilterList(listOf("pdb"))
	// TEMP: pdb.gz?

	override fun menu(imgui: Commands, win: WindowCommands) = imgui.run {
		if (menuItem("Import PDB")) {
			FileDialog.openFile(filterList, UserSettings.openSaveDir)?.let { path ->
				UserSettings.openSaveDir = path.parent
				open(win, path)
			}
		}
	}

	private fun open(win: WindowCommands, path: Path) = win.showExceptions {

		// start a new prep
		MoleculePrep(win, listOf(Molecule.fromPDB(path.read())))
	}
}
