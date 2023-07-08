package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.window.FileDialog
import cuchaz.kludge.window.FilterList
import edu.duke.cs.osprey.molscope.gui.Alert
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.*
import edu.duke.cs.osprey.gui.prep.ConfSpace
import edu.duke.cs.osprey.molscope.molecule.Molecule


/**
 * A cache for the built-in conflib metadata.
 */
object ConfLibs {

	private val builtInConflibPaths = listOf(
		"conflib/lovell.conflib",
		"conflib/D-lovell.conflib"
	)

	data class ConfLibInfo(
		val path: String,
		val id: String,
		val name: String,
		val description: String?,
		val citation: String?
	) {
		val lib : ConfLib by lazy {
			ConfLib.from(OspreyGui.getResourceAsString(path))
		}
	}

	val infos: List<ConfLibInfo> by lazy {
		ArrayList<ConfLibInfo>().apply {
			for (path in builtInConflibPaths) {
				val conflib = ConfLib.from(OspreyGui.getResourceAsString(path))
				add(ConfLibInfo(
					path,
					conflib.id,
					conflib.name,
					conflib.description,
					conflib.citation
				))
			}
		}
	}
}

class ConfLibPicker(val confSpace: ConfSpace) {

	private val conflibFilter = FilterList(listOf("conflib"))

	private val alert = Alert()

	var onAdd: ((ConfLib) -> Unit)? = null

	fun render(imgui: Commands, mol: Molecule) = imgui.run {

		fun conflibTooltip(name: String?, desc: String?, citation: String?) {
			if (isItemHovered()) {
				beginTooltip()
				if (name != null) {
					text(name)
				}
				if (desc != null) {
					text(desc)
				}
				if (citation != null) {
					text(citation)
				}
				endTooltip()
			}
		}

		// show available libraries
		text("Conformation Libraries")
		child("libs", 300f, 100f, true) {
			for (conflib in confSpace.getConflibsByMol(mol)) {
				text(conflib.name)
				conflibTooltip(conflib.name, conflib.description, conflib.citation)
			}
		}

		if (button("Add")) {
			openPopup("addlib")
		}
		popup("addlib") {
			for (info in ConfLibs.infos) {
				if (menuItem(info.name)) {
					addLib(mol, info.lib)
				}
				conflibTooltip(null, info.description, info.citation)
			}
		}

		sameLine()

		if (button("Add from file")) {
			addLibFromFile(mol)
		}

		alert.render(this)
	}

	private fun addLibFromFile(mol: Molecule) {
		FileDialog.openFiles(
			conflibFilter,
			defaultPath = UserSettings.openSaveDir
		)?.let { paths ->
			paths.firstOrNull()?.parent?.let { UserSettings.openSaveDir = it }
			for (path in paths) {
				addLib(mol, ConfLib.from(path.read()))
			}
		}
	}

	private fun addLib(mol: Molecule, conflib: ConfLib) {

		if (confSpace.getConflibsByMol(mol).contains(conflib)) {
			alert.show("Skipped adding duplicate Conformation Library:\n${conflib.name} for ${mol.name}")
			return
		}

		confSpace.addConflibByMol(mol, conflib)
		onAdd?.invoke(conflib)
	}
}
