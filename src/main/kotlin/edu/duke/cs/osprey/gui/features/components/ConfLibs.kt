package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.window.FileDialog
import cuchaz.kludge.window.FilterList
import edu.duke.cs.osprey.molscope.gui.Alert
import edu.duke.cs.osprey.gui.OspreyGui
import edu.duke.cs.osprey.gui.io.ConfLib
import edu.duke.cs.osprey.gui.io.UserSettings
import edu.duke.cs.osprey.gui.io.read
import edu.duke.cs.osprey.gui.prep.ConfSpace


/**
 * A cache for the built-in conflib metadata.
 */
object ConfLibs {

	private val builtInConflibPaths = listOf(
		"conflib/lovell.conflib"
	)

	data class ConfLibInfo(
		val path: String,
		val id: String,
		val name: String,
		val description: String?,
		val citation: String?
	) {
		fun read(): String =
			OspreyGui.getResourceAsString(path)

		fun load(): ConfLib =
			ConfLib.from(read())
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

	fun render(imgui: Commands) = imgui.run {

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
			for (conflib in confSpace.conflibs) {
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
					addLib(info.read())
				}
				conflibTooltip(null, info.description, info.citation)
			}
		}

		sameLine()

		if (button("Add from file")) {
			addLibFromFile()
		}

		alert.render(this)
	}

	private fun addLibFromFile() {
		FileDialog.openFiles(
			conflibFilter,
			defaultPath = UserSettings.openSaveDir
		)?.let { paths ->
			paths.firstOrNull()?.parent?.let { UserSettings.openSaveDir = it }
			for (path in paths) {
				addLib(path.read())
			}
		}
	}

	private fun addLib(toml: String) {

		val conflib = ConfLib.from(toml)

		if (confSpace.conflibs.contains(conflib)) {
			alert.show("Skipped adding duplicate Conformation Library:\n${conflib.name}")
			return
		}

		confSpace.conflibs.add(conflib)
		onAdd?.invoke(conflib)
	}
}
