package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.render.LoadedImage
import edu.duke.cs.osprey.gui.OspreyGui


enum class Icon(val id: String) {

	SignCheck("sign-check"),
	SignInfo("sign-info"),
	SignWarning("sign-warning"),
	SignError("sign-error");

	val resourcePath get() = "icons/$id.png"
}

object Icons {

	private val loaded = HashMap<Icon,LoadedImage>()

	fun get(win: WindowCommands, icon: Icon) =
		loaded.getOrPut(icon) {
			win.loadImage(OspreyGui.getResourceAsBytes(icon.resourcePath))
		}

	fun get(slidewin: SlideCommands, icon: Icon) =
		loaded.getOrPut(icon) {
			slidewin.loadImage(OspreyGui.getResourceAsBytes(icon.resourcePath))
		}
}

fun Commands.icon(win: WindowCommands, icon: Icon) {
	image(Icons.get(win, icon).descriptor)
}

fun Commands.icon(slidewin: SlideCommands, icon: Icon) {
	image(Icons.get(slidewin, icon).descriptor)
}
