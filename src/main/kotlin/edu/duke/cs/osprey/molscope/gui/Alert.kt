package edu.duke.cs.osprey.molscope.gui

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags


class Alert {

	private val id = "Alert##${System.identityHashCode(this)}"

	private var text: String? = null

	fun show(text: String) {
		this.text = text
	}

	fun render(imgui: Commands) = imgui.run {

		if (text != null && !isPopupOpen(id)) {
			openPopup(id)
		}

		val text = text ?: return

		popupModal(id, null, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

			text(text)

			spacing()

			if (button("Ok")) {
				this@Alert.text = null
				closeCurrentPopup()
			}
		}
	}
}