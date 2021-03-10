package edu.duke.cs.osprey.molscope.gui

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import java.awt.Toolkit
import java.awt.datatransfer.StringSelection
import java.io.PrintWriter
import java.io.StringWriter


class ExceptionViewer {

	private val exceptions = ArrayList<Throwable>()
	private var textBuffer: Commands.TextBuffer? = null

	private val pOpen = Ref.of(false)

	fun add(t: Throwable) {
		exceptions.add(t)
		updateText()
		pOpen.value = true
	}

	fun clear() {
		exceptions.clear()
		updateText()
		pOpen.value = false
	}

	private fun updateText() {

		val sw = StringWriter()
		val pw = PrintWriter(sw)
		exceptions.forEachIndexed { i, t ->
			if (i > 0) {
				pw.print("\n")
			}
			t.printStackTrace(pw)
		}

		textBuffer = Commands.TextBuffer.of(sw.toString())
	}

	fun gui(imgui: Commands) = imgui.run {
		if (pOpen.value) {
			window("Exceptions", pOpen, IntFlags.of(Commands.BeginFlags.NoResize)) {

				val buf = textBuffer
				if (buf != null) {
					inputTextMultiline(
						"",
						buf,
						600f, 400f,
						IntFlags.of(Commands.InputTextFlags.ReadOnly)
					)
					// yes, the horizontal scroll bar doesn't appear automatically for multiline input text widgets
					// this is an ongoing enhancement in ImGUI, but I wouldn't wait for it
				} else {
					text("No exceptions")
				}

				popupContextItem("exception-text") {
					if (button("Copy all")) {
						copyToClipboard()
						closeCurrentPopup()
					}
				}

				if (button("Copy to Clipboard")) {
					copyToClipboard()
				}

				sameLine()

				if (button("Clear")) {
					clear()
				}
			}
		}
	}

	private fun copyToClipboard() {
		val buf = textBuffer ?: return
		Toolkit.getDefaultToolkit()
			.systemClipboard
			.setContents(StringSelection(buf.toString()), null)
	}
}
