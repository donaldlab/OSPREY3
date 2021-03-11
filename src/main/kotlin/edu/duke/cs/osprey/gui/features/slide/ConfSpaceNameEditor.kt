package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.gui.prep.ConfSpace


class ConfSpaceNameEditor(val confSpace: ConfSpace) : SlideFeature {

	override val id = FeatureId("edit.confSpaceName")

	private val winState = WindowState()

	private val buf = Commands.TextBuffer(1024)

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Name")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		winState.render(
			onOpen = {

				// init the text buffer
				buf.text = confSpace.name
			},
			whenOpen = {

				// draw the window
				window("Name Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					if (inputText("Name", buf)) {
						confSpace.name = buf.text
					}
				}
			},
			onClose = {

				// update the slide name too
				slide.name = buf.text
			}
		)
	}
}