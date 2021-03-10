package edu.duke.cs.osprey.molscope.gui

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.vulkan.ColorRGBA
import edu.duke.cs.osprey.molscope.view.ColorsMode


private var nextId = 0

/**
 * Renders a "i" icon which when hovered will show a popup.
 */
fun Commands.infoTip(block: () -> Unit) {

	pushStyleVar(Commands.StyleVar.WindowPadding, 5f, 2f)
	pushStyleVar(Commands.StyleVar.ChildRounding, 6f)
	pushStyleColor(Commands.StyleColor.Border, ColorRGBA.Int(128, 128, 128))
	pushStyleColor(Commands.StyleColor.Text, ColorRGBA.Int(128, 128, 128))
	pushStyleColor(Commands.StyleColor.ChildBg, when (ColorsMode.current) {
		ColorsMode.Dark -> ColorRGBA.Int(32, 32, 32)
		ColorsMode.Light -> ColorRGBA.Int(220, 220, 220)
	})

	child(
		"infoTip-${nextId++}",
		width = 17f,
		height = 18f,
		border = true,
		flags= IntFlags.of(Commands.BeginFlags.NoScrollbar)
	) {
		text("i")
	}

	popStyleVar(2)
	popStyleColor(3)

	if (isItemHovered()) {
		beginTooltip()
		block()
		endTooltip()
	}
}

/**
 * Calls `infoTip(block)` to display text in the tooltip.
 * All newlines in the text are discarded and the text is wrapped at the desired width.
 */
fun Commands.infoTip(text: String, wrapWidth: Float = 200f) = infoTip {
	pushTextWrapPos(wrapWidth)
	textWrapped(text.replace("\n"," "))
	popTextWrapPos()
}
