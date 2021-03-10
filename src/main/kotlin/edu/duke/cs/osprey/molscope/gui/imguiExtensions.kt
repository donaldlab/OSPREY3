package edu.duke.cs.osprey.molscope.gui

import cuchaz.kludge.imgui.Commands
import edu.duke.cs.osprey.molscope.render.LoadedImage


/**
 * Calls `pushStyleVar` to temporarily modify the ImGUI style,
 * and `pushItemFlag` to disable interactions with controls.
 * Make sure to nest other calls to push/pop StyleVar  and push/pop ItemFlag correctly.
 *
 * ImGUI doesn't have an (official) enabled/disabled system for most controls (yet),
 * so this uses internal/BETA functionality for now.
 * See: https://github.com/ocornut/imgui/issues/211
 */
fun Commands.pushDisabled() {
	pushStyleVar(Commands.StyleVar.Alpha, 0.4f)
	pushItemFlag(Commands.ItemFlags.Disabled, true)
}

fun Commands.popDisabled(num: Int = 1) {
	popItemFlag()
	popStyleVar(num)
}

inline fun <R> Commands.disabledIf(isDisabled: Boolean, block: () -> R): R {
	if (isDisabled) {
		pushDisabled()
	}
	val ret = block()
	if (isDisabled) {
		popDisabled()
	}
	return ret
}

inline fun <R> Commands.enabledIf(isEnabled: Boolean, block: () -> R) =
	disabledIf(!isEnabled, block)


interface WithColumns {
	fun column(width: Float? = null, block: () -> Unit)
}
fun Commands.columns(num: Int, border: Boolean = false, block: WithColumns.() -> Unit) {
	columns(num, border = border)
	try {
		var offset = 0f
		var i = 0
		object : WithColumns {
			override fun column(width: Float?, block: () -> Unit) {
				child("column$i", block = block)
				i += 1
				if (width != null) {
					offset += width
					setColumnOffset(i, offset)
				}
				nextColumn()
			}
		}.block()
	} finally {
		columns(1)
	}
}

fun Commands.image(img: LoadedImage) {
	image(img.descriptor, img.width.toFloat(), img.height.toFloat())
}

// TODO: move this into Kludge?
inline fun Commands.listBox(label: String, itemsCount: Int, heightInItems: Int = itemsCount, block: () -> Unit) {
	if (listBoxHeader(label, itemsCount, heightInItems)) {
		try {
			block()
		} finally {
			listBoxFooter()
		}
	}
}
