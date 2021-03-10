package edu.duke.cs.osprey.molscope.view

import cuchaz.kludge.vulkan.ColorRGBA
import java.util.*


enum class ColorsMode {

	Dark,
	Light;

	companion object {

		// start with dark mode by default (can change with window main menu)
		var current = Dark
	}
}

object ColorPalette {

	val darkGrey = Color(
		"Dark Grey",
		dark = ColorRGBA.Int(60, 60, 60),
		light = ColorRGBA.Int(90, 90, 90)
	)

	val lightGrey = Color(
		"Light Grey",
		dark = ColorRGBA.Int(200, 200, 200),
		light = ColorRGBA.Int(220, 220, 220)
	)

	val red = Color(
		"Red",
		dark = ColorRGBA.Int(237, 0, 0),
		light = ColorRGBA.Int(204, 0, 0)
	)

	val blue = Color(
		"Blue",
		dark = ColorRGBA.Int(0, 35, 254),
		light = ColorRGBA.Int(0, 51, 251)
	)

	val yellow = Color(
		"Yellow",
		dark = ColorRGBA.Int(255, 255, 0),
		light = ColorRGBA.Int(232, 205, 0)
	)

	val orange = Color(
		"Orange",
		dark = ColorRGBA.Int(255, 85, 0),
		light = ColorRGBA.Int(255, 102, 25)
	)
}


class Color(
	val name: String,
	val shades: Map<ColorsMode,ColorRGBA>
) {

	constructor(name: String, dark: ColorRGBA, light: ColorRGBA) : this(
		name,
		EnumMap<ColorsMode,ColorRGBA>(ColorsMode::class.java).apply {
			this[ColorsMode.Dark] = dark
			this[ColorsMode.Light] = light
		}
	)

	operator fun get(mode: ColorsMode): ColorRGBA =
		shades[mode] ?: throw NoSuchElementException("no color for mode: $mode")
}
