package edu.duke.cs.osprey.molscope.gui

import org.joml.Vector2f


/**
 * Tracks the mouse down and mouse up positions during a click or drag.
 * If the two positions match (within tolerance), a clicked returns true.
 */
class ClickTracker(val tolerance: Int = 5) {

	private var pos: Vector2f? = null

	fun clicked(slidewin: SlideCommands): Boolean {

		// track the down position
		if (slidewin.mouseLeftClick) {
			pos = Vector2f(slidewin.mouseOffset)
		}

		var isClick = false

		if (slidewin.mouseLeftRelease) {

			// register the click if the mouse up pos is near the mouse down pos
			pos?.let { pos ->
				if (slidewin.mouseOffset.distance(pos) <= tolerance) {
					isClick = true
				}
			}

			pos = null
		}

		return isClick
	}
}
