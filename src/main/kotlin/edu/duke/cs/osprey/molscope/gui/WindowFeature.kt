package edu.duke.cs.osprey.molscope.gui

import cuchaz.kludge.imgui.Commands
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.features.HasFeatureId
import edu.duke.cs.osprey.molscope.render.LoadedImage


interface WindowFeature : HasFeatureId {

	/**
	 * Renders the menu item of your feature.
	 *
	 * Called from the window thread, not the creating thread,
	 * so be careful about synchronization of shared memory.
	 */
	fun menu(imgui: Commands, win: WindowCommands) {
		// do nothing by default
	}

	/**
	 * Renders the GUI of your feature.
	 *
	 * Called from the window thread, not the creating thread,
	 * so be careful about synchronization of shared memory.
	 */
	fun gui(imgui: Commands, win: WindowCommands) {
		// do nothing by default
	}
}

/**
 * Interface to talk to the window that is hosting this feature.
 *
 * Should only access this from the window thread.
 */
interface WindowCommands {

	fun showExceptions(block: () -> Unit)
	var shouldClose: Boolean
	fun addSlide(slide: Slide)
	fun removeSlide(slide: Slide): Boolean

	/**
	 * Loads an image for rendering.
	 * The image and all associated resources will be automatically
	 * cleaned up when the window is cleaned up, so callers will
	 * not need to worry about cleaning up images loaded with this method.
	 */
	fun loadImage(bytes: ByteArray): LoadedImage
}
