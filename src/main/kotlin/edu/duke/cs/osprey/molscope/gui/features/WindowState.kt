package edu.duke.cs.osprey.molscope.gui.features

import cuchaz.kludge.tools.Ref


class WindowState {

	val pOpen = Ref.of(false)
	var isOpen
		get() = pOpen.value
		set(value) { pOpen.value = value }

	private var wasOpen = false

	fun render(
		onOpen: () -> Unit = {},
		whenOpen: () -> Unit,
		onClose: () -> Unit = {}
	) {
		if (isOpen && !wasOpen) {
			wasOpen = true
			onOpen()
		}

		if (isOpen) {
			whenOpen()
		}

		if (!isOpen && wasOpen) {
			wasOpen = false
			onClose()
		}
	}
}