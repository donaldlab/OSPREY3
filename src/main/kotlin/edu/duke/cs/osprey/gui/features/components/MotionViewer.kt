package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import edu.duke.cs.osprey.molscope.molecule.Atom
import kotlin.random.Random


interface MotionViewer {

	val label: String
	fun gui(imgui: Commands)
	fun jiggle(rand: Random)
	fun reset()
	fun mapAtomToOriginalMol(atom: Atom): Atom
}
