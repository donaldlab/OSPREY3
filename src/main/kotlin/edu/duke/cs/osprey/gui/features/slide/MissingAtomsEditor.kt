package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.tools.toFloat
import edu.duke.cs.osprey.gui.forcefield.amber.MissingAtom
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.gui.infoTip
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.forcefield.amber.inferMissingAtomsAmberBlocking


class MissingAtomsEditor : SlideFeature {

	override val id = FeatureId("edit.missingAtoms")

	private val winState = WindowState()

	private data class AddedAtom(
		val view: MoleculeRenderView,
		val missingAtom: MissingAtom
	) {

		val pIncluded = Ref.of(true)
	}

	/** keep track of the atoms we added so we can remove them later */
	private var addedAtoms: MutableList<AddedAtom>? = null

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Missing Atoms")) {
			winState.isOpen = true
		}
	}

	private fun Slide.Locked.molViews() = views.mapNotNull { it as? MoleculeRenderView }

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val molViews = slide.molViews()

		winState.render(
			onOpen = {
				// clear any previous state
				addedAtoms = null
			},
			whenOpen = {

				// draw the window
				window("Missing Atom Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					text("Tools:")
					indent(10f)

					if (button("Add missing atoms automatically")) {
						slidewin.showExceptions {
							addedAtoms = guessAtoms(molViews)
						}
					}
					sameLine()
					infoTip("""
						|This tool infers the positions of missing heavy atoms for molecules
						|with regular structure, such as proteins and nucleic acids,
						|based on definitions in the AmberTools package.
					""".trimMargin())

					unindent(10f)

					addedAtoms?.let { addedAtoms ->

						text("Added Atoms:")
						child("addedAtoms", 300f, 200f, true) {

							if (addedAtoms.isNotEmpty()) {

								var atomi = 0

								// show a checkbox to toggle the atom on/off
								for (addedAtom in addedAtoms) {

									val label = "${addedAtom.missingAtom}##$atomi"
									if (checkbox(label, addedAtom.pIncluded)) {
										updateAtom(molViews, addedAtom)
									}

									// show a context menu for the checkbox
									popupContextItem("addedAtom:$atomi") {

										if (button("Center Camera")) {
											slidewin.camera.lookAt(addedAtom.missingAtom.atom.pos.toFloat(), slide.views)
											slidewin.camera.changed()
											closeCurrentPopup()
										}
									}
								}

								atomi += 1

							} else {

								text("no missing atoms were added")
							}
						}

						if (addedAtoms.isNotEmpty()) {
							if (button("Clear all automatically-added atoms")) {
								clearAtoms(addedAtoms)
								this@MissingAtomsEditor.addedAtoms = null
							}
						}
					}
				}
			},
			onClose = {
				// clear any previous state
				addedAtoms = null
			}
		)
	}

	private fun guessAtoms(views: List<MoleculeRenderView>) = ArrayList<AddedAtom>().apply {

		for (view in views) {

			val mol = view.molStack.originalMol
			val missingAtoms = mol.inferMissingAtomsAmberBlocking()
			for (missingAtom in missingAtoms) {

				// add the atoms to the molecule
				missingAtom.add()

				// add them to the internal state too, so we can see what was added
				add(AddedAtom(view, missingAtom))
			}

			view.moleculeChanged()
		}
	}

	private fun updateAtom(views: List<MoleculeRenderView>, addedAtom: AddedAtom) {

		val view = views.find { it.molStack.originalMol == addedAtom.missingAtom.mol } ?: return

		// add or remove the atom
		if (addedAtom.pIncluded.value) {
			addedAtom.missingAtom.add()
		} else {
			addedAtom.missingAtom.remove()
		}

		view.moleculeChanged()
	}

	private fun clearAtoms(addedAtoms: MutableList<AddedAtom>) {

		for (addedAtom in addedAtoms) {
			addedAtom.missingAtom.remove()
		}

		// update the affected views
		addedAtoms
			.map { it.view }
			.toSet()
			.forEach { it.moleculeChanged() }
	}
}
