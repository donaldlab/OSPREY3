package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.tools.toFloat
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.gui.infoTip
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.forcefield.amber.inferMissingAtomsAmber
import kotlinx.coroutines.runBlocking

class MissingAtomsEditor : SlideFeature {

	override val id = FeatureId("edit.missingAtoms")

	private val winState = WindowState()

	private data class AddedAtom(
		val mol: Molecule,
		val res: Polymer.Residue?,
		val atom: Atom,
		val location: String?
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

								// show a checkbox to toggle the atom on/off
								for (addedAtom in addedAtoms) {

									val location = addedAtom.location
									val label = if (location != null) {
										"${addedAtom.atom.name} @ $location"
									} else {
										addedAtom.atom.name
									}

									if (checkbox(label, addedAtom.pIncluded)) {
										updateAtom(molViews, addedAtom)
									}

									// show a context menu for the checkbox
									popupContextItem("addedAtom:$label") {

										if (button("Center Camera")) {
											slidewin.camera.lookAt(addedAtom.atom.pos.toFloat(), slide.views)
											slidewin.camera.changed()
											closeCurrentPopup()
										}
									}
								}

							} else {

								text("no missing atoms were added")
							}
						}

						if (addedAtoms.isNotEmpty()) {
							if (button("Clear all automatically-added atoms")) {
								clearAtoms(addedAtoms, molViews)
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
			val missingAtoms = runBlocking { mol.inferMissingAtomsAmber() }
			for ((atom, res) in missingAtoms) {

				// add the atoms to the molecule
				mol.atoms.add(atom)
				res?.atoms?.add(atom)

				// make a friendly description for the location
				val resId = res?.id
				val chainId = (mol as Polymer).chains.find { res in it.residues }?.id
				val location = if (resId != null && chainId != null) {
					"$chainId$resId"
				} else {
					null
				}

				// add them to the internal state too, so we can see what was added
				add(AddedAtom(mol, res, atom, location))
			}

			view.moleculeChanged()
		}
	}

	private fun updateAtom(views: List<MoleculeRenderView>, addedAtom: AddedAtom) {

		val view = views.find { it.molStack.originalMol == addedAtom.mol } ?: return

		// add or remove the atom
		addedAtom.run {
			if (pIncluded.value) {
				mol.atoms.add(atom)
				res?.atoms?.add(atom)
			} else {
				mol.atoms.remove(atom)
				res?.atoms?.remove(atom)
			}
		}

		view.moleculeChanged()
	}

	private fun clearAtoms(addedAtoms: MutableList<AddedAtom>, views: List<MoleculeRenderView>) {

		for (addedAtom in addedAtoms) {
			addedAtom.run {
				mol.atoms.remove(atom)
				res?.atoms?.remove(atom)
			}
		}

		// update the affected views
		addedAtoms
			.map { it.mol }
			.toSet()
			.mapNotNull { mol -> views.find { it.molStack.originalMol == mol } }
			.forEach { it.moleculeChanged() }
	}
}
