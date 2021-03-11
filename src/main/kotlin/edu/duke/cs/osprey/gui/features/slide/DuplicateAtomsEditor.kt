package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.tools.toFloat
import cuchaz.kludge.tools.toString
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


class DuplicateAtomsEditor : SlideFeature {

	override val id = FeatureId("edit.duplicateAtoms")

	private val winState = WindowState()

	private data class DupeAtoms(
		val view: MoleculeRenderView,
		val mol: Molecule,
		val res: Polymer.Residue?,
		val atoms: List<Atom>,
		val location: String?
	) {

		val pIncluded = atoms.map { Ref.of(true) }

		val name get() = atoms.first().name

		fun update(atomi: Int) {

			// add or remove the atom
			val atom = atoms[atomi]
			val pIncluded = pIncluded[atomi]
			if (pIncluded.value) {
				mol.atoms.add(atom)
				res?.atoms?.add(atom)
			} else {
				mol.atoms.remove(atom)
				res?.atoms?.remove(atom)
			}

			view.moleculeChanged()
		}

		fun resolve() =
			pIncluded.count { it.value } == 1
	}

	private val dupeAtoms = ArrayList<DupeAtoms>()
	private val toResolve = ArrayList<Int>()

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Duplicated Atoms")) {
			winState.isOpen = true
		}
	}

	private fun Slide.Locked.molViews() = views.mapNotNull { it as? MoleculeRenderView }

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val molViews = slide.molViews()

		winState.render(
			onOpen = {
				detectDupeAtoms(molViews)
			},
			whenOpen = {

				toResolve.clear()

				// draw the window
				window("Duplicated Atom Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					if (dupeAtoms.isNotEmpty()) {
						for ((dupei, dupeAtoms) in dupeAtoms.withIndex()) {

							spacing()
							spacing()
							spacing()

							text("Duplicated: ${dupeAtoms.name} @ ${dupeAtoms.location}:")
							indent(20f)

							// show checkboxes to toggle the atoms on/off
							for ((atomi, pair) in dupeAtoms.atoms.zip(dupeAtoms.pIncluded).withIndex()) {
								val (atom, pIncluded) = pair

								val label = "${atom.element.symbol} @ ${atom.pos.toString(6)}##$atomi"
								if (checkbox(label, pIncluded)) {
									dupeAtoms.update(atomi)
								}

								// show a context menu for the checkbox
								popupContextItem("dupeAtoms:$label") {

									if (button("Center Camera")) {
										slidewin.camera.lookAt(atom.pos.toFloat(), slide.views)
										slidewin.camera.changed()
										closeCurrentPopup()
									}
								}
							}

							if (button("Resolve##$dupei")) {
								toResolve.add(dupei)
							}
							sameLine()
							infoTip("""
								|If the atom is still duplicated, this button does nothing.
								|Otherwise, this button will remove this entry from the duplicated atoms list.
							""".trimMargin())

							unindent(20f)
						}
					} else {
						text("(No duplicate atoms detected)")
					}
				}

				// handle the resolutions
				toResolve
					.map { dupeAtoms[it] }
					.forEach {
						if (it.resolve()) {
							dupeAtoms.remove(it)
						}
					}
			},
			onClose = {
				// clear any previous state
				dupeAtoms.clear()
			}
		)
	}

	private fun detectDupeAtoms(views: List<MoleculeRenderView>) {

		dupeAtoms.clear()

		for (view in views) {
			val mol = view.currentMol

			// look for duplicated atoms
			when (mol) {

				// check atoms for each residue
				is Polymer -> {
					for (chain in mol.chains) {
						for (res in chain.residues) {
							detectDupeAtoms(view, mol, chain, res, res.atoms)
						}
					}
				}

				// check all the atoms at once
				else -> {
					detectDupeAtoms(view, mol, null, null, mol.atoms)
				}
			}
		}
	}

	private fun detectDupeAtoms(view: MoleculeRenderView, mol: Molecule, chain: Polymer.Chain?, res: Polymer.Residue?, atoms: List<Atom>) {

		// look for duplicates by atom name
		val dupesByName = HashMap<String,MutableList<Atom>>()
		for (atom in atoms) {
			dupesByName.getOrPut(atom.name) { ArrayList() }.add(atom)
		}

		// find the actual duplicates and flag them for user review
		dupesByName.values
			.filter { it.size > 1 }
			.forEach {

				// make a friendly description for the location
				val resId = res?.id
				val chainId = chain?.id
				val location = if (resId != null && chainId != null) {
					"$chainId$resId"
				} else {
					null
				}

				dupeAtoms.add(DupeAtoms(view, mol, res, it, location))
			}
	}
}
