package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.tools.toFloat
import cuchaz.kludge.tools.toString
import edu.duke.cs.osprey.gui.prep.DuplicateAtoms
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.gui.infoTip
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import java.util.*


class DuplicateAtomsEditor : SlideFeature {

	override val id = FeatureId("edit.duplicateAtoms")

	private val winState = WindowState()

	private data class MolInfo(val view: MoleculeRenderView) {

		val mol = view.currentMol
		val dupes = DuplicateAtoms(mol)
		val groupInfos = dupes.map { GroupInfo(it) }.toMutableList()

		inner class GroupInfo(val group: DuplicateAtoms.AtomGroup) {

			val pIncluded = group.atoms.map { Ref.of(true) }

			fun update(atomi: Int) {

				// add or remove the atom
				val pIncluded = pIncluded[atomi]
				if (pIncluded.value) {
					group.remove(atomi)
				} else {
					group.add(atomi)
				}

				view.moleculeChanged()
			}

			fun resolve() {
				if (pIncluded.count { it.value } == 1) {
					groupInfos.remove(this)
				}
			}
		}

	}

	private val molInfos = IdentityHashMap<MoleculeRenderView,MolInfo>()
	private val toResolve = ArrayList<MolInfo.GroupInfo>()

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

					var groupi = 0

					for (molInfo in molInfos.values) {

						spacing()
						spacing()
						spacing()
						text(molInfo.mol.toString())

						indent(20f)

						if (molInfo.groupInfos.isNotEmpty()) {
							for (groupInfo in molInfo.groupInfos) {

								spacing()
								spacing()
								spacing()

								text("Duplicated: ${groupInfo.group}:")
								indent(20f)

								// show checkboxes to toggle the atoms on/off
								for ((atomi, pair) in groupInfo.group.atoms.zip(groupInfo.pIncluded).withIndex()) {
									val (atom, pIncluded) = pair

									val label = "${atom.element.symbol} @ ${atom.pos.toString(6)}##$atomi"
									if (checkbox(label, pIncluded)) {
										groupInfo.update(atomi)
									}

									// show a context menu for the checkbox
									popupContextItem("group:$groupi") {

										if (button("Center Camera")) {
											slidewin.camera.lookAt(atom.pos.toFloat(), slide.views)
											slidewin.camera.changed()
											closeCurrentPopup()
										}
									}
								}

								if (button("Resolve##$groupi")) {
									toResolve.add(groupInfo)
								}
								sameLine()
								infoTip("""
									|If the atom is still duplicated, this button does nothing.
									|Otherwise, this button will remove this entry from the duplicated atoms list.
								""".trimMargin())

								unindent(20f)

								groupi += 1
							}

						} else {
							text("(No duplicate atoms detected)")
						}

						unindent(20f)
					}
				}

				// handle the resolutions outside of the loop to avoid concurrent modifications
				for (groupInfo in toResolve) {
					groupInfo.resolve()
				}
				toResolve.clear()
			},
			onClose = {
				// clear any previous state
				molInfos.clear()
			}
		)
	}

	private fun detectDupeAtoms(views: List<MoleculeRenderView>) {

		molInfos.clear()

		// look for duplicated atoms
		for (view in views) {
			molInfos[view] = MolInfo(view)
		}
	}
}
