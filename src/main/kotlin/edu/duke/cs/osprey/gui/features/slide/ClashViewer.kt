package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.filter
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.combineForPDB
import edu.duke.cs.osprey.gui.io.toPDB
import edu.duke.cs.osprey.gui.io.toVector3d
import edu.duke.cs.osprey.gui.view.ProbeView
import edu.duke.cs.osprey.service.services.ClashesRequest
import kotlinx.coroutines.runBlocking


class ClashViewer : SlideFeature {

	override val id = FeatureId("view.clashes")

	private val winState = WindowState()

	private val view = ProbeView()

	private class Counts(val pVisible: Ref<Boolean>) {
		var dots = 0
		var vectors = 0
	}
	private val countsByType = HashMap<String,Counts>()

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Clashes")) {
			winState.isOpen = true
		}
	}

	private fun Slide.Locked.molViews() = views.mapNotNull { it as? MoleculeRenderView }

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val molViews = slide.molViews()

		winState.render(
			onOpen = {
				slide.views.add(view)
				slidewin.showExceptions {
					loadClashes(molViews)
				}
			},
			whenOpen = {

				// draw the window
				setNextWindowSize(300f, 0f)
				window("Clashes##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					if (button("Refresh")) {
						slidewin.showExceptions {
							unloadClashes()
							loadClashes(molViews)
						}
					}

					// show the counts, with toggles to show/hide
					if (countsByType.isNotEmpty()) {

						columns(2)
						for ((type, counts) in countsByType) {

							checkbox(type, counts.pVisible)
							nextColumn()

							if (counts.dots > 0) {
								text("${counts.dots} dots")
							} else if (counts.vectors > 0) {
								text("${counts.vectors} vectors")
							}
							nextColumn()
						}
						columns(1)

					} else {
						text("(no clash information)")
					}
				}

			},
			onClose = {
				unloadClashes()
				slide.views.remove(view)
			}
		)
	}

	private fun loadClashes(views: List<MoleculeRenderView>) {

		// combine all the molecules into one PDB
		val pdb = views.map { view ->
				// grab just the selected part of the molecule
				view.selector.filter(view.currentMol)
			}
			.combineForPDB("combined").first
			.toPDB()

		// run probe
		val results = runBlocking { OspreyService.clashes(ClashesRequest(pdb)) }
		view.groups = results.groups
			.mapValues { (id, group) ->
				ProbeView.Group(
					id,
					group.dots.mapValues { (_, dots) -> dots.map { it.toVector3d() } },
					group.vectors.mapValues { (_, vectors) ->
						vectors.map { (a, b) -> a.toVector3d() to b.toVector3d() }
					}
				)
			}

		// update the counts
		countsByType.clear()
		for ((type, group) in results.groups) {
			countsByType
				.getOrPut(type) {
					Counts(Ref.of(
						getter = { view.visibility[type] },
						setter = { value -> view.visibility[type] = value }
					))
				}
				.apply {
					dots += group.dots.values.sumBy { it.size }
					vectors += group.vectors.values.sumBy { it.size }
				}
		}
	}

	private fun unloadClashes() {

		// clear probe results
		view.groups = null
	}
}