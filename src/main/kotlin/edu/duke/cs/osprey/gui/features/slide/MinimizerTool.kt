package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.tools.toFloat
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.forcefield.amber.*
import org.joml.Vector3d
import java.util.*
import java.util.concurrent.atomic.AtomicBoolean


class MinimizerTool : SlideFeature {

	override val id = FeatureId("edit.minimize")

	private val winState = WindowState()

	private class MolInfo(val view: MoleculeRenderView) {
		val mol = view.molStack.originalMol
		// keep the heavy atoms from wandering, but let the hydrogens adjust freely
		val restrainedAtoms = mol.atoms.filter { it.element != Element.Hydrogen }
		val minInfo = MinimizerInfo(mol, restrainedAtoms)
		val pSelected = Ref.of(true)
	}

	private val molInfos = IdentityHashMap<Molecule,MolInfo>()
	private val pNumSteps = Ref.of(100)
	private var job: Job? = null

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Minimize")) {
			winState.isOpen = true
		}
	}

	private fun Slide.Locked.molViews() = views.mapNotNull { it as? MoleculeRenderView }

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val molViews = slide.molViews()

		winState.render(
			onOpen = {
				// clear any old mol infos
				molInfos.clear()
				for (view in molViews) {
					molInfos[view.molStack.originalMol] = MolInfo(view)
				}
			},
			whenOpen = {

				// draw the window
				window("Minimize##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					for ((index, info) in molInfos.values.withIndex()) {
						withId(index) {

							val mol = info.mol

							// show the molecule type
							checkbox("$mol", info.pSelected)

							// show a context menu to center the camera on the molecule
							popupContextItem("centerCamera") {

								if (button("Center Camera")) {

									// center on the molecule centroid
									val center = Vector3d().apply {
										mol.atoms.forEach { add(it.pos) }
										div(mol.atoms.size.toDouble())
									}
									slidewin.camera.lookAt(center.toFloat(), slide.views)
									slidewin.camera.changed()

									closeCurrentPopup()
								}
							}

							// add buttons to set coords
							val unminimizedCoords = info.minInfo.unminimizedCoords
							val minimizedCoords = info.minInfo.minimizedCoords
							text("Set Coords:")
							sameLine()
							disabledIf(unminimizedCoords == null) {
								if (button("Unminimized") && unminimizedCoords != null) {
									info.minInfo.setCoords(unminimizedCoords)
									info.view.moleculeChanged()
								}
							}
							sameLine()
							disabledIf(minimizedCoords == null) {
								if (button("Minimized") && minimizedCoords != null) {
									info.minInfo.setCoords(minimizedCoords)
									info.view.moleculeChanged()
								}
							}

							// let the entries breathe a little
							spacing()
							spacing()
							separator()
							spacing()
							spacing()
						}
					}

					val job = job
					if (job == null) {

						sliderInt("Num Steps", pNumSteps, 1, 1000)

						if (button("Minimize Selected Molecules")) {
							this@MinimizerTool.job = Job(
								molInfos.values.filter { it.pSelected.value },
								pNumSteps.value
							)
						}

					} else {

						text("Minimizing...")
						// TODO: show minimization progress
						// TODO: cancel button?

						if (job.isFinished.get()) {

							// report any errors
							job.throwable?.let { t ->
								slidewin.showExceptions { throw t }
							}

							// set the minmized coords now
							for (info in job.infos) {
								info.minInfo.minimizedCoords?.let { info.minInfo.setCoords(it) }
								info.view.moleculeChanged()
							}

							// cleanup the finished job
							this@MinimizerTool.job = null
						}
					}
				}
			},
			onClose = {
				// cleanup our mol infos
				molInfos.clear()
			}
		)
	}

	private class Job(
		val infos: List<MolInfo>,
		val numSteps: Int
	) {

		val isFinished = AtomicBoolean(false)
		var throwable: Throwable? = null

		val thread =
			Thread {

				try {
					infos
						.map { it.minInfo }
						.minimizeBlocking(numSteps)
				} catch (t: Throwable) {
					throwable = t
				}

				// signal the thread is done
				isFinished.set(true)
			}
			.apply {
				name = "Minimizer"
			}
			.start()
	}
}
