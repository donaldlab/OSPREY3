package edu.duke.cs.osprey.molscope.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.*
import cuchaz.kludge.vulkan.Extent2D
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.molscope.render.Camera
import edu.duke.cs.osprey.molscope.render.HoverEffects
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import org.joml.Vector2fc
import kotlin.math.abs
import kotlin.math.atan2


class CameraTool : SlideFeature {

	override val id = FeatureId("camera")

	private val winState = WindowState()
	private var hoverEffects = null as HoverEffects.Writer?

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Camera")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		// handle mouse interactions, even if the window isn't open
		if (slidewin.mouseLeftClick) {
			click(slidewin)
		}
		if (slidewin.mouseLeftDrag) {
			drag(slide, slidewin)
		}
		if (slidewin.mouseWheelDelta != 0f) {
			wheel(slidewin)
		}

		winState.render(
			onOpen = {
				// add the hover effect
				hoverEffects = slidewin.hoverEffects.writer().apply {
					effect = hoverEffect
				}
			},
			whenOpen = {

				// draw the window
				window("Camera##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					// show camera properties
					var changed = false
					changed = sliderFloat("Magnification", Ref.of(slidewin.camera::magnification), 1f, 200f, "%.1fx", power=4f) || changed
					changed = sliderFloat("Closeness", Ref.of(slidewin.camera::closeness), 0.001f, 1f, "%.3f") || changed
					changed = sliderFloat("View Distance", Ref.of(slidewin.camera::viewDistance), 1f, 400f, "%.1f", power=4f) || changed
					if (changed) {
						slidewin.camera.changed()
					}
				}
			},
			onClose = {

				// cleanup effects
				hoverEffects?.close()
				hoverEffects = null
			}
		)
	}

	override fun contextMenu(contextMenu: ContextMenu, slide: Slide.Locked, slidewin: SlideCommands, target: ViewIndexed) {

		if (!winState.isOpen) {
			return
		}

		// get the atom, if any
		val view = target.view as? MoleculeRenderView ?: return
		val mol = view.currentMol
		val atom = target.target as? Atom ?: return

		contextMenu.add {

			// show details about the atom
			text("Atom: ${atom.name} ${atom.element}")
			indent(10f)

			// show the position
			text("pos: ")
			sameLine()
			val coordsText = Commands.TextBuffer.of("%.3f,%.3f,%.3f".format(atom.pos.x, atom.pos.y, atom.pos.z()))
			inputText("", coordsText, IntFlags.of(Commands.InputTextFlags.ReadOnly))

			// show the molecule
			text("mol: ${mol.name}")

			// show the residue, if any
			if (mol is Polymer) {
				mol.chains
					.mapNotNull { chain ->
						chain.residues.find { atom in it.atoms }
					}
					.firstOrNull()
					?.let { res ->
						text("res: ${res.id} ${res.type}")
					}
			}

			unindent(10f)

			// show a button to center the camera on the atom
			if (button("Center")) {
				closeCurrentPopup()
				slidewin.camera.apply {
					lookAt(atom.pos.toFloat(), slide.views)
					changed()
				}
			}
		}
	}

	private var cameraRotator: Camera.Rotator? = null
	private var dragStartAngle = 0f
	private var dragMode: DragMode = DragMode.RotateXY

	private enum class DragMode {
		RotateXY,
		RotateZ
	}

	private fun getDragAngle(mouseOffset: Vector2fc, extent: Extent2D): Float {
		return atan2(
			// if the camera z points away, and y points up, then x points left,
			// so flip the screen x coords where x points right
			mouseOffset.y - extent.height.toFloat()/2,
			extent.width.toFloat()/2 - mouseOffset.x
		)
	}

	private fun click(slidewin: SlideCommands) {

		// reset the camera rotator
		if (cameraRotator?.camera != slidewin.camera) {
			cameraRotator = slidewin.camera.Rotator()
		}
		cameraRotator?.capture()

		// get the normalized click dist from center
		val w = slidewin.extent.width.toFloat()
		val h = slidewin.extent.height.toFloat()
		val dx = abs(slidewin.mouseOffset.x*2f/w - 1f)
		val dy = abs(slidewin.mouseOffset.y*2f/h - 1f)

		// pick the drag mode based on the click pos
		// if we're near the center, rotate about xy
		// otherwise, rotate about z
		val cutoff = 0.8
		dragMode = if (dx < cutoff && dy < cutoff) {
			DragMode.RotateXY
		} else {
			dragStartAngle = getDragAngle(slidewin.mouseOffset, slidewin.extent)
			DragMode.RotateZ
		}
	}

	private fun drag(slide: Slide.Locked, slidewin: SlideCommands) {

		// apply the drag rotations
		cameraRotator?.apply {
			q.identity()
			when (dragMode) {
				DragMode.RotateXY -> {
					// if the camera z points away, and y points up, then x points left,
					// so flip the screen x coords where x points right
					q.rotateAxis(-slidewin.mouseLeftDragDelta.x/100f, up)
					q.rotateAxis(slidewin.mouseLeftDragDelta.y/100f, side)
				}
				DragMode.RotateZ -> {
					q.rotateAxis(getDragAngle(slidewin.mouseOffset, slidewin.extent) - dragStartAngle, look)
				}
			}
			update()

			slidewin.camera.changed()
		}
	}

	private fun wheel(slidewin: SlideCommands) {

		// adjust the magnification
		slidewin.camera.magnification *= 1f + slidewin.mouseWheelDelta/10f

		slidewin.camera.changed()
	}
}

private val hoverEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
	200u, 200u, 200u
)
