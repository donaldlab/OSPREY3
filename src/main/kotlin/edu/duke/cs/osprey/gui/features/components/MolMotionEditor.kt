package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.ClickTracker
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.gui.infoTip
import edu.duke.cs.osprey.molscope.gui.enabledIf
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.view.MoleculeRenderStack
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.features.slide.ConformationEditor
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.motions.MolMotion
import edu.duke.cs.osprey.gui.motions.TranslationRotation


class MolMotionEditor(
	val molInfo: ConformationEditor.MolInfo,
	var desc: MolMotion.Description?,
	val onClose: () -> Unit
) {

	/**
	 * Since we're editing a thing that might not exist yet (ie desc might be null),
	 * we can't use that thing to establish identity for the GUI window.
	 * So just pick an unique number to use as the identity instead.
	 * We probably won't need to edit more than ~2 billion motions in single session, right?
	 */
	companion object {
		private var nextId = 0
		fun makeId() = nextId++
	}
	private val id = makeId()

	private val winState = WindowState()
		.apply { pOpen.value = true }

	private var info: MotionInfo? = null
	private var viewer: MotionViewer? = null
	private var stackedMol: MoleculeRenderStack.StackedMol? = null


	private fun Slide.Locked.findView() =
		views
			.filterIsInstance<MoleculeRenderView>()
			.find { it.molStack.originalMol === molInfo.mol }
			?: throw Error("can't edit molecule motion, molecule has no render view")

	fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		winState.render(
			onOpen = {

				val view = slide.findView()

				// init the motion info
				when (desc) {

					// make info for an existing motion
					is DihedralAngle.MolDescription -> resetInfo(view, MotionInfo.DihedralAngleInfo(this@MolMotionEditor))
					is TranslationRotation.MolDescription -> resetInfo(view, MotionInfo.TranslationRotationInfo(this@MolMotionEditor))

					// by default, make a new translation/rotation motion and its info
					null -> resetInfo(view, MotionInfo.TranslationRotationInfo(this@MolMotionEditor)).initDefault(view)
				}
			},
			whenOpen = {

				val view = slide.findView()
				val desc = desc

				// draw the window
				window("Edit Molecule Motion##$id", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					// show options to pick the motion type
					text("Motion type:")
					indent(20f)

					enabledIf(molInfo.mol.atoms.size >= 4) {
						if (radioButton("Dihedral Angle", desc is DihedralAngle.MolDescription)) {
							resetInfo(view, MotionInfo.DihedralAngleInfo(this@MolMotionEditor)).initDefault(view)
						}
					}

					if (radioButton("Translation & Rotation", desc is TranslationRotation.MolDescription)) {
						resetInfo(view, MotionInfo.TranslationRotationInfo(this@MolMotionEditor)).initDefault(view)
					}

					unindent(20f)

					spacing()
					spacing()
					spacing()

					// show the editing gui for the motion
					text("Settings:")
					indent(20f)
					info?.gui(imgui, slidewin, view)
					unindent(20f)

					// show the motion viewer, if possible
					val viewer = viewer
					if (viewer != null) {

						spacing()
						spacing()
						spacing()

						text("Try it out:")
						indent(20f)
						viewer.gui(imgui)
						unindent(20f)
					}
				}
			},
			onClose = {

				// cleanup
				viewer?.reset()
				viewer = null
				info?.close()
				info = null
				stackedMol?.pop()
				stackedMol = null

				onClose()
			}
		)
	}

	private fun resetInfo(view: MoleculeRenderView, info: MotionInfo): MotionInfo {

		// cleanup the old info, if any
		this.info?.close()

		// set the new info
		this.info = info

		// initialize the viewer
		resetViewer(view)

		return info
	}

	fun resetMotion(view: MoleculeRenderView, desc: MolMotion.Description) {

		val index = molInfo.motions.indexOfFirst { it === this.desc }
		if (index >= 0) {
			// if the old motion exists, replace it
			molInfo.motions[index] = desc
		} else {
			// otherwise, append to the list
			molInfo.motions.add(desc)
		}

		this.desc = desc

		resetViewer(view)
	}

	private fun resetViewer(view: MoleculeRenderView) {

		// cleanup the old viewer, if needed
		viewer?.reset()
		stackedMol?.pop()

		// make a copy of the molecule so the viewer can modify it
		val (molCopy, molMaps) = molInfo.mol.copyWithMaps()
		stackedMol = view.molStack.push(molCopy)

		// start the new viewer
		val desc = desc
		viewer = when(desc) {
			is DihedralAngle.MolDescription -> DihedralAngleViewer.make(desc, molCopy, molMaps, view)
			is TranslationRotation.MolDescription -> TranslationRotationViewer.make(desc, molCopy, molMaps, view)
			else -> null
		}
	}

	private fun mapAtomToOriginalMol(atom: Atom) =
		viewer?.mapAtomToOriginalMol(atom) ?: atom

	private sealed class MotionInfo : AutoCloseable {

		abstract fun initDefault(view: MoleculeRenderView)
		abstract fun gui(imgui: Commands, slidewin: SlideCommands, view: MoleculeRenderView)

		override fun close() {
			// nothing to cleanup by default
		}

		class DihedralAngleInfo(val editor: MolMotionEditor) : MotionInfo() {

			companion object {
				const val defaultRadiusDegrees = 9.0
			}

			private val desc get() = editor.desc as? DihedralAngle.MolDescription

			override fun initDefault(view: MoleculeRenderView) {

				// start with arbitrary atoms, we'll change them later
				val mol = editor.molInfo.mol
				val a = mol.atoms[0]
				val b = mol.atoms[1]
				val c = mol.atoms[2]
				val d = mol.atoms[3]
				val initialDegrees = DihedralAngle.measureDegrees(a.pos, b.pos, c.pos, d.pos)

				editor.resetMotion(view, DihedralAngle.MolDescription(
					mol,
					a, b, c, d,
					initialDegrees - defaultRadiusDegrees,
					initialDegrees + defaultRadiusDegrees
				))
			}

			private enum class AtomSel {

				A, B, C, D;

				fun next() =
					values()[(ordinal + 1) % 4]
			}

			private val pRadiusDegrees = Ref.of((desc?.radiusDegrees ?: defaultRadiusDegrees).toFloat())
			private var currentAtom = AtomSel.A
			private val clickTracker = ClickTracker()

			override fun gui(imgui: Commands, slidewin: SlideCommands, view: MoleculeRenderView) = imgui.run {

				val desc = desc ?: return

				text("Angle radius (in degrees)")
				sameLine()
				infoTip("""
					|This value specifies the half-width (ie radius) of the interval of allowed angles.
				""".trimMargin())
				if (sliderFloat("##dihedralRadius", pRadiusDegrees, 0f, 180f,"%.1f")) {
					editor.resetMotion(view, DihedralAngle.MolDescription(
						desc.mol, desc.a, desc.b, desc.c, desc.d, radiusDegrees = pRadiusDegrees.value.toDouble()
					))
				}
				sameLine()
				infoTip("Ctrl-click to type a precise value")

				// atom pickers
				if (radioButton("Atom A: ${desc.a.name}###atoma", currentAtom == AtomSel.A)) {
					currentAtom = AtomSel.A
				}
				if (radioButton("Atom B: ${desc.b.name}###atomb", currentAtom == AtomSel.B)) {
					currentAtom = AtomSel.B
				}
				if (radioButton("Atom C: ${desc.c.name}###atomc", currentAtom == AtomSel.C)) {
					currentAtom = AtomSel.C
				}
				if (radioButton("Atom D: ${desc.d.name}###atomd", currentAtom == AtomSel.D)) {
					currentAtom = AtomSel.D
				}

				// are we hovering over an atom?
				val hoverAtom: Atom? =
					(slidewin.mouseTarget
						?.takeIf { it.view === view }
						?.target as? Atom)
						// this atom can be from the viewer's copy of the molecule
						// so translate the atom to the original molecule
						?.let { editor.mapAtomToOriginalMol(it) }

				// if we clicked an atom, update the motion
				if (hoverAtom != null && clickTracker.clicked(slidewin)) {
					val newdesc = when (currentAtom) {
						AtomSel.A -> DihedralAngle.MolDescription(
							desc.mol, hoverAtom, desc.b, desc.c, desc.d, radiusDegrees = pRadiusDegrees.value.toDouble()
						)
						AtomSel.B -> DihedralAngle.MolDescription(
							desc.mol, desc.a, hoverAtom, desc.c, desc.d, radiusDegrees = pRadiusDegrees.value.toDouble()
						)
						AtomSel.C -> DihedralAngle.MolDescription(
							desc.mol, desc.a, desc.b, hoverAtom, desc.d, radiusDegrees = pRadiusDegrees.value.toDouble()
						)
						AtomSel.D -> DihedralAngle.MolDescription(
							desc.mol, desc.a, desc.b, desc.c, hoverAtom, radiusDegrees = pRadiusDegrees.value.toDouble()
						)
					}
					editor.resetMotion(view, newdesc)
					currentAtom = currentAtom.next()
				}
			}
		}

		class TranslationRotationInfo(val editor: MolMotionEditor) : MotionInfo() {

			companion object {
				const val defaultDist = 1.2 // angstroms
				const val defaultDegrees = 5.0
			}

			private val desc get() = editor.desc as? TranslationRotation.MolDescription

			override fun initDefault(view: MoleculeRenderView) {
				editor.resetMotion(view, TranslationRotation.MolDescription(
					editor.molInfo.mol,
					defaultDist,
					defaultDegrees
				))
			}

			private val pDist = Ref.of(desc?.maxTranslationDist ?: defaultDist)
			private val pDegrees = Ref.of(desc?.maxRotationDegrees ?: defaultDegrees)

			override fun gui(imgui: Commands, slidewin: SlideCommands, view: MoleculeRenderView) = imgui.run {

				val desc = editor.desc as? TranslationRotation.MolDescription ?: return

				if (inputDouble("Max Distance (Angstroms)", pDist)) {
					editor.resetMotion(view, desc.copy(maxTranslationDist = pDist.value))
				}

				if (inputDouble("Max Rotation (Degrees)", pDegrees)) {
					editor.resetMotion(view, desc.copy(maxRotationDegrees = pDegrees.value))
				}
			}
		}
	}
}
