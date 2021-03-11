package edu.duke.cs.osprey.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.*
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.*
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Element
import edu.duke.cs.osprey.molscope.molecule.Polymer
import edu.duke.cs.osprey.molscope.render.HoverEffects
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.dof.DihedralRotation
import edu.duke.cs.osprey.tools.Protractor
import edu.duke.cs.osprey.gui.forcefield.amber.*
import kotlinx.coroutines.runBlocking
import org.joml.Vector3d
import java.util.*
import kotlin.math.cos
import kotlin.math.sin


class ProtonationEditor : SlideFeature {

	override val id = FeatureId("edit.protonation")

	private val winState = WindowState()
	private val clickTracker = ClickTracker()
	private var selection: Selection? = null
	private var hoverEffects = null as HoverEffects.Writer?
	private val renderEffects = IdentityHashMap<MoleculeRenderView,MoleculeRenderEffects.Writer>()

	private fun Slide.Locked.molViews() = views.filterIsInstance<MoleculeRenderView>()

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Protonation")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		val views = slide.molViews()

		winState.render(
			onOpen = {

				// add the hover effect
				hoverEffects = slidewin.hoverEffects.writer().apply {
					effect = hoverEffect
				}

				// init the render effects
				for (view in views) {
					renderEffects[view] = view.renderEffects.writer()
				}
			},
			whenOpen = {

				// did we click anything?
				if (clickTracker.clicked(slidewin)) {

					// clear any previous selection
					renderEffects.values.forEach { it.clear() }
					selection = null

					// select the heavy atom from the click, if any
					slidewin.mouseTarget?.let { target ->
						(target.view as? MoleculeRenderView)?.let { view ->
							(target.target as? Atom)?.let { atom ->
								if (atom.element != Element.Hydrogen) {
									selection = Selection(view, atom)
								}
							}
						}
					}
				}

				// draw the window
				window("Protonation Editor##${slide.name}", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					text("Tools:")
					indent(10f)

					if (button("Clear all Hydrogens")) {
						clearAll(views)
					}

					// TODO: Reduce

					if (button("Add Hydrogens automatically")) {
						slidewin.showExceptions {
							autoProtonate(views)
						}
					}
					sameLine()
					infoTip("""
						|This tool adds Hydrogen atoms to heavy atoms based on inferred
						|forcefield atom types and bond types. Atom and bond type inference
						|on unprotonated molecules is very error-prone, so feel free to use
						|the fine-grained editing tools to add any missing Hydrogen atoms,
						|removea any extraneous Hydrogen atoms, or change hybridizations and geometry.
					""".trimMargin())

					unindent(10f)

					val selection = selection

					// show the selected atom
					text("Selected:")
					child("selected", 300f, 30f, true) {
						if (selection != null) {
							text(selection.atom.name)
						} else {
							text("(Click a heavy atom to show protonation options.)")
						}
					}

					// show the nearby atoms
					if (selection != null) {

						// show a list box to pick the protonation state
						text("Protonation states:")

						val numItems = selection.protonations.size + 1
						listBox("", numItems) {

							// always add an option for no hydrogens
							if (selectable("0 H", selection.current == null)) {
								selection.set(null)
							}

							for (protonation in selection.protonations) {
								if (selectable("${protonation.numH} H, ${protonation.hybridization}", selection.current == protonation)) {
									slidewin.showExceptions {
										selection.set(protonation)
									}
								}
							}
						}

						// if hydrogens are rotatable, show a slider to pick the dihedral angle
						if (selection.rotator?.rotationH != null) {
							val rotator = selection.rotator

							spacing()

							text("Rotation")
							if (sliderFloat("Dihedral angle", rotator.pDihedral, -180f, 180f, "%.1f")) {
								rotator.set()
							}
							sameLine()
							infoTip("Ctrl-click to type angle exactly")

							if (rotator.rotationHeavies.size > 1) {
								if (listBoxHeader("Anchor Heavy Atom", rotator.rotationHeavies.size)) {
									for (atom in rotator.rotationHeavies) {
										if (selectable(atom.name, rotator.rotationHeavy == atom)) {
											rotator.rotationHeavy = atom
										}
									}
									listBoxFooter()
								}
							}
						}
					}
				}
			},
			onClose = {

				// remove the hover effect
				hoverEffects?.close()
				hoverEffects = null

				// clear any leftover selections when the window closes
				renderEffects.values.forEach { it.close() }
				renderEffects.clear()
				selection = null
			}
		)
	}

	private fun clearAll(views: List<MoleculeRenderView>) {
		for (view in views) {
			val mol = view.molStack.originalMol
			mol.deprotonate()
			view.moleculeChanged()
		}
	}

	private fun autoProtonate(views: List<MoleculeRenderView>) {
		for (view in views) {
			val mol = view.molStack.originalMol
			val protonation = runBlocking { mol.inferProtonation() }
			for ((heavyAtom, hAtom) in protonation) {
				mol.apply {
					atoms.add(hAtom)
					bonds.add(heavyAtom, hAtom)
					if (this is Polymer) {
						findResidue(heavyAtom)?.atoms?.add(hAtom)
					}
				}
			}
			view.moleculeChanged()
		}
	}

	private inner class Selection(val view: MoleculeRenderView, val atom: Atom) {

		val mol = view.molStack.originalMol

		// get the protonations
		val protonations = mol.protonations(atom)

		var current: Protonation? = run {
			val bondedAtoms = mol.bonds.bondedAtoms(atom)
			val numHeavy = bondedAtoms.count { it.element != Element.Hydrogen }
			val numH = bondedAtoms.count { it.element == Element.Hydrogen }
			protonations
				.filter { it.numHeavy == numHeavy && it.numH == numH }
				// TODO: try to detect the hybridization?
				//  (so we don't have to return null when we get multiple answers)
				.takeIf { it.size == 1 }
				?.first()
		}

		inner class Rotator(val bondedHeavy: Atom, val rotationHeavies: List<Atom>) {

			var rotationHeavy: Atom = rotationHeavies.first()
				set(value) {
					field = value
					update()
					updateSelectionEffects()
				}

			var rotationH: Atom? = pickHydrogen()

			// pick the hydrogen atom to define the dihedral angle
			private fun pickHydrogen() =
				mol.bonds.bondedAtoms(atom)
					.filter { it.element == Element.Hydrogen }
					.minBy { it.toInt() }

			val pDihedral = Ref.of(measureDihedral())

			private fun measureDihedral(): Float {
				val rotationH = rotationH ?: return 0f
				return measureDihedral(
					rotationHeavy,
					bondedHeavy,
					atom,
					rotationH
				).toFloat()
			}

			fun Atom.toInt() = name.filter { it.isDigit() }.toIntOrNull() ?: 1

			fun Vector3d.toArray() = doubleArrayOf(x, y, z)
			fun Vector3d.fromArray(array: DoubleArray) = set(array[0], array[1], array[2])

			fun measureDihedral(a: Atom, b: Atom, c: Atom, d: Atom) =
				Protractor.measureDihedral(arrayOf(
					a.pos.toArray(),
					b.pos.toArray(),
					c.pos.toArray(),
					d.pos.toArray()
				))

			fun update() {

				// update the rotation hydrogen
				rotationH = pickHydrogen()

				// measure the current dihedral angle
				pDihedral.value = measureDihedral()
			}

			fun set() {

				// compute the target dihedral
				val (targetSin, targetCos) = pDihedral.value.toDouble().toRadians()
					.let { sin(it) to cos(it) }

				// measure the current dihedral
				val (currentSin, currentCos) = measureDihedral().toDouble().toRadians()
					.let { sin(it) to cos(it) }

				// calc the dihedral rotation as a rigid body transformation relative to the current pose
				val dsin = targetSin*currentCos - targetCos*currentSin
				val dcos = targetCos*currentCos + targetSin*currentSin
				val rotation = DihedralRotation(
					bondedHeavy.pos.toArray(),
					atom.pos.toArray(),
					dsin, dcos
				)

				// rotate all the hydrogens
				mol.bonds.bondedAtoms(atom)
					.filter { it.element == Element.Hydrogen }
					.forEach { h ->
						val coords = h.pos.toArray()
						rotation.transform(coords, 0)
						h.pos.fromArray(coords)
					}

				view.moleculeChanged()
			}
		}

		// if the atom has only one bonded heavy atom,
		// and that heavy atom has at least one other bonded heavy atom,
		// then the hydrogens are rotatable
		val rotator: Rotator? = run {

			val bondedHeavies = mol.bonds.bondedAtoms(atom)
				.filter { it.element != Element.Hydrogen }
			if (bondedHeavies.size == 1) {
				val bondedHeavy = bondedHeavies.first()

				val rotationHeavies = mol.bonds.bondedAtoms(bondedHeavy)
					.filter { it.element != Element.Hydrogen && it != atom }
				if (rotationHeavies.isNotEmpty()) {
					return@run Rotator(bondedHeavy, rotationHeavies)
				}
			}

			return@run null
		}

		init {
			updateSelectionEffects()
		}

		fun set(protonation: Protonation?) {

			// update the selection
			if (current == protonation) {
				return
			}
			current = protonation

			// update the molecule
			if (protonation != null) {
				runBlocking { mol.protonate(atom, protonation) }
			} else {
				mol.deprotonate(atom)
			}
			view.moleculeChanged()

			rotator?.update()

			updateSelectionEffects()
		}

		private fun updateSelectionEffects() {

			val effects = renderEffects[view] ?: return

			effects.clear()

			// add the selection effect
			effects[atom] = selectedEffect

			// highlight the dihedral atoms if needed
			rotator?.rotationH?.let { rotationH ->
				effects[listOf(
					rotator.rotationHeavy,
					rotator.bondedHeavy,
					atom,
					rotationH
				)] = rotationEffect
			}
		}
	}
}


private val hoverEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
	200u, 200u, 200u
)
private val selectedEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
	255u, 255u, 255u
)
private val rotationEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Inset),
	100u, 200u, 100u
)
