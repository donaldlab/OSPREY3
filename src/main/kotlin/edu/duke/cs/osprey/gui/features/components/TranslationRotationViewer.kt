package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.ByteFlags
import cuchaz.kludge.tools.Ref
import cuchaz.kludge.tools.toRadians
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.MoleculeMaps
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.motions.TranslationRotation
import edu.duke.cs.osprey.gui.tools.nextFloatIn
import kotlin.random.Random


class TranslationRotationViewer private constructor(
	val transRot: TranslationRotation,
	val maxTranslationDist: Float,
	val maxRotationDegrees: Float,
	val view: MoleculeRenderView,
	val mapsToCopy: MoleculeMaps
) : MotionViewer {

	companion object {

		fun make(desc: TranslationRotation.MolDescription, molCopy: Molecule, molMaps: MoleculeMaps, view: MoleculeRenderView): TranslationRotationViewer {
			return TranslationRotationViewer(
				desc.copyTo(molCopy, molMaps.atoms).make(),
				desc.maxTranslationDist.toFloat(),
				desc.maxRotationDegrees.toFloat(),
				view,
				molMaps
			)
		}
	}

	private val pPsi = Ref.of(0.0f)
	private val pTheta = Ref.of(0.0f)
	private val pPhi = Ref.of(0.0f)
	private val px = Ref.of(0.0f)
	private val py = Ref.of(0.0f)
	private val pz = Ref.of(0.0f)

	private val rmax = maxRotationDegrees
	private val rmin = -rmax
	private val tmax = maxTranslationDist
	private val tmin = -tmax

	override val label = "Translation and Rotation of ${transRot.mol}"

	private var renderEffects = view.renderEffects.writer().apply {

		// show all the translatable atoms
		for (atom in transRot.mol.atoms) {
			this[atom] = selectedEffect
		}
	}

	private fun updateMol(view: MoleculeRenderView) {
		transRot.set(
			pPsi.value.toDouble().toRadians(),
			pTheta.value.toDouble().toRadians(),
			pPhi.value.toDouble().toRadians(),
			px.value.toDouble(),
			py.value.toDouble(),
			pz.value.toDouble()
		)
		view.moleculeChanged()
	}

	override fun gui(imgui: Commands) = imgui.run {

		text("Tait-Bryan Rotation:")
		if (sliderFloat("Psi (X)", pPsi, rmin, rmax, "%.3f")) {
			updateMol(view)
		}
		if (sliderFloat("Theta (Y)", pTheta, rmin, rmax, "%.3f")) {
			updateMol(view)
		}
		if (sliderFloat("Phi (Z)", pPhi, rmin, rmax, "%.3f")) {
			updateMol(view)
		}

		text("Cartesian Translation:")
		if (sliderFloat("X", px, tmin, tmax, "%.3f")) {
			updateMol(view)
		}
		if (sliderFloat("Y", py, tmin, tmax, "%.3f")) {
			updateMol(view)
		}
		if (sliderFloat("Z", pz, tmin, tmax, "%.3f")) {
			updateMol(view)
		}
	}

	override fun jiggle(rand: Random) {

		pPsi.value = rand.nextFloatIn(rmin, rmax)
		pTheta.value = rand.nextFloatIn(rmin, rmax)
		pPhi.value = rand.nextFloatIn(rmin, rmax)

		px.value = rand.nextFloatIn(tmin, tmax)
		py.value = rand.nextFloatIn(tmin, tmax)
		pz.value = rand.nextFloatIn(tmin, tmax)

		updateMol(view)
	}

	override fun reset() {

		// remove any render effects
		renderEffects.close()
	}

	override fun mapAtomToOriginalMol(atom: Atom): Atom = mapsToCopy.atoms.getAOrThrow(atom)
}

private val selectedEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Outset),
	0u, 255u, 0u
)
