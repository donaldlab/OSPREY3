package edu.duke.cs.osprey.gui.features.components

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.ByteFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.MoleculeMaps
import edu.duke.cs.osprey.molscope.render.RenderEffect
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.gui.motions.DihedralAngle
import edu.duke.cs.osprey.gui.prep.Assignments
import edu.duke.cs.osprey.gui.tools.nextFloatIn
import kotlin.random.Random


class DihedralAngleViewer private constructor(
	val dihedral: DihedralAngle,
	val minDegrees: Float,
	val maxDegrees: Float,
	val view: MoleculeRenderView,
	val mapsToCopy: MoleculeMaps
) : MotionViewer {

	companion object {

		fun make(desc: DihedralAngle.ConfDescription, assignmentInfo: Assignments.AssignmentInfo, view: MoleculeRenderView): DihedralAngleViewer {
			return DihedralAngleViewer(
				desc.make(assignmentInfo.molInfo.assignedMol, assignmentInfo.confSwitcher.atomResolverOrThrow),
				desc.minDegrees.toFloat(),
				desc.maxDegrees.toFloat(),
				view,
				assignmentInfo.confSwitcher.molInfo.maps
			)
		}

		fun make(desc: DihedralAngle.MolDescription, molCopy: Molecule, molMaps: MoleculeMaps, view: MoleculeRenderView): DihedralAngleViewer {
			return DihedralAngleViewer(
				desc.copyTo(molCopy, molMaps.atoms).make(),
				desc.minDegrees.toFloat(),
				desc.maxDegrees.toFloat(),
				view,
				molMaps
			)
		}
	}

	val radius = (maxDegrees - minDegrees)/2f

	// start in the center of the interval
	val pValue = Ref.of((minDegrees + maxDegrees)/2f)

	override val label =
		"Dihedral Angle: ${listOf(dihedral.a,  dihedral.b, dihedral.c, dihedral.d).joinToString(", ") { it.name }}"

	private var renderEffects = view.renderEffects.writer().apply {

		// show the rotatable atoms
		for (atom in dihedral.rotatedAtoms) {
			this[atom] = rotatableEffect
		}

		// show the dihedral atoms
		this[dihedral.a] = selectedEffect
		this[dihedral.b] = selectedEffect
		this[dihedral.c] = selectedEffect
		this[dihedral.d] = selectedEffect
	}

	override fun gui(imgui: Commands) = imgui.run {

		// show a slider to manipulate the dihedral angle
		text("Radius: %.1f degrees".format(radius))
		text("Range: %.1f to %.1f degrees".format(minDegrees, maxDegrees))
		if (minDegrees.isFinite() && maxDegrees.isFinite() && pValue.value.isFinite()) {
			if (sliderFloat("Angle", pValue, minDegrees, maxDegrees, format = "%.1f")) {
				dihedral.setDegrees(pValue.value.toDouble())
				view.moleculeChanged()
			}
		} else {
			text("(invalid angle)")
		}
	}

	override fun jiggle(rand: Random) {
		pValue.value = rand.nextFloatIn(minDegrees, maxDegrees)
		dihedral.setDegrees(pValue.value.toDouble())
		view.moleculeChanged()
	}

	override fun reset() {

		// remove any render effects
		renderEffects.close()
	}

	override fun mapAtomToOriginalMol(atom: Atom): Atom = mapsToCopy.atoms.getAOrThrow(atom)
}

private val selectedEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Inset, RenderEffect.Flags.Outset),
	0u, 255u, 0u
)

private val rotatableEffect = RenderEffect(
	ByteFlags.of(RenderEffect.Flags.Highlight, RenderEffect.Flags.Outset),
	255u, 255u, 255u
)