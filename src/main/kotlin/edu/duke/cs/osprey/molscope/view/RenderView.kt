package edu.duke.cs.osprey.molscope.view

import cuchaz.kludge.tools.indexOfOrNull
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.MoleculeSelector
import edu.duke.cs.osprey.molscope.molecule.MoleculeSelectors
import edu.duke.cs.osprey.molscope.render.CylinderRenderable
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.SphereRenderable
import org.joml.AABBf


/**
 * Represents a renderable view of a Thing.
 */
interface RenderView {
	var isVisible: Boolean
	fun calcBoundingBox(): AABBf? = null
	fun getIndexed(index: Int): Any? = null
	val spheres: SphereRenderable? get() = null
	val cylinders: CylinderRenderable? get() = null
}


/**
 * Represents a renderable view of a Molecule.
 */
interface MoleculeRenderView : RenderView {

	/** The currently-visible molecule, can be overridden by the molecule stack */
	var currentMol: Molecule

	var selector: MoleculeSelector

	/**
	 * Call to update the render view after making changes to the molecule.
	 * If the original molecule was updated, but the molecule stack has overrides, then the changes will not be visible.
	 * If the top molecule on the stack was updated, then the changes will be visible.
	 */
	fun moleculeChanged()

	val renderEffects: MoleculeRenderEffects
	val molStack: MoleculeRenderStack
}


/**
 * Tracks changes to molecules so renderers can update their render buffers.
 */
class MoleculeRenderInfo(
	mol: Molecule,
	selector: MoleculeSelector = MoleculeSelectors.all,
	val onChange: (Int) -> Unit = {}
) {

	var mol: Molecule = mol
		set(value) {
			field = value
			renderEffects.mol = value
			selectedAtoms = selector(value)
			changed()
		}

	var selector: MoleculeSelector = selector
		set(value) {
			field = value
			selectedAtoms = value(mol)
			changed()
		}

	fun moleculeChanged() {
		selectedAtoms = selector(mol)
		changed()
	}

	var sequence = 0
		private set

	private fun changed() {
		sequence += 1
		onChange(sequence)
	}

	var selectedAtoms = selector(mol)
		private set

	var renderEffects = MoleculeRenderEffects(mol)
		private set

	var isVisible = true
}


/**
 * A mechanism to temporarily override a molecule being rendered
 */
class MoleculeRenderStack(
	val originalMol: Molecule,
	val info: MoleculeRenderInfo
) {
	private val overrideMols = ArrayList<Molecule>()

	inner class StackedMol(mol: Molecule) {

		var mol: Molecule = mol
			private set

		fun replace(mol: Molecule) {
			val i = overrideMols.indexOfOrNull(this.mol)
				?: throw NoSuchElementException("molecule was already popped from the stack")
			overrideMols[i] = mol
			this.mol = mol
			info.mol = mol
		}

		fun pop() {
			overrideMols.remove(mol)
			info.mol = overrideMols.lastOrNull() ?: originalMol
		}
	}

	fun push(mol: Molecule): StackedMol {
		overrideMols.add(mol)
		info.mol = mol
		return StackedMol(mol)
	}

	val size get() = overrideMols.size
}
