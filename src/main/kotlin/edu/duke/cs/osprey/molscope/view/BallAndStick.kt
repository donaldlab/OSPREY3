package edu.duke.cs.osprey.molscope.view

import cuchaz.kludge.tools.*
import cuchaz.kludge.vulkan.putColor4Bytes
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.render.CylinderRenderable
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.SphereRenderable
import edu.duke.cs.osprey.molscope.render.put
import org.joml.AABBf
import java.nio.ByteBuffer


/**
 * views a molecule using the ball and stick convention
 */
class BallAndStick(
	mol: Molecule,
	selector: MoleculeSelector = MoleculeSelectors.all
): MoleculeRenderView {

	private val info = MoleculeRenderInfo(mol, selector) {
		updateBonds()
	}

	// delegate some interface members to the info
	// TODO: delegation here gets MUCH less verbose if we upgrade to kotlin 1.4
	// see: https://kotlinlang.org/docs/reference/delegated-properties.html#delegating-to-another-property
	override var currentMol: Molecule
		get() = info.mol
		set(value) { info.mol = value }
	override var selector: MoleculeSelector
		get() = info.selector
		set(value) { info.selector = value }
	override fun moleculeChanged() = info.moleculeChanged()
	override val renderEffects: MoleculeRenderEffects
		get() = info.renderEffects
	override var isVisible: Boolean
		get() = info.isVisible
		set(value) { info.isVisible = value }

	override val molStack = MoleculeRenderStack(mol, info)


	// copy all the bonds (in the selection) as a list
	data class Bond(val i1: Int, val i2: Int)
	private val bonds = HashSet<Bond>()

	private fun updateBonds() {
		bonds.clear()
		info.selectedAtoms.forEachIndexed { i1, a1 ->
			for (a2 in currentMol.bonds.bondedAtoms(a1)) {
				val i2 = info.selectedAtoms.indexOf(a2)
				if (i2 < 0) {
					continue
				}

				// sort the indices to normalize the bond
				bonds.add(if (i1 < i2) {
					Bond(i1, i2)
				} else {
					Bond(i2, i1)
				})
			}
		}
	}

	// initialize the bonds
	init {
		updateBonds()
	}


	// render the atoms as spheres
	override val spheres = object : SphereRenderable {

		override val numVertices get() = info.selectedAtoms.size*4
		override val verticesSequence get() = info.sequence + renderEffects.sequence

		override fun fillVertexBuffer(buf: ByteBuffer, colorsMode: ColorsMode) {

			info.selectedAtoms.forEachIndexed { atomIndex, atom ->

				// write all vertex data 4 times
				for (i in 0 until 4) {

					// downgrade atom pos to floats for rendering
					buf.putFloat(atom.pos.x().toFloat())
					buf.putFloat(atom.pos.y().toFloat())
					buf.putFloat(atom.pos.z().toFloat())

					buf.putFloat(atomRadius)
					ElementProps[atom].apply {
						buf.putColor4Bytes(color[colorsMode])
					}

					// TODO: allow different indexing strategies (eg residue, molecule)
					buf.putInt(atomIndex)
					buf.put(renderEffects[atom])
				}
			}
		}

		override val boundingBox get() = calcBoundingBox()

		override val numOccluders get() = info.selectedAtoms.size

		override fun fillOcclusionBuffer(buf: ByteBuffer) {

			for (atom in info.selectedAtoms) {

				// downgrade atom pos to floats for rendering
				buf.putFloat(atom.pos.x.toFloat())
				buf.putFloat(atom.pos.y.toFloat())
				buf.putFloat(atom.pos.z.toFloat())
				buf.putFloat(atomRadius)
			}
		}
	}

	// render the bonds as cylinders
	override val cylinders = object : CylinderRenderable {

		override val numVertices get() = bonds.size*4
		override val verticesSequence get() = info.sequence + renderEffects.sequence

		override fun fillVertexBuffer(buf: ByteBuffer, colorsMode: ColorsMode) {
			for (bond in bonds) {
				val atomIndices = listOf(bond.i1, bond.i2)

				// write all vertex data 4 times
				for (i in 0 until 4) {

					for (atomIndex in atomIndices) {
						val atom = info.selectedAtoms[atomIndex]

						// downgrade atom pos to floats for rendering
						buf.putFloat(atom.pos.x().toFloat())
						buf.putFloat(atom.pos.y().toFloat())
						buf.putFloat(atom.pos.z().toFloat())
					}

					for (atomIndex in atomIndices) {
						buf.putFloat(bondRadius)
					}

					for (atomIndex in atomIndices) {
						val atom = info.selectedAtoms[atomIndex]

						buf.putColor4Bytes(ElementProps[atom].color[colorsMode])
					}

					for (atomIndex in atomIndices) {
						// TODO: allow different indexing strategies (eg residue, molecule)
						buf.putInt(atomIndex)
					}

					for (atomIndex in atomIndices) {
						val atom = info.selectedAtoms[atomIndex]

						buf.put(renderEffects[atom])
					}
				}
			}
		}

		override val boundingBox get() = calcBoundingBox()

		override val numOccluders get() = bonds.size

		override fun fillOcclusionBuffer(buf: ByteBuffer) {
			for (bond in bonds) {
				val atom1 = info.selectedAtoms[bond.i1]
				val atom2 = info.selectedAtoms[bond.i2]

				// downgrade atom pos to floats for rendering
				buf.putFloat(atom1.pos.x().toFloat())
				buf.putFloat(atom1.pos.y().toFloat())
				buf.putFloat(atom1.pos.z().toFloat())
				buf.putFloat(0f) // padding
				buf.putFloat(atom2.pos.x().toFloat())
				buf.putFloat(atom2.pos.y().toFloat())
				buf.putFloat(atom2.pos.z().toFloat())
				buf.putFloat(bondRadius)
			}
		}
	}

	override fun calcBoundingBox() =
		AABBf().apply {

			val r = atomRadius

			info.selectedAtoms.forEachIndexed { i, atom ->

				val x = atom.pos.x.toFloat()
				val y = atom.pos.y.toFloat()
				val z = atom.pos.z.toFloat()

				if (i == 0) {
					setMin(x - r, y - r, z - r)
					setMax(x + r, y + r, z + r)
				} else {
					expandToInclude(x - r, y - r, z - r)
					expandToInclude(x + r, y + r, z + r)
				}
			}
		}


	override fun getIndexed(index: Int) = info.selectedAtoms.getOrNull(index)
	// TODO: allow indexing other things?

	companion object {

		// TODO: make these parameters
		private const val atomRadius = 0.2f
		private const val bondRadius = 0.2f
	}

	// TODO: allow overriding these in constructor args
	private data class ElementProps(
		val color: Color
	) {

		companion object {

			operator fun get(atom: Atom) =
				when (atom.element) {
					Element.Hydrogen -> ElementProps(ColorPalette.lightGrey)
					Element.Carbon -> ElementProps(ColorPalette.darkGrey)
					Element.Nitrogen -> ElementProps(ColorPalette.blue)
					Element.Oxygen -> ElementProps(ColorPalette.red)
					Element.Sulfur -> ElementProps(ColorPalette.yellow)
					Element.Phosphorus -> ElementProps(ColorPalette.orange)
					else -> ElementProps(ColorPalette.darkGrey)
				}
		}
	}
}
