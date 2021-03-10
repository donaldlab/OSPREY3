package edu.duke.cs.osprey.molscope.view

import cuchaz.kludge.tools.*
import cuchaz.kludge.vulkan.putColor4Bytes
import edu.duke.cs.osprey.molscope.molecule.*
import edu.duke.cs.osprey.molscope.render.MoleculeRenderEffects
import edu.duke.cs.osprey.molscope.render.SphereRenderable
import edu.duke.cs.osprey.molscope.render.put
import org.joml.AABBf
import java.nio.ByteBuffer


/**
 * views a molecule using the space-filling sphere convention
 */
// TODO: optimize molecule transformations so we don't have to re-create the whole view for a large molecule?
class SpaceFilling(
	mol: Molecule,
	selector: MoleculeSelector = MoleculeSelectors.all
): MoleculeRenderView {

	private val info = MoleculeRenderInfo(mol, selector)

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


	override val spheres = object : SphereRenderable {
		
		override val numVertices get() = info.selectedAtoms.size*4
		override val verticesSequence get() = info.sequence + renderEffects.sequence

		override fun fillVertexBuffer(buf: ByteBuffer, colorsMode: ColorsMode) {

			info.selectedAtoms.forEachIndexed { atomIndex, atom ->

				// write all vertex data 4 times
				for (i in 0 until 4) {

					// downgrade atom pos to floats for rendering
					buf.putFloat(atom.pos.x.toFloat())
					buf.putFloat(atom.pos.y.toFloat())
					buf.putFloat(atom.pos.z.toFloat())

					ElementProps[atom].apply {
						buf.putFloat(radius)
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

				ElementProps[atom].apply {
					buf.putFloat(radius)
				}
			}
		}
	}

	override fun calcBoundingBox() =
		AABBf().apply {
			info.selectedAtoms.forEachIndexed { i, atom ->

				val x = atom.pos.x.toFloat()
				val y = atom.pos.y.toFloat()
				val z = atom.pos.z.toFloat()
				val r = ElementProps[atom].radius

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


	// TODO: allow overriding these in constructor args
	private data class ElementProps(
		val radius: Float,
		val color: Color
	) {

		companion object {

			operator fun get(atom: Atom) =
				when (atom.element) {
					Element.Hydrogen -> ElementProps(1f, ColorPalette.lightGrey)
					Element.Carbon -> ElementProps(1.75f, ColorPalette.darkGrey)
					Element.Nitrogen -> ElementProps(1.55f, ColorPalette.blue)
					Element.Oxygen -> ElementProps(1.4f, ColorPalette.red)
					Element.Sulfur -> ElementProps(1.8f, ColorPalette.yellow)
					Element.Phosphorus -> ElementProps(1.8f, ColorPalette.orange)
					else -> ElementProps(2.0f, ColorPalette.darkGrey)
				}
		}
	}
}
