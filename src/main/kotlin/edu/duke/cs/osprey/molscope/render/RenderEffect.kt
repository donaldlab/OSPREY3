package edu.duke.cs.osprey.molscope.render

import cuchaz.kludge.tools.ByteFlags
import cuchaz.kludge.tools.put
import cuchaz.kludge.vulkan.ClearValue
import edu.duke.cs.osprey.molscope.molecule.Atom
import edu.duke.cs.osprey.molscope.molecule.Molecule
import edu.duke.cs.osprey.molscope.molecule.MoleculeSelector
import java.nio.ByteBuffer
import java.util.*


data class RenderEffect(val flags: ByteFlags<Flags>, val r: UByte, val g: UByte, val b: UByte) {

	enum class Flags(override val value: Byte) : ByteFlags.Bit {

		// these must match the render effects in post.frag

		/**
		 * Shows the selected items in brighter colors:
		 * color.rgb *= (1 + rgb)
		 */
		Highlight(1 shl 0),

		/**
		 * Draws an inset border around the selected items:
		 * color.rgb = rgb
		 */
		Inset(1 shl 1),

		/**
		 * Draws an outset border around the selected items:
		 * color.rgb = rgb
		 */
		Outset(1 shl 2),
	}

	companion object {
		val clearColor = ClearValue.Color.Int(0, 0, 0, 0)
	}
}

fun ByteBuffer.put(effect: RenderEffect?) {
	if (effect != null) {
		put(effect.r)
		put(effect.g)
		put(effect.b)
		put(effect.flags.value)
	} else {
		putInt(0)
	}
}

class MoleculeRenderEffects(mol: Molecule) {

	var mol: Molecule = mol
		set(value) {
			field = value
			sequence += 1
			for (writer in writers) {
				writer.clear()
			}
		}

	inner class Writer : AutoCloseable {

		val mol get() = this@MoleculeRenderEffects.mol

		init {
			writers.add(this)
		}

		override fun close() {
			writers.remove(this)
			sequence += 1
		}

		private val effectsByAtom = IdentityHashMap<Atom,RenderEffect>()

		fun clear() {
			if (effectsByAtom.isNotEmpty()) {
				effectsByAtom.clear()
				sequence += 1
			}
		}

		operator fun set(atom: Atom, effect: RenderEffect) {
			if (effectsByAtom[atom] != effect) {
				effectsByAtom[atom] = effect
				sequence += 1
			}
		}

		operator fun set(atoms: Collection<Atom>, effect: RenderEffect) {
			var isChanged = false
			for (atom in atoms) {
				if (effectsByAtom[atom] != effect) {
					effectsByAtom[atom] = effect
					isChanged = true
				}
			}
			if (isChanged) {
				sequence += 1
			}
		}

		fun remove(atom: Atom) {
			val oldVal = effectsByAtom.remove(atom)
			if (oldVal != null) {
				sequence += 1
			}
		}

		fun remove(atoms: Collection<Atom>) {
			var isChanged = false
			for (atom in atoms) {
				val oldVal = effectsByAtom.remove(atom)
				if (oldVal != null) {
					isChanged = true
				}
			}
			if (isChanged) {
				sequence += 1
			}
		}

		operator fun set(selector: MoleculeSelector, effect: RenderEffect) {
			set(selector(mol), effect)
		}

		fun remove(selector: MoleculeSelector) {
			remove(selector(mol))
		}

		operator fun get(atom: Atom): RenderEffect? = effectsByAtom[atom]
	}

	private val writers = ArrayList<Writer>()

	fun writer() = Writer()

	var sequence: Int = 0
		private set

	/**
	 * Gets the effect from the most recently-created writer that has set an effect
	 */
	operator fun get(atom: Atom): RenderEffect? {
		for (i in (0 until writers.size).reversed()) {
			val effect = writers[i][atom]
			if (effect != null) {
				return effect
			}
		}
		return null
	}
}
