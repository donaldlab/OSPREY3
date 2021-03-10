package edu.duke.cs.osprey.molscope.render

import cuchaz.kludge.tools.AutoCloser
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.SIZE_BYTES
import cuchaz.kludge.tools.diff
import cuchaz.kludge.vulkan.*
import edu.duke.cs.osprey.molscope.shaders.Shaders
import edu.duke.cs.osprey.molscope.view.ColorsMode
import org.joml.AABBf
import java.nio.ByteBuffer
import java.util.*
import kotlin.math.max


internal class CylinderRenderer(
	val slideRenderer: SlideRenderer
): AutoCloseable {

	private val closer = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose(replace: R? = null) = apply { closer.add(this, replace) }
	private fun <R> R.autoClose(block: R.() -> Unit): R = apply { closer.add { block() } }
	override fun close() = closer.close()

	private val device get() = slideRenderer.device

	// make the graphics pipeline
	val graphicsPipeline = slideRenderer
		.graphicsPipeline(
			listOf(
				device.shaderModule(Shaders["cylinder.vert"])
					.autoClose()
					.stage("main", ShaderStage.Vertex),
				device.shaderModule(Shaders["cylinder.frag"])
					.autoClose()
					.stage("main", ShaderStage.Fragment)
			),
			vertexInput = VertexInput {
				binding(stride = Float.SIZE_BYTES*8 + Byte.SIZE_BYTES*8 + Int.SIZE_BYTES*4) {
					// pos
					attribute(
						location = 0,
						format = Image.Format.R32G32B32_SFLOAT,
						offset = 0
					)
					attribute(
						location = 1,
						format = Image.Format.R32G32B32_SFLOAT,
						offset = Float.SIZE_BYTES*3
					)
					// radius
					attribute(
						location = 2,
						format = Image.Format.R32_SFLOAT,
						offset = Float.SIZE_BYTES*6
					)
					attribute(
						location = 3,
						format = Image.Format.R32_SFLOAT,
						offset = Float.SIZE_BYTES*7
					)
					// color
					attribute(
						location = 4,
						format = Image.Format.R8G8B8A8_UNORM,
						offset = Float.SIZE_BYTES*8
					)
					attribute(
						location = 5,
						format = Image.Format.R8G8B8A8_UNORM,
						offset = Float.SIZE_BYTES*8 + Byte.SIZE_BYTES*4
					)
					// target index
					attribute(
						location = 6,
						format = Image.Format.R32_SINT,
						offset = Float.SIZE_BYTES*8 + Byte.SIZE_BYTES*8
					)
					attribute(
						location = 7,
						format = Image.Format.R32_SINT,
						offset = Float.SIZE_BYTES*8 + Byte.SIZE_BYTES*8 + Int.SIZE_BYTES
					)
					// effects
					attribute(
						location = 8,
						format = Image.Format.R8G8B8A8_UINT,
						offset = Float.SIZE_BYTES*8 + Byte.SIZE_BYTES*8 + Int.SIZE_BYTES*2
					)
					attribute(
						location = 9,
						format = Image.Format.R8G8B8A8_UINT,
						offset = Float.SIZE_BYTES*8 + Byte.SIZE_BYTES*8 + Int.SIZE_BYTES*3
					)
				}
			},
			inputAssembly = InputAssembly(InputAssembly.Topology.TriangleStrip)
		)
		.autoClose()

	inner class Entry(val src: CylinderRenderable): AutoCloseable {

		private val closer = AutoCloser()
		private fun <R:AutoCloseable> R.autoClose(replace: R? = null) = apply { closer.add(this, replace) }
		override fun close() = closer.close()

		// Vulkan won't allow 0-sized buffers, so use at least one byte
		private fun bufSize(size: Long) = max(1L, size)

		inner class VBO(val numVertices: Int) : AutoCloseable {

			val buf = device
				.buffer(
					size = bufSize(numVertices*graphicsPipeline.vertexInput.size),
					usage = IntFlags.of(Buffer.Usage.VertexBuffer, Buffer.Usage.TransferDst)
				)
			val allocated = buf.allocateDevice()

			override fun close() {
				buf.close()
				allocated.close()
			}
		}
		var vertexBuf = VBO(src.numVertices).autoClose()

		private val dirtyness = Dirtyness()

		fun update(colorsMode: ColorsMode) {

			// track state changes
			dirtyness.update(colorsMode, src.verticesSequence)
			if (!dirtyness.isDirty) {
				return
			}

			// reallocate buffers if needed
			if (vertexBuf.numVertices < src.numVertices) {
				vertexBuf = VBO(src.numVertices).autoClose(replace = vertexBuf)
			}

			// update buffers
			vertexBuf.allocated.transferHtoD { buf ->
				src.fillVertexBuffer(buf, colorsMode)
				buf.flip()
			}
		}
	}

	private val entries = IdentityHashMap<CylinderRenderable,Entry>()
		.autoClose {
			for (entry in values) {
				entry.close()
			}
			clear()
		}

	fun update(sources: List<CylinderRenderable>) {

		// add/remove entries
		sources.diff(
			entries,
			added = { src ->
				entries[src] = Entry(src)
			},
			removed = { src, entry ->
				entry.close()
				entries.remove(src)
			}
		)

		// update entries
		for (entry in entries.values) {
			entry.update(ColorsMode.current)
		}
	}

	fun render(cmdbuf: CommandBuffer, src: CylinderRenderable, viewIndex: Int) = cmdbuf.apply {

		if (src.numVertices > 0) {

			val entry = entries[src] ?: throw NoSuchElementException("call update() with this source, before render()")

			// draw geometry
			bindPipeline(graphicsPipeline)
			bindDescriptorSet(slideRenderer.mainDescriptorSet, graphicsPipeline)
			bindVertexBuffer(entry.vertexBuf.buf)
			pushConstants(graphicsPipeline, IntFlags.of(ShaderStage.Fragment),
				viewIndex, 0, 0, 0
			)
			draw(vertices = src.numVertices)
		}
	}
}


interface CylinderRenderable {

	val numVertices: Int
	val verticesSequence: Int
	fun fillVertexBuffer(buf: ByteBuffer, colorsMode: ColorsMode)

	val boundingBox: AABBf? get() = null
	val numOccluders: Int get() = 0
	fun fillOcclusionBuffer(buf: ByteBuffer) {}
}
