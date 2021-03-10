package edu.duke.cs.osprey.molscope.render

import cuchaz.kludge.imgui.Imgui
import cuchaz.kludge.tools.AutoCloser
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.position
import cuchaz.kludge.vulkan.*
import java.awt.image.BufferedImage
import javax.imageio.ImageIO


class LoadedImage internal constructor(
	internal val queue: Queue,
	/** 4-byte pixels, tightly-packed, top row first */
	rgba: ByteArray,
	val width: Int,
	val height: Int
) : AutoCloseable {

	internal constructor(queue: Queue, image: BufferedImage) : this(
		queue,
		image.toRGBA(),
		image.width,
		image.height
	)

	private val closer = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose(replace: R? = null) = apply { closer.add(this, replace) }
	override fun close() = closer.close()

	private val commandPool = queue.device
		.commandPool(queue.family)
		.autoClose()

	// allocate the image on the GPU
	private val gpuImage = queue.device.
		image(
			type = Image.Type.TwoD,
			extent = Extent3D(width, height, 1),
			format = Image.Format.R8G8B8A8_UNORM,
			usage = IntFlags.of(Image.Usage.TransferDst, Image.Usage.Sampled),
			tiling = Image.Tiling.Optimal
			// TODO: Optimal tiling is better for rendering,
			//  but the HtoD transfer doesn't re-tile from Linear to Optimal
			//  so would need to find a way to do that somehow to support Optimal tiling
			//  also, OSX seems to work with Optimal tiling correctly somehow, but not linear tiling
		)
		.autoClose()
		.allocateDevice()
		.autoClose()
		.apply {

			// upload image to the GPU
			transitionImage(IntFlags.of(Access.HostWrite), Image.Layout.General)
			transferHtoD { buf ->

				// upload the image to the GPU,
				// but respect the row pitch the GPU wants

				val srcPitch = width*4
				val dstPitch = image.getSubresourceLayout().rowPitch.toInt()

				for (y in 0 until height) {

					// copy the row
					buf.put(
						rgba,
						y*srcPitch,
						srcPitch
					)

					// skip to the next row
					buf.position += dstPitch - srcPitch
				}
				buf.flip()
			}

			// prep the image for the fragment shader
			transitionImage(IntFlags.of(Access.ShaderRead), Image.Layout.ShaderReadOnlyOptimal)
		}

	private fun Image.Allocated.transitionImage(access: IntFlags<Access>, layout: Image.Layout) {
		queue.submit(commandPool.buffer().apply {
			begin(IntFlags.of(CommandBuffer.Usage.OneTimeSubmit))

			pipelineBarrier(
				srcStage = IntFlags.of(PipelineStage.AllCommands),
				dstStage = IntFlags.of(PipelineStage.AllCommands),
				images = listOf(
					image.barrier(
						dstAccess = access,
						newLayout = layout
					)
				)
			)

			end()
		})
	}

	val view = gpuImage.image
		.view()
		.autoClose()

	// make a sampler
	val sampler = queue.device
		.sampler()
		.autoClose()

	val descriptor: Imgui.ImageDescriptor by lazy {
		Imgui.imageDescriptor(view, sampler)
			.autoClose()
	}
}

/**
 * Convert the image bytes into an ImageIO Buffer
 */
fun ByteArray.toBuffer(): BufferedImage =
	inputStream().use { stream ->
		ImageIO.read(stream)
	}

/**
 * Convert an ImageIO buffer to a tightly-packed RGBA array.
 *
 * As far as image-processing functions go, this will be super duper slow.
 * Hopefully, that won't be a bottleneck for us though.
 */
fun BufferedImage.toRGBA(): ByteArray {

	val bytes = ByteArray(width*height*4)
	var i = 0

	for (y in 0 until height) {
		for (x in 0 until width) {
			val argb = getRGB(x, y)
			val a = (argb shr 8*3) and 0xff
			val r = (argb shr 8*2) and 0xff
			val g = (argb shr 8*1) and 0xff
			val b = (argb shr 8*0) and 0xff
			bytes[i++] = r.toByte()
			bytes[i++] = g.toByte()
			bytes[i++] = b.toByte()
			bytes[i++] = a.toByte()
		}
	}

	return bytes
}
