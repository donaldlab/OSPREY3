package edu.duke.cs.osprey.molscope.render

import cuchaz.kludge.tools.*
import cuchaz.kludge.vulkan.*
import cuchaz.kludge.vulkan.Queue
import edu.duke.cs.osprey.molscope.shaders.Shaders
import edu.duke.cs.osprey.molscope.tools.SphereGrid
import org.joml.AABBf
import kotlin.math.max
import kotlin.math.min


/**
 * An implementation of world-space ambient occlusion for lighting
 *
 * This technique is based loosely on the "Voxel Field Ambient Occlusion" technique presented in:
 * https://lambdacube3d.wordpress.com/2016/05/15/ambient-occlusion-fields/
 * https://www.reddit.com/r/gamedev/comments/4jhqot/ambient_occlusion_fields/
 */
internal class OcclusionCalculator(
	val queue: Queue
) : AutoCloseable {

	val device: Device get() = queue.device

	private val closer = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose() = apply { closer.add(this) }
	override fun close() = closer.close()

	private val descriptorPool = device.descriptorPool(
		maxSets = 3,
		sizes = DescriptorType.Counts(
			DescriptorType.StorageBuffer to 6,
			DescriptorType.StorageImage to 3,
			DescriptorType.CombinedImageSampler to 1,
			DescriptorType.UniformBuffer to 1
		)
	).autoClose()

	// make the compute pipeline
	private val linesBinding = DescriptorSetLayout.Binding(
		binding = 0,
		type = DescriptorType.StorageBuffer,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val spheresBinding = DescriptorSetLayout.Binding(
		binding = 1,
		type = DescriptorType.StorageBuffer,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val cylindersBinding = DescriptorSetLayout.Binding(
		binding = 2,
		type = DescriptorType.StorageBuffer,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val occlusionBinding = DescriptorSetLayout.Binding(
		binding = 3,
		type = DescriptorType.StorageImage,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val descriptorSetLayout = device.descriptorSetLayout(listOf(
		linesBinding, spheresBinding, cylindersBinding, occlusionBinding
	)).autoClose()
	private val descriptorSet = descriptorPool.allocate(descriptorSetLayout)
	private val pipeline = device
		.computePipeline(
			stage = device.shaderModule(Shaders["ambientOcclusion.comp"])
				.autoClose()
				.stage("main", ShaderStage.Compute),
			descriptorSetLayouts = listOf(descriptorSetLayout),
			pushConstantRanges = listOf(
				PushConstantRange(IntFlags.of(ShaderStage.Compute), 16*2)
			)
		).autoClose()

	// make the blur pipeline
	private val occlusionInBinding = DescriptorSetLayout.Binding(
		binding = 0,
		type = DescriptorType.StorageImage,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val occlusionOutBinding = DescriptorSetLayout.Binding(
		binding = 1,
		type = DescriptorType.StorageImage,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val blurDescriptorSetLayout = device.descriptorSetLayout(listOf(
		occlusionInBinding, occlusionOutBinding
	)).autoClose()
	private val blurDescriptorSet = descriptorPool.allocate(blurDescriptorSetLayout)
	private val blurPipeline = device
		.computePipeline(
			stage = device.shaderModule(Shaders["occlusionBlur.comp"])
				.autoClose()
				.stage("main", ShaderStage.Compute),
			descriptorSetLayouts = listOf(blurDescriptorSetLayout),
			pushConstantRanges = listOf(
				PushConstantRange(IntFlags.of(ShaderStage.Compute), 16*3)
			)
		).autoClose()

	// make the min pipeline
	private val minOcclusionImageBinding = DescriptorSetLayout.Binding(
		binding = 1,
		type = DescriptorType.CombinedImageSampler,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val minBoundsBinding = DescriptorSetLayout.Binding(
		binding = 2,
		type = DescriptorType.UniformBuffer,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val minSpheresBinding = DescriptorSetLayout.Binding(
		binding = 4,
		type = DescriptorType.StorageBuffer,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val minCylindersBinding = DescriptorSetLayout.Binding(
		binding = 5,
		type = DescriptorType.StorageBuffer,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val minOutBinding = DescriptorSetLayout.Binding(
		binding = 6,
		type = DescriptorType.StorageBuffer,
		stages = IntFlags.of(ShaderStage.Compute)
	)
	private val minDescriptorSetLayout = device.descriptorSetLayout(listOf(
		minOcclusionImageBinding, minBoundsBinding, minSpheresBinding, minCylindersBinding, minOutBinding
	)).autoClose()
	private val minDescriptorSet = descriptorPool.allocate(minDescriptorSetLayout)
	private val minPipeline = device
		.computePipeline(
			stage = device.shaderModule(Shaders["occlusionMin.comp"])
				.autoClose()
				.stage("main", ShaderStage.Compute),
			descriptorSetLayouts = listOf(minDescriptorSetLayout)
		).autoClose()

	// make a command buffer
	private val commandPool = device
		.commandPool(
			queue.family,
			flags = IntFlags.of(CommandPool.Create.ResetCommandBuffer)
		)
		.autoClose()

	private fun cmdbuf(block: CommandBuffer.() -> Unit) {
		queue.submit(commandPool.buffer().apply {
			begin(IntFlags.of(CommandBuffer.Usage.OneTimeSubmit))
			this.block()
			end()
		})
	}

	// make an interpolating sampler for the occlusion image
	private val sampler = device
		.sampler(
			minFilter = Sampler.Filter.Linear,
			magFilter = Sampler.Filter.Linear,
			addressU = Sampler.Address.ClampToEdge,
			addressV = Sampler.Address.ClampToEdge,
			addressW = Sampler.Address.ClampToEdge
		)
		.autoClose()


	inner class Field(
		val extent: Extent3D,
		val gridSubdivisions: Int,
		val renderables: ViewRenderables
	) : AutoCloseable {

		private val closer = AutoCloser()
		private fun <R:AutoCloseable> R.autoClose() = apply { closer.add(this) }
		override fun close() = closer.close()

		// allocate an image for the final blurred occlusion data
		val blurredOcclusionImage = device
			.image(
				type = Image.Type.ThreeD,
				extent = extent,
				format = Image.Format.R32_SINT,
				usage = IntFlags.of(Image.Usage.Storage, Image.Usage.Sampled)
			)
			.autoClose()
			.allocateDevice()
			.autoClose()
		val blurredOcclusionView = blurredOcclusionImage.image.view().autoClose()

		// use a sphere grid to define which rays to shoot
		val sphereGrid = SphereGrid(gridSubdivisions)
			.filter {
				// keep only the directions facing the "front"
				// so we ignore antipodal directions (the sphere grid has some symmetry)
				if (it.z == 0.0) {
					if (it.y == 0.0) {
						it.x >= 0.0
					} else {
						it.y >= 0.0
					}
				} else {
					it.z >= 0.0
				}
			}
			.map { it.toFloat() }

		// each line supports two rays
		val maxOcclusion = sphereGrid.size*2

		// allocate the lines buffer from the sphere grid
		val linesBuf = device
			.buffer(
				size = 16L + sphereGrid.size*16L, // sizeof(ivec4), sizeof(vec4)
				usage = IntFlags.of(Buffer.Usage.StorageBuffer, Buffer.Usage.TransferDst)
			)
			.autoClose()
			.allocateDevice()
			.autoClose()
			.apply {
				transferHtoD { buf ->

					// write the header
					buf.putInts(extent.width, extent.height, extent.depth)
					buf.putInt(sphereGrid.size)

					// write the lines
					for (v in sphereGrid) {
						buf.putFloats(v.x, v.y, v.z)
						buf.putFloat(0f) // padding
					}
				}
			}

		// compute the bounding box for the occlusion field
		val box = (
				renderables.spheres.mapNotNull { it.boundingBox }
				+ renderables.cylinders.mapNotNull { it.boundingBox }
			)
			.run {
				// if there's nothing to render, just use a dummy box
				if (isEmpty()) {
					listOf(AABBf())
				} else {
					this
				}
			}
			.reduce { a, b -> AABBf(a).union(b) }
			.apply {

				// pad the box a little to give us some breathing room
				expand(0.1f)
			}

		// Vulkan won't allow 0-sized buffers, so use at least one byte
		private fun bufSize(size: Long) = max(1L, size)

		// upload the spheres
		val numSpheres = renderables.spheres.sumBy { it.numOccluders }
		val spheresBuf = device
			.buffer(
				size = bufSize(16L + numSpheres*16L), // sizeof(struct Sphere)
				usage = IntFlags.of(Buffer.Usage.StorageBuffer, Buffer.Usage.TransferDst)
			)
			.autoClose()
			.allocateDevice()
			.autoClose()
			.apply {
				transferHtoD { buf ->

					buf.putInt(numSpheres)
					buf.putInts(0, 0, 0) // padding

					// fill the buffer
					renderables.spheres.forEach { it.fillOcclusionBuffer(buf) }
					buf.flip()
				}
			}

		// upload the cylinders
		val numCylinders = renderables.cylinders.sumBy { it.numOccluders }
		val cylindersBuf = device
			.buffer(
				size = bufSize(16L + numCylinders*32L), // sizeof(struct Cylinder)
				usage = IntFlags.of(Buffer.Usage.StorageBuffer, Buffer.Usage.TransferDst)
			)
			.autoClose()
			.allocateDevice()
			.autoClose()
			.apply {
				transferHtoD { buf ->

					buf.putInt(numCylinders)
					buf.putInts(0, 0, 0) // padding

					// fill the buffer
					renderables.cylinders.forEach { it.fillOcclusionBuffer(buf) }
					buf.flip()
				}
			}

		// allocate the buffer for the min occlusion
		val minBuf = device
			.buffer(
				size = bufSize(Float.SIZE_BYTES*(numSpheres + numCylinders).toLong()),
				usage = IntFlags.of(Buffer.Usage.StorageBuffer, Buffer.Usage.TransferSrc)
			)
			.autoClose()
			.allocateDevice()
			.autoClose()

		// allocate the occlusion image
		val occlusionImage = device
			.image(
				type = blurredOcclusionImage.image.type,
				extent = extent,
				format = blurredOcclusionImage.image.format,
				usage = IntFlags.of(Image.Usage.Storage, Image.Usage.TransferDst)
			)
			.autoClose()
			.allocateDevice()
			.autoClose()
		val occlusionView = occlusionImage.image.view().autoClose()

		// upload the occlusion field info for fragment shaders
		val boundsBuf = device
			.buffer(
				size = Int.SIZE_BYTES*4L + Float.SIZE_BYTES*8L,
				usage = IntFlags.of(Buffer.Usage.UniformBuffer)
			)
			.autoClose()
			.allocateDevice()
			.autoClose()
			.apply {
				transferHtoD { buf ->
					buf.putInts(
						extent.width,
						extent.height,
						extent.depth
					)
					buf.putInt(maxOcclusion)
					buf.putFloats(box.minX, box.minY, box.minZ)
					buf.putFloat(1f) // initial min is 1
					buf.putFloats(box.maxX, box.maxY, box.maxZ)
					buf.putFloat(0f) // padding
					buf.flip()
				}
			}

		init {
			// update the descriptor sets
			device.updateDescriptorSets(
				writes = listOf(

					descriptorSet.address(linesBinding).write(
						DescriptorSet.BufferInfo(linesBuf.buffer)
					),
					descriptorSet.address(spheresBinding).write(
						DescriptorSet.BufferInfo(spheresBuf.buffer)
					),
					descriptorSet.address(cylindersBinding).write(
						DescriptorSet.BufferInfo(cylindersBuf.buffer)
					),
					descriptorSet.address(occlusionBinding).write(
						DescriptorSet.ImageInfo(null, occlusionView, Image.Layout.General)
					),

					blurDescriptorSet.address(occlusionInBinding).write(
						DescriptorSet.ImageInfo(null, occlusionView, Image.Layout.General)
					),
					blurDescriptorSet.address(occlusionOutBinding).write(
						DescriptorSet.ImageInfo(null, blurredOcclusionView, Image.Layout.General)
					),

					minDescriptorSet.address(minOcclusionImageBinding).write(
						DescriptorSet.ImageInfo(sampler, blurredOcclusionView, Image.Layout.General)
					),
					minDescriptorSet.address(minBoundsBinding).write(
						DescriptorSet.BufferInfo(boundsBuf.buffer)
					),
					minDescriptorSet.address(minSpheresBinding).write(
						DescriptorSet.BufferInfo(spheresBuf.buffer)
					),
					minDescriptorSet.address(minCylindersBinding).write(
						DescriptorSet.BufferInfo(cylindersBuf.buffer)
					),
					minDescriptorSet.address(minOutBinding).write(
						DescriptorSet.BufferInfo(minBuf.buffer)
					)
				)
			)

			// clear the occlusion field
			cmdbuf {
				pipelineBarrier(
					srcStage = IntFlags.of(PipelineStage.TopOfPipe),
					dstStage = IntFlags.of(PipelineStage.ComputeShader),
					images = listOf(
						occlusionImage.image.barrier(
							dstAccess = IntFlags.of(Access.ShaderRead, Access.ShaderWrite),
							newLayout = Image.Layout.TransferDstOptimal
						)
					)
				)
				clearImage(occlusionImage.image, Image.Layout.TransferDstOptimal, ClearValue.Color.Int(0, 0, 0, 0))
			}
			queue.waitForIdle()
		}

		fun updateDescriptorSet(slideRenderer: SlideRenderer) {

			slideRenderer.device.updateDescriptorSets(
				writes = listOf(
					slideRenderer.mainDescriptorSet.address(slideRenderer.occlusionImageBinding).write(
						DescriptorSet.ImageInfo(sampler, blurredOcclusionView, Image.Layout.ShaderReadOnlyOptimal)
					),
					slideRenderer.mainDescriptorSet.address(slideRenderer.boundsBinding).write(
						DescriptorSet.BufferInfo(boundsBuf.buffer)
					)
				)
			)
		}

		private var linesProcessed = 0

		var active: Boolean = false

		val processingProgress get() =
			if (active) {
				linesProcessed.toDouble()/sphereGrid.size.toDouble()
			} else {
				1.0
			}

		val needsProcessing get() =
			active && linesProcessed < sphereGrid.size

		fun process() {

			if (!needsProcessing) {
				return
			}

			// arbitrarily chosen, seems to work well on my laptop
			val frameCostBudget = 500

			// how many lines should we process this time?
			// (cylinders are a bit more expensive to process than spheres)
			val lineCost = 1 + numSpheres + numCylinders*2
			val numLines = max(1, frameCostBudget/lineCost)
			val iLines = linesProcessed until min(sphereGrid.size, linesProcessed + numLines)

			// call the compute shader
			cmdbuf {
				// prep the occlusion image
				pipelineBarrier(
					srcStage = IntFlags.of(PipelineStage.TopOfPipe),
					dstStage = IntFlags.of(PipelineStage.ComputeShader),
					images = listOf(
						occlusionImage.image.barrier(
							dstAccess = IntFlags.of(Access.ShaderWrite),
							newLayout = Image.Layout.General
						)
					)
				)

				memstack { mem ->

					// run the occlusion kernel
					bindPipeline(pipeline)
					bindDescriptorSet(descriptorSet, pipeline)
					pushConstants(pipeline, IntFlags.of(ShaderStage.Compute), mem.malloc(16*2).apply {
						putFloats(box.minX, box.minY, box.minZ)
						putInt(iLines.first)
						putFloats(box.maxX, box.maxY, box.maxZ)
						putInt(iLines.last)
						flip()
					})
					dispatch(
						extent.width,
						extent.height,
						extent.depth
						// no 4d kernels in Vulkan, so multiplex the z and line params to keep it 3d
					)
				}
			}

			linesProcessed += numLines

			// if all the processing is done, apply the post processing
			if (!needsProcessing) {

				cmdbuf {

					// prep the occlusion image
					pipelineBarrier(
						srcStage = IntFlags.of(PipelineStage.TopOfPipe),
						dstStage = IntFlags.of(PipelineStage.ComputeShader),
						images = listOf(
							occlusionImage.image.barrier(
								dstAccess = IntFlags.of(Access.ShaderRead),
								newLayout = Image.Layout.General
							),
							blurredOcclusionImage.image.barrier(
								dstAccess = IntFlags.of(Access.ShaderWrite),
								newLayout = Image.Layout.General
							)
						)
					)

					memstack { mem ->

						// run the blur kernel
						bindPipeline(blurPipeline)
						bindDescriptorSet(blurDescriptorSet, blurPipeline)
						pushConstants(blurPipeline, IntFlags.of(ShaderStage.Compute), mem.malloc(16*3).apply {
							putFloats(box.minX, box.minY, box.minZ)
							putFloat(0f) // padding
							putFloats(box.maxX, box.maxY, box.maxZ)
							putFloat(0f) // padding
							putInts(
								extent.width,
								extent.height,
								extent.depth
							)
							putInt(maxOcclusion)
							flip()
						})
						dispatch(extent)

						// run the min kernel
						bindPipeline(minPipeline)
						bindDescriptorSet(minDescriptorSet, minPipeline)
						dispatch(numSpheres + numCylinders)
					}
				}

				// wait for the kernels to finish, so we can read the min occlusion
				queue.waitForIdle()

				// read the min occlusion
				var minOcclusion = 1f
				minBuf.transferDtoH { buf ->
					for (i in 0 until numSpheres + numCylinders) {
						minOcclusion = min(minOcclusion, buf.float)
					}
				}

				// write the min occlusion to the bounds buf
				boundsBuf.transferHtoD { buf ->
					buf.position = 28
					buf.putFloat(minOcclusion)
				}
			}
		}
	}
}


/**
 * A renderer for the occlusion field
 */
internal class OcclusionRenderer(
	val slideRenderer: SlideRenderer
) : AutoCloseable {

	private val closer = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose(replace: R? = null) = apply { closer.add(this, replace) }
	override fun close() = closer.close()

	private val device get() = slideRenderer.device

	// make the render pipeline
	private val renderPipeline = slideRenderer
		.graphicsPipeline(
			listOf(
				device.shaderModule(Shaders["ambientOcclusion.vert"])
					.autoClose()
					.stage("main", ShaderStage.Vertex),
				device.shaderModule(Shaders["ambientOcclusion.geom"])
					.autoClose()
					.stage("main", ShaderStage.Geometry),
				device.shaderModule(Shaders["ambientOcclusion.frag"])
					.autoClose()
					.stage("main", ShaderStage.Fragment)
			),
			inputAssembly = InputAssembly(InputAssembly.Topology.PointList)
		)
		.autoClose()

	fun barriers(cmdbuf: CommandBuffer, field: OcclusionCalculator.Field) = cmdbuf.run {

		// prep the occlusion images
		pipelineBarrier(
			srcStage = IntFlags.of(PipelineStage.TopOfPipe),
			dstStage = IntFlags.of(PipelineStage.FragmentShader),
			images = listOf(
				field.blurredOcclusionImage.image.barrier(
					dstAccess = IntFlags.of(Access.ShaderRead),
					newLayout = Image.Layout.ShaderReadOnlyOptimal
				)
			)
		)
	}

	fun render(cmdbuf: CommandBuffer, field: OcclusionCalculator.Field) = cmdbuf.run {

		// render it!
		bindPipeline(renderPipeline)
		bindDescriptorSet(slideRenderer.mainDescriptorSet, renderPipeline)
		draw(field.extent.run { width*height*depth })
	}
}
