package edu.duke.cs.osprey.molscope.render

import cuchaz.kludge.tools.*
import cuchaz.kludge.vulkan.*
import cuchaz.kludge.vulkan.Queue
import edu.duke.cs.osprey.molscope.CameraCommand
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.shaders.Shaders
import edu.duke.cs.osprey.molscope.view.MoleculeRenderView
import edu.duke.cs.osprey.molscope.view.RenderView
import org.joml.Vector3f
import java.nio.ByteBuffer
import kotlin.reflect.KProperty


internal class SlideRenderer(
	val queue: Queue,
	val width: Int,
	val height: Int,
	oldRenderer: SlideRenderer? = null,
	backgroundColor: ColorRGBA = ColorRGBA.Float(0f, 0f, 0f)
) : AutoCloseable {

	private val closer = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose() = apply { closer.add(this) }
	override fun close() = closer.close()

	val device get() = queue.device

	val extent = Extent2D(width, height)
	val rect = Rect2D(Offset2D(0, 0), extent)

	var cursorPos: Offset2D? = null
	val cursorIndices = Indices()
	var cursorEffect: RenderEffect? = null

	// get the old settings (and mark them dirty) or make new settings
	val settings: RenderSettings = oldRenderer
		?.settings
		?.apply { dirty = true }
		?: RenderSettings()

	var backgroundColor: ColorRGBA = backgroundColor
		set(value) {
			field = value

			// sync the settings with the background color
			settings.backgroundColor = Vector3f(backgroundColor.rf, backgroundColor.gf, backgroundColor.bf)
		}

	// make the main render pass
	val colorAttachment =
		Attachment(
			format = Image.Format.R8G8B8A8_UNORM,
			loadOp = LoadOp.Clear,
			storeOp = StoreOp.Store,
			finalLayout = Image.Layout.ShaderReadOnlyOptimal
		)
	val indexAttachment =
		Attachment(
			format = Image.Format.R32G32_SINT,
			loadOp = LoadOp.Clear,
			storeOp =  StoreOp.Store,
			finalLayout = Image.Layout.General
		)
	val effectsAttachment =
		Attachment(
			format = Image.Format.R8G8B8A8_UINT,
			loadOp = LoadOp.Clear,
			storeOp = StoreOp.Store,
			finalLayout = Image.Layout.General
		)
	val depthAttachment =
		Attachment(
			format = Image.Format.D32_SFLOAT,
			loadOp = LoadOp.Clear,
			storeOp = StoreOp.Store,
			finalLayout = Image.Layout.DepthStencilAttachmentOptimal
		)
	val subpass =
		Subpass(
			pipelineBindPoint = PipelineBindPoint.Graphics,
			colorAttachments = listOf(
				colorAttachment to Image.Layout.ColorAttachmentOptimal,
				indexAttachment to Image.Layout.ColorAttachmentOptimal,
				effectsAttachment to Image.Layout.ColorAttachmentOptimal
			),
			depthStencilAttachment = depthAttachment to Image.Layout.DepthStencilAttachmentOptimal
		)
	val renderPass = device
		.renderPass(
			attachments = listOf(colorAttachment, indexAttachment, effectsAttachment, depthAttachment),
			subpasses = listOf(subpass),
			subpassDependencies = listOf(
				SubpassDependency(
					src = Subpass.External.dependency(
						stage = IntFlags.of(PipelineStage.ColorAttachmentOutput)
					),
					dst = subpass.dependency(
						stage = IntFlags.of(PipelineStage.ColorAttachmentOutput),
						access = IntFlags.of(Access.ColorAttachmentRead, Access.ColorAttachmentWrite)
					)
				)
			)
		)
		.autoClose()

	// make the post processing render pass
	val postSubpass =
		Subpass(
			pipelineBindPoint = PipelineBindPoint.Graphics,
			colorAttachments = listOf(colorAttachment to Image.Layout.ColorAttachmentOptimal)
		)
	val postRenderPass = device
		.renderPass(
			attachments = listOf(colorAttachment),
			subpasses = listOf(postSubpass),
			subpassDependencies = listOf(
				SubpassDependency(
					src = Subpass.External.dependency(
						stage = IntFlags.of(PipelineStage.ColorAttachmentOutput)
					),
					dst = postSubpass.dependency(
						stage = IntFlags.of(PipelineStage.ColorAttachmentOutput),
						access = IntFlags.of(Access.ColorAttachmentRead, Access.ColorAttachmentWrite)
					)
				)
			)
		)
		.autoClose()

	// make the render images
	val colorImage = device
		.image(
			Image.Type.TwoD,
			extent.to3D(1),
			colorAttachment.format,
			IntFlags.of(Image.Usage.ColorAttachment, Image.Usage.Sampled)
		)
		.autoClose()
		.allocateDevice()
		.autoClose()
	val colorView = colorImage.image.view().autoClose()

	val indexImage = device
		.image(
			Image.Type.TwoD,
			extent.to3D(1),
			indexAttachment.format,
			IntFlags.of(Image.Usage.ColorAttachment, Image.Usage.Sampled)
		)
		.autoClose()
		.allocateDevice()
		.autoClose()
	val indexView = indexImage.image.view().autoClose()

	val effectsImage = device
		.image(
			Image.Type.TwoD,
			extent.to3D(1),
			effectsAttachment.format,
			IntFlags.of(Image.Usage.ColorAttachment, Image.Usage.Sampled)
		)
		.autoClose()
		.allocateDevice()
		.autoClose()
	val effectsView = effectsImage.image.view().autoClose()

	val postImage = device
		.image(
			Image.Type.TwoD,
			extent.to3D(1),
			colorAttachment.format,
			IntFlags.of(Image.Usage.ColorAttachment, Image.Usage.Sampled)
		)
		.autoClose()
		.allocateDevice()
		.autoClose()
	val postView = postImage.image.view().autoClose()

	val sampler = device.sampler().autoClose()

	// make the depth buffer
	val depth = device
		.image(
			Image.Type.TwoD,
			extent.to3D(1),
			depthAttachment.format,
			IntFlags.of(Image.Usage.DepthStencilAttachment)
		)
		.autoClose()
		.allocateDevice()
		.autoClose()
	val depthView = depth.image.view(
		range = Image.SubresourceRange(aspectMask = IntFlags.of(Image.Aspect.Depth))
	).autoClose()

	// make a framebuffer
	val framebuffer = device
		.framebuffer(
			renderPass,
			imageViews = listOf(colorView, indexView, effectsView, depthView),
			extent = extent
		)
		.autoClose()

	// make a framebuffer for the post pass too
	val postFramebuffer = device
		.framebuffer(
			postRenderPass,
			imageViews = listOf(postView),
			extent = extent
		)
		.autoClose()

	// make a camera
	val camera: Camera = Camera(device)
		.autoClose()
		.apply {

			// if we have an old renderer, copy the camera state
			if (oldRenderer != null) {

				set(oldRenderer.camera)
				resize(width, height)
			}
		}
	var cameraSequence = -1

	// make the descriptor pool
	val descriptorPool = device.descriptorPool(
		maxSets = 3,
		sizes = DescriptorType.Counts(
			DescriptorType.UniformBuffer to 3,
			DescriptorType.StorageBuffer to 2,
			DescriptorType.CombinedImageSampler to 5
		)
	).autoClose()

	// make the main descriptor set
	val viewBufBinding = DescriptorSetLayout.Binding(
		binding = 0,
		type = DescriptorType.UniformBuffer,
		stages = IntFlags.of(ShaderStage.Vertex, ShaderStage.Geometry, ShaderStage.Fragment)
	)
	val occlusionImageBinding = DescriptorSetLayout.Binding(
		binding = 1,
		type = DescriptorType.CombinedImageSampler,
		stages = IntFlags.of(ShaderStage.Fragment)
	)
	val boundsBinding = DescriptorSetLayout.Binding(
		binding = 2,
		type = DescriptorType.UniformBuffer,
		stages = IntFlags.of(ShaderStage.Vertex, ShaderStage.Fragment)
	)
	val settingsBinding = DescriptorSetLayout.Binding(
		binding = 3,
		type = DescriptorType.UniformBuffer,
		stages = IntFlags.of(ShaderStage.Fragment)
	)
	val mainDescriptorSetLayout = device.descriptorSetLayout(listOf(
		viewBufBinding, occlusionImageBinding, boundsBinding, settingsBinding
	)).autoClose()
	val mainDescriptorSet = descriptorPool.allocate(mainDescriptorSetLayout)

	// make the cursor descriptor set
	val cursorBufBinding = DescriptorSetLayout.Binding(
		binding = 0,
		type = DescriptorType.StorageBuffer,
		// TODO: effects in geometry shader too? (eg, expand billboards for fades/blurs?)
		stages = IntFlags.of(ShaderStage.Compute, ShaderStage.Fragment)
	)
	val indexImageBinding = DescriptorSetLayout.Binding(
		binding = 2,
		type = DescriptorType.CombinedImageSampler,
		stages = IntFlags.of(ShaderStage.Compute, ShaderStage.Fragment)
	)
	val cursorDescriptorSetLayout = device.descriptorSetLayout(listOf(
		cursorBufBinding, indexImageBinding
	)).autoClose()
	val cursorDescriptorSet = descriptorPool.allocate(cursorDescriptorSetLayout)

	// make the post descriptor set
	val colorImageBinding = DescriptorSetLayout.Binding(
		binding = 1,
		type = DescriptorType.CombinedImageSampler,
		stages = IntFlags.of(ShaderStage.Fragment)
	)
	val effectsImageBinding = DescriptorSetLayout.Binding(
		binding = 3,
		type = DescriptorType.CombinedImageSampler,
		stages = IntFlags.of(ShaderStage.Compute, ShaderStage.Fragment)
	)
	val postDescriptorSetLayout = device.descriptorSetLayout(listOf(
		cursorBufBinding, colorImageBinding, indexImageBinding, effectsImageBinding
	)).autoClose()
	val postDescriptorSet = descriptorPool.allocate(postDescriptorSetLayout)

	// make a graphics command buffer
	val commandPool = device
		.commandPool(
			queue.family,
			flags = IntFlags.of(CommandPool.Create.ResetCommandBuffer)
		)
		.autoClose()
	val commandBuffer = commandPool.buffer()

	// allocate the settings buffer
	val settingsBuf = device
		.buffer(
			size = RenderSettings.bufferSize,
			usage = IntFlags.of(Buffer.Usage.UniformBuffer, Buffer.Usage.TransferDst)
		)
		.autoClose()
		.allocateDevice()
		.autoClose()

	// allocate the cursor buffer on the device
	val cursorBufDevice = device
		.buffer(
			size = Int.SIZE_BYTES*12L,
			usage = IntFlags.of(Buffer.Usage.StorageBuffer, Buffer.Usage.TransferDst, Buffer.Usage.TransferSrc)
		)
		.autoClose()
		.allocateDevice()
		.autoClose()

	// allocate another cursor buffer on the host
	val cursorBufHost = device
		.buffer(
			size = cursorBufDevice.buffer.size,
			usage = IntFlags.of(Buffer.Usage.TransferDst, Buffer.Usage.TransferSrc)
		)
		.autoClose()
		.allocateHost()
		.autoClose()

	init {
		// update the descriptor sets
		device.updateDescriptorSets(
			writes = listOf(
				mainDescriptorSet.address(viewBufBinding).write(
					DescriptorSet.BufferInfo(camera.buf.buffer)
				),
				// AmbientOcclusion updates some of the main bindings
				mainDescriptorSet.address(settingsBinding).write(
					DescriptorSet.BufferInfo(settingsBuf.buffer)
				),
				cursorDescriptorSet.address(cursorBufBinding).write(
					DescriptorSet.BufferInfo(cursorBufDevice.buffer)
				),
				cursorDescriptorSet.address(indexImageBinding).write(
					DescriptorSet.ImageInfo(sampler, indexView, Image.Layout.General)
				),
				postDescriptorSet.address(cursorBufBinding).write(
					DescriptorSet.BufferInfo(cursorBufDevice.buffer)
				),
				postDescriptorSet.address(colorImageBinding).write(
					DescriptorSet.ImageInfo(sampler, colorView, Image.Layout.ShaderReadOnlyOptimal)
				),
				postDescriptorSet.address(indexImageBinding).write(
					DescriptorSet.ImageInfo(sampler, indexView, Image.Layout.General)
				),
				postDescriptorSet.address(effectsImageBinding).write(
					DescriptorSet.ImageInfo(sampler, effectsView, Image.Layout.General)
				)
			)
		)
	}

	/**
	 * next descriptor set layout binding = 3
	 * next push constant range offset = 16
	 */
	fun graphicsPipeline(
		stages: List<ShaderModule.Stage>,
		vertexInput: VertexInput = VertexInput(),
		inputAssembly: InputAssembly,
		descriptorSetLayouts: List<DescriptorSetLayout> = emptyList(),
		pushConstantRanges: List<PushConstantRange> = emptyList()
	) = device.graphicsPipeline(
		renderPass,
		stages,
		descriptorSetLayouts = listOf(mainDescriptorSetLayout) + descriptorSetLayouts,
		pushConstantRanges = listOf(
				PushConstantRange(IntFlags.of(ShaderStage.Fragment), 16)
			) + pushConstantRanges,
		vertexInput = vertexInput,
		inputAssembly = inputAssembly,
		rasterizationState = RasterizationState(
			cullMode = IntFlags.of(CullMode.Back),
			frontFace = FrontFace.Counterclockwise
		),
		viewports = listOf(Viewport(
			0.0f,
			0.0f,
			extent.width.toFloat(),
			extent.height.toFloat()
		)),
		scissors = listOf(rect),
		colorAttachmentBlends = listOf(

			// use typical alpha blending
			colorAttachment to ColorBlendState.Attachment(
				color = ColorBlendState.Attachment.Part(
					src = BlendFactor.SrcAlpha,
					dst = BlendFactor.OneMinusSrcAlpha,
					op = BlendOp.Add
				),
				alpha = ColorBlendState.Attachment.Part(
					src = BlendFactor.One,
					dst = BlendFactor.One,
					op = BlendOp.Max
				)
			),

			// disable blending to always overwrite the dest (framebuf) values
			indexAttachment to ColorBlendState.Attachment(),
			effectsAttachment to ColorBlendState.Attachment()
		),
		depthStencilState = DepthStencilState()
	)

	// make a compute shader to download the index under the cursor
	private val cursorPipeline = device
		.computePipeline(
			stage = device.shaderModule(Shaders["cursorIndex.comp"])
				.autoClose()
				.stage("main", ShaderStage.Compute),
			descriptorSetLayouts = listOf(cursorDescriptorSetLayout)
		).autoClose()

	private val postPipeline = device
		.graphicsPipeline(
			postRenderPass,
			stages = listOf(
				device.shaderModule(Shaders["post.vert"])
					.autoClose()
					.stage("main", ShaderStage.Vertex),
				device.shaderModule(Shaders["post.frag"])
					.autoClose()
					.stage("main", ShaderStage.Fragment)
			),
			descriptorSetLayouts = listOf(postDescriptorSetLayout),
			inputAssembly = InputAssembly(InputAssembly.Topology.TriangleStrip),
			rasterizationState = RasterizationState(
				cullMode = IntFlags.of(CullMode.Back),
				frontFace = FrontFace.Counterclockwise
			),
			viewports = listOf(Viewport(
				0.0f,
				0.0f,
				extent.width.toFloat(),
				extent.height.toFloat()
			)),
			scissors = listOf(rect),
			colorAttachmentBlends = listOf(
				colorAttachment to ColorBlendState.Attachment(
					color = ColorBlendState.Attachment.Part(
						src = BlendFactor.SrcAlpha,
						dst = BlendFactor.OneMinusSrcAlpha,
						op = BlendOp.Add
					),
					alpha = ColorBlendState.Attachment.Part(
						src = BlendFactor.One,
						dst = BlendFactor.One,
						op = BlendOp.Max
					)
				)
			)
		)
		.autoClose()

	private val sphereRenderer = SphereRenderer(this).autoClose()
	private val cylinderRenderer = CylinderRenderer(this).autoClose()

	private val occlusionRenderer = OcclusionRenderer(this).autoClose()

	fun render(slide: Slide.Locked, renderables: ViewRenderables, occlusionField: OcclusionCalculator.Field? = null, renderFinished: Semaphore? = null) {

		sphereRenderer.update(renderables.spheres)
		cylinderRenderer.update(renderables.cylinders)

		// update the camera
		while (slide.camera.queue.isNotEmpty()) {
			val cmd = slide.camera.queue.pollFirst() ?: break
			when (cmd) {
				is CameraCommand.LookAtBox -> {
					camera.lookAtBox(
						width, height,
						focalLength = 200f,
						look = Vector3f(0f, 0f, 1f),
						up = Vector3f(0f, 1f, 0f),
						box = cmd.aabb
					)
					camera.changed()
				}
			}
		}

		// did the camera change?
		if (camera.sequence > cameraSequence) {

			// yup, upload the new camera to the GPU
			cameraSequence = camera.sequence
			camera.upload()

			// update the depth settings, if needed
			if (settings.autoDepth) {
				settings.updateDepth(camera, slide.views)
			}
		}

		// update the hover buffer
		cursorBufHost.memory.map { buf ->

			// update the pos
			val cursorPos = cursorPos
			if (cursorPos != null) {
				buf.putInt(1)
				buf.skip(Int.SIZE_BYTES*3)
				buf.putInt(cursorPos.x)
				buf.putInt(cursorPos.y)
			} else {
				buf.putInt(0)
				buf.skip(Int.SIZE_BYTES*5)
			}

			// skip the indices (the GPU writes those, the CPU just reads them)
			buf.skip(Int.SIZE_BYTES*2)

			// update the cursor effect
			val cursorEffect = cursorEffect
			if (cursorEffect != null) {
				// need to expand bytes to ints for the cursor buffer
				buf.putInt(cursorEffect.r.toInt())
				buf.putInt(cursorEffect.g.toInt())
				buf.putInt(cursorEffect.b.toInt())
				buf.putInt(cursorEffect.flags.value.toInt())
			} else {
				buf.putInts(0, 0, 0, 0)
			}

			buf.flip()
		}

		// did the render settings change?
		if (settings.dirty) {

			// yup, upload them to the GPU
			settingsBuf.transferHtoD { buf ->
				buf.putSettings(settings)
				buf.flip()
			}
			settings.dirty = false
		}

		// record the command buffer
		commandBuffer.apply {
			begin(IntFlags.of(CommandBuffer.Usage.OneTimeSubmit))

			// upload the cursor buffer
			copyBuffer(cursorBufHost.buffer, cursorBufDevice.buffer)

			if (occlusionField != null) {
				occlusionRenderer.barriers(this, occlusionField)
			}

			// draw all the views
			beginRenderPass(
				renderPass,
				framebuffer,
				rect,
				clearValues = mapOf(
					colorAttachment to backgroundColor.toClearColor(),
					indexAttachment to Indices.clearColor,
					effectsAttachment to RenderEffect.clearColor,
					depthAttachment to ClearValue.DepthStencil(depth = 1f)
				)
			)
			slide.views.forEachIndexed { i, view ->
				if (view.isVisible) {
					view.spheres?.let { sphereRenderer.render(this, it, i) }
				}
			}
			slide.views.forEachIndexed { i, view ->
				if (view.isVisible) {
					view.cylinders?.let { cylinderRenderer.render(this, it, i) }
				}
			}
			if (occlusionField != null && settings.showOcclusionField) {
				occlusionRenderer.render(this, occlusionField)
			}
			endRenderPass()

			// figure out what was under the cursor, if needed
			if (cursorPos != null) {

				// wait for the render pass to finish
				pipelineBarrier(
					srcStage = IntFlags.of(PipelineStage.ColorAttachmentOutput),
					dstStage = IntFlags.of(PipelineStage.ComputeShader)
				)

				bindPipeline(cursorPipeline)
				bindDescriptorSet(cursorDescriptorSet, cursorPipeline)
				dispatch(1)

				// download the cursor buffer so we can read the index on the host side
				copyBuffer(cursorBufDevice.buffer, cursorBufHost.buffer)

				// wait for the compute shader to finish before running the post shader
				pipelineBarrier(
					srcStage = IntFlags.of(PipelineStage.ComputeShader),
					dstStage = IntFlags.of(PipelineStage.FragmentShader)
				)
			}

			// do the post-processing pass
			beginRenderPass(
				postRenderPass,
				postFramebuffer,
				rect,
				clearValues = mapOf(colorAttachment to ClearValue.Color.Int(0, 0, 0, 0))
			)
			bindPipeline(postPipeline)
			bindDescriptorSet(postDescriptorSet, postPipeline)
			draw(4)
			endRenderPass()

			end()
		}

		// render the frame
		queue.submit(
			commandBuffer,
			waitFor = listOf(),
			signalTo =
				if (renderFinished != null) {
					listOf(renderFinished)
				} else {
					emptyList()
				}
		)

		// read the cursor indices, if any
		if (cursorPos != null) {

			// TODO: is there a more efficient wait mechanism here?
			// we only need to wait for the buffer download, not the whole post-processing step
			queue.waitForIdle()

			cursorBufHost.memory.map { buf ->
				buf.position = 24
				buf.getIndices(cursorIndices)
			}

		} else {
			cursorIndices.clear()
		}
	}
}

data class Indices(var target: Int, var view: Int) {

	companion object {
		const val noTarget = -1
		const val noView = -1
		val clearColor = ClearValue.Color.Int(noTarget, noView, 0, 0)
	}

	constructor() : this(noTarget, noView)

	fun clear() {
		target = noTarget
		view = noView
	}

	val isEmpty get() = target == noTarget && view == noView
}

fun ByteBuffer.put(indices: Indices) {
	putInt(indices.target)
	putInt(indices.view)
}

fun ByteBuffer.getIndices(indices: Indices) {
	indices.target = int
	indices.view = int
}


internal data class ViewRenderables(
	val spheres: List<SphereRenderable>,
	val cylinders: List<CylinderRenderable>
)

class RenderSettings {

	/** set true if we need to upload the render settings buffer */
	internal var dirty = true

	private inner class Dirtyable<T>(var value: T) {

		operator fun getValue(thisRef: Any?, property: KProperty<*>) = value

		operator fun setValue(thisRef: Any?, property: KProperty<*>, value: T) {
			dirty = dirty || value != this.value
			this.value = value
		}
	}

	companion object {
		val bufferSize: Long = Float.SIZE_BYTES*10L
	}

	var backgroundColor: Vector3f by Dirtyable(Vector3f(0f, 0f, 0f))
	var colorWeight: Float by Dirtyable(1f)
	var lightWeight: Float by Dirtyable(1f)
	var shadingWeight: Float by Dirtyable(1f)
	var depthWeight: Float by Dirtyable(0.2f)
	var depthZMin: Float by Dirtyable(0f)
	var depthZMax: Float by Dirtyable(1f)
	var ambientOcclusionWeight: Float by Dirtyable(0f) // default to 0, so occlusion is off

	var showOcclusionField: Boolean = false
	var autoDepth: Boolean = true

	fun set(other: RenderSettings) {

		this.backgroundColor = other.backgroundColor
		this.colorWeight = other.colorWeight
		this.lightWeight = other.lightWeight
		this.shadingWeight = other.shadingWeight
		this.depthWeight = other.depthWeight
		this.depthZMin = other.depthZMin
		this.depthZMax = other.depthZMax
		this.ambientOcclusionWeight = other.ambientOcclusionWeight

		this.showOcclusionField = other.showOcclusionField
		this.autoDepth = other.autoDepth
	}

	fun updateDepth(camera: Camera, views: List<RenderView>, margin: Float = 0.05f) {

		// get atom distances to the camera
		val distances = views
			.filterIsInstance<MoleculeRenderView>()
			.flatMap { it.currentMol.atoms }
			.map { it.pos.toFloat().sub(camera.pos).parallelTo(camera.look).lengthSquared() }

		if (distances.isEmpty()) {
			return
		}

		// get the min,max distances
		var minDist = distances.min()!!.sqrt()
		var maxDist = distances.max()!!.sqrt()

		// add the margin
		val pad = (maxDist - minDist)*margin/2
		minDist -= pad
		maxDist += pad

		// normalize the distances into the [zNear,zFar] range
		val zLen = camera.zFar - camera.zNear
		depthZMin = ((minDist - camera.zNear)/zLen).atLeast(0f)
		depthZMax = ((maxDist - camera.zNear)/zLen).atMost(1f)
	}
}

private fun ByteBuffer.putSettings(settings: RenderSettings) {
	putFloats(settings.backgroundColor[0], settings.backgroundColor[1], settings.backgroundColor[2])
	putFloat(settings.colorWeight)
	putFloat(settings.shadingWeight)
	putFloat(settings.lightWeight)
	putFloat(settings.depthWeight)
	putFloat(settings.depthZMin)
	putFloat(settings.depthZMax)
	putFloat(settings.ambientOcclusionWeight)
}
