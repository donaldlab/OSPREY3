package edu.duke.cs.osprey.molscope.gui

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.imgui.Imgui
import cuchaz.kludge.tools.*
import cuchaz.kludge.vulkan.*
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.render.*
import edu.duke.cs.osprey.molscope.tools.IdentityChangeTracker
import edu.duke.cs.osprey.molscope.view.*
import org.joml.Vector2f
import org.lwjgl.system.Platform
import kotlin.NoSuchElementException


internal class SlideWindow(
	val slide: Slide,
	val queue: Queue,
	val exceptionViewer: ExceptionViewer
): AutoCloseable {

	private val closer = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose(replace: R? = null) = also { closer.add(this@autoClose, replace) }
	override fun close() = closer.close()

	private inner class RendererInfo(
		val renderer: SlideRenderer,
		var imageDesc: Imgui.ImageDescriptor
	) {
		init {
			// if we already have an occlusion field, update the renderer right away
			occlusionField?.updateDescriptorSet(renderer)
		}
	}

	private var rendererInfo: RendererInfo? = null
	private val rendererInfoOrThrow: RendererInfo get() = rendererInfo ?: throw NoSuchElementException("no renderer yet, this is a bug")

	val renderFinished = queue.device.semaphore().autoClose()

	fun resizeIfNeeded() {

		// how big is the window content area?
		val width = (contentMax.x - contentMin.x).toInt()
		val height = (contentMax.y - contentMin.y).toInt()

		if (width <= 0 || height <= 0) {

			// not big enough, don't bother resizing
			return
		}

		val rendererInfo = this.rendererInfo
		if (rendererInfo == null) {

			// make a new renderer
			val renderer = SlideRenderer(queue, width, height).autoClose()
			val imageDesc = Imgui.imageDescriptor(renderer.postView, renderer.sampler).autoClose()
			this.rendererInfo = RendererInfo(renderer, imageDesc)

		} else {

			if (width == rendererInfo.renderer.width && height == rendererInfo.renderer.height) {

				// same size as before, don't resize
				return
			}

			// replace the old renderer
			val renderer = SlideRenderer(queue, width, height, rendererInfo.renderer).autoClose(replace = rendererInfo.renderer)
			val imageDesc = Imgui.imageDescriptor(renderer.postView, renderer.sampler).autoClose(replace = rendererInfo.imageDesc)
			this.rendererInfo = RendererInfo(renderer, imageDesc)
		}
	}

	fun updateImageDesc() {
		rendererInfo?.let {
			it.imageDesc = Imgui.imageDescriptor(it.renderer.postView, it.renderer.sampler)
				.autoClose(replace = it.imageDesc)
		}
	}

	private val renderablesTracker = RenderablesTracker()

	// OSX doesn't support some functions in the compute shader, so we can't do ambient occlusion on that platform
	private val occlusionCalculator: OcclusionCalculator? = when(Platform.get()) {
		Platform.MACOSX -> null
		else -> OcclusionCalculator(queue).autoClose()
	}

	private var occlusionField: OcclusionCalculator.Field? = null

	fun render(slide: Slide.Locked, renderFinished: Semaphore): Boolean {

		val rendererInfo = rendererInfo ?: return false

		rendererInfo.apply {

			// update the background color based on settings
			renderer.backgroundColor = backgroundColors.getValue(ColorsMode.current)

			// gather all the renderables by type
			val renderables = ViewRenderables(
				spheres = slide.views.mapNotNull { it.spheres },
				cylinders = slide.views.mapNotNull { it.cylinders }
			)

			// did any renderables change?
			renderablesTracker.update(renderables)
			if (renderablesTracker.changed) {

				// update the occlusion field
				occlusionField = occlusionCalculator?.run {
					Field(
						// TODO: make configurable?
						extent = Extent3D(16, 16, 16),
						gridSubdivisions = 2,
						renderables = renderables
					)
					.autoClose(replace = occlusionField)
					.apply {
						updateDescriptorSet(renderer)
					}
				}
			}

			// process the occlusion field if needed
			occlusionField?.run {
				active = rendererInfo.renderer.settings.ambientOcclusionWeight > 0f
				if (needsProcessing) {
					process()
				}
			}

			renderer.render(slide, renderables, occlusionField, renderFinished)
		}
		return true
	}


	// GUI state
	private val contentMin = Vector2f()
	private val contentMax = Vector2f()
	private var contextMenu: ContextMenu? = null

	private val commands = object : SlideCommands {

		override fun showExceptions(block: () -> Unit) {
			try {
				block()
			} catch (t: Throwable) {
				t.printStackTrace(System.err)
				exceptionViewer.add(t)
			}
		}

		override val renderSettings get() = rendererInfoOrThrow.renderer.settings
		override val extent get() = rendererInfoOrThrow.renderer.extent
		override val hoverEffects = HoverEffects()

		override var mouseTarget: ViewIndexed? = null
		override var mouseLeftClick = false
		override var mouseLeftRelease = false
		override var mouseLeftDrag = false
		override val mouseOffset = Vector2f()
		override val mouseLeftDragDelta = Vector2f()
		override var mouseWheelDelta = 0f

		override val camera get() = rendererInfoOrThrow.renderer.camera

		override fun loadImage(bytes: ByteArray) =
			LoadedImage(queue, bytes.toBuffer())
				.autoClose()
	}

	fun gui(imgui: Commands) = imgui.run {

		// to start, the window title is the slide name
		var title = slide.name

		// append rendering progress to the window title
		occlusionField?.let { field ->
			val progress = field.processingProgress
			if (progress < 1.0) {
				title += " (lighting ${"%.0f".format(progress*100.0)}%)"
			}
		}

		// add a unique id for this window
		title += "###${System.identityHashCode(slide)}"

		// start the window
		val initialSize = slide.initialSize
		if (initialSize != null) {
			setNextWindowSize(initialSize, Commands.Cond.Once)
		} else {
			setNextWindowSizeConstraints(
				320f, 240f,
				Float.POSITIVE_INFINITY, Float.POSITIVE_INFINITY
			)
		}
		if (!begin(title, flags = IntFlags.of(Commands.BeginFlags.MenuBar, Commands.BeginFlags.NoBringToFrontOnFocus))) {
			end()
			return
		}

		// render the slide feature menus
		slide.lock { slide ->
			if (beginMenuBar()) {

				if (rendererInfo != null) {

					// render feature menus
					for (menu in slide.features.features.menus) {
						if (beginMenu(menu.name)) {
							for (feature in menu.features) {
								feature.menu(this, slide, commands)
							}
							endMenu()
						}
					}
				}

				endMenuBar()
			}
		}

		// track the window content area
		getWindowContentRegionMin(contentMin)
		getWindowContentRegionMax(contentMax)

		val rendererInfo = rendererInfo ?: run {
			end()
			return
		}

		// draw the slide image
		setCursorPos(contentMin)
		image(rendererInfo.imageDesc)

		// draw a big invisible button over the image so we can capture mouse events
		setCursorPos(contentMin)
		invisibleButton("button", rendererInfo.renderer.extent)

		// what's the mouse looking at?
		commands.mouseTarget =
			if (rendererInfo.renderer.cursorIndices.isEmpty) {
				null
			} else {
				slide.lock { slide ->
					slide.views.getOrNull(rendererInfo.renderer.cursorIndices.view)?.let { view ->
						view.getIndexed(rendererInfo.renderer.cursorIndices.target)?.let { target ->
							ViewIndexed(view, target)
						}
					}
				}
			}

		// translate ImGUI mouse inputs into mouse state for features to consume
		commands.mouseLeftClick = isItemClicked(0)
		commands.mouseLeftRelease = isMouseReleased(0)
		commands.mouseLeftDrag = isItemActive() && Imgui.io.mouse.buttonDown[0]
		commands.mouseOffset
			.apply { getMousePos(this) }
			.sub(Vector2f().apply { getItemRectMin(this) })
		getMouseDragDelta(0, commands.mouseLeftDragDelta)
		commands.mouseWheelDelta = if (isItemHovered()) Imgui.io.mouse.wheel else 0f

		// handle context menus
		val isContextMenuOpen = isPopupOpen(ContextMenu.id)
		if (!isContextMenuOpen) {
			contextMenu = null

			// update the cursor in the renderer
			if (isItemHovered()) {
				rendererInfo.renderer.cursorPos = commands.mouseOffset.toOffset()
				rendererInfo.renderer.cursorEffect = commands.hoverEffects.get()
			} else {
				rendererInfo.renderer.cursorPos = null
				rendererInfo.renderer.cursorEffect = null
			}
		}
		popupContextItem(ContextMenu.id) {

			if (contextMenu == null) {
				val contextMenu = ContextMenu()

				// did we click on anything?
				commands.mouseTarget?.let { target ->

					// add slide features to the menu
					slide.lock { slide ->
						for (menu in slide.features.features.menus) {
							for (feature in menu.features) {
								feature.contextMenu(contextMenu, slide, commands, target)
							}
						}
					}
				}

				this@SlideWindow.contextMenu = contextMenu
			}

			// render the context menu if we have one
			contextMenu?.render(imgui)
		}

		end()

		// render the slide feature windows
		slide.lock { slide ->
			for (menu in slide.features.features.menus) {
				for (feature in menu.features) {
					feature.gui(this, slide, commands)
				}
			}
		}
	}
}


private val backgroundColors = mapOf(
	ColorsMode.Dark to ColorRGBA.Int(0, 0, 0),
	ColorsMode.Light to ColorRGBA.Int(255, 255, 255)
)


class ViewIndexed(val view: RenderView, val target: Any)


class ContextMenu {

	private val features = ArrayList<Commands.() -> Unit>()

	fun add(block: Commands.() -> Unit) {
		features.add(block)
	}

	fun render(imgui: Commands) = imgui.run {

		if (features.isEmpty()) {

			// no features, close the popup
			closeCurrentPopup()
		}

		for (feature in features) {
			feature()
		}
	}

	companion object {
		const val id = "contextMenu"
	}
}


internal class RenderablesTracker {

	private val spheres = IdentityChangeTracker<SphereRenderable>()
	private val cylinders = IdentityChangeTracker<CylinderRenderable>()

	var changed: Boolean = false
		private set

	fun update(renderables: ViewRenderables) {

		spheres.update(renderables.spheres)
		cylinders.update(renderables.cylinders)

		changed = spheres.changed || cylinders.changed
	}
}
