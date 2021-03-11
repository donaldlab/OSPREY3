package edu.duke.cs.osprey.molscope

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.expand
import cuchaz.kludge.tools.toFloat
import cuchaz.kludge.vulkan.Extent2D
import edu.duke.cs.osprey.Osprey
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.Features
import edu.duke.cs.osprey.molscope.gui.features.slide.DevOcclusionField
import edu.duke.cs.osprey.molscope.view.RenderView
import org.joml.AABBd
import org.joml.AABBf
import java.util.*
import java.util.concurrent.atomic.AtomicInteger
import kotlin.collections.ArrayList


/**
 * a place to stage molecules and geoemtry for viewing
 *
 * slides must be thread-safe since they are directly accessed by the renderer
 */
class Slide(
	name: String,
	val initialSize: Extent2D? = null
) {

	var name: String = name
		private set

	inner class Locked {

		val unlocked: Slide = this@Slide

		var name: String
			get() = this@Slide.name
			set(value) { this@Slide.name = value }

		val views: MutableList<RenderView> = ArrayList()

		fun calcBoundingBox(): AABBf? {

			if (views.isEmpty()) {
				return null
			}

			return views
				.mapNotNull { view -> view.calcBoundingBox() }
				.reduce { a, b -> AABBf().apply { a.union(b, this) } }
		}

		inner class Camera {

			internal val queue = ArrayDeque<CameraCommand>()

			fun lookAtEverything(padding: Double = 1.0) {
				calcBoundingBox()?.let { lookAtBox(it, padding.toFloat()) }
			}

			fun lookAtBox(aabb: AABBf, padding: Float = 1f) {
				queue.addLast(CameraCommand.LookAtBox(aabb.apply { expand(padding) }))
			}

			fun lookAtBox(aabb: AABBd, padding: Double = 1.0) =
				lookAtBox(aabb.toFloat(), padding.toFloat())
		}

		val camera = Camera()

		inner class SlideFeatures {

			internal val features = Features<SlideFeature>()
			private var nextSpacingId = AtomicInteger(0)
			private var nextSeparatorId = AtomicInteger(0)

			fun <R> menu(name: String, block: SlideMenu.() -> R): R {
				val menu = features.menu(name)
				return object : SlideMenu {
					override fun add(feature: SlideFeature) = menu.add(feature)
					override fun addSpacing(num: Int) = menu.add(SlideSpacing(nextSpacingId.getAndIncrement(), num))
					override fun addSeparator() = menu.add(SlideSeparator(nextSeparatorId.getAndIncrement()))
				}.block()
			}
		}
		val features = SlideFeatures().apply {

			// add dev-only features if needed
			if (Osprey.dev) {
				menu("Dev") {
					add(DevOcclusionField())
				}
			}
		}
	}

	/**
	 * lock this slide to prevent races with the renderer
	 *
	 * don't talk to the window while the slide is locked, or you might deadlock!
	 */
	inline fun <R> lock(block: (Locked) -> R): R =
		synchronized(this) {
			block(_locked)
		}

	/**
	 * DO NOT USE THIS WITHOUT SYNCHRONIZING ON SLIDE OR YOU WILL RACE THE RENDER THREAD
	 * Use the lock() convenience method rather than accessing this directly
	 *
	 * If Kotlin/JVM were smarter, this would be private or internal
	 * but as things are now, we can't inline lock() when this is private or internal ;_;
	 */
	val _locked = Locked()
}


internal interface CameraCommand {
	data class LookAtBox(val aabb: AABBf) : CameraCommand
}


interface SlideMenu {
	fun add(feature: SlideFeature)
	fun addSpacing(num: Int = 1)
	fun addSeparator()
}

private class SlideSpacing(id: Int, val num: Int) : SlideFeature {

	override val id = FeatureId("SlideSpace_$id")

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		for (i in 0 until num) {
			spacing()
		}
	}
}

private class SlideSeparator(id: Int) : SlideFeature {

	override val id = FeatureId("SlideSeparator_$id")

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		separator()
	}
}
