package edu.duke.cs.osprey.molscope.gui.features.slide

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.molscope.Slide
import edu.duke.cs.osprey.molscope.gui.SlideFeature
import edu.duke.cs.osprey.molscope.gui.SlideCommands
import edu.duke.cs.osprey.molscope.gui.disabledIf
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.render.RenderSettings


class MenuRenderSettings(private val initialSettings: RenderSettings) : SlideFeature {

	override val id = FeatureId("rendersettings")

	val pOpen = Ref.of(false)

	var isInitialized = false
		private set

	override fun menu(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {
		if (menuItem("Render Settings")) {
			pOpen.value = true
		}
	}

	override fun gui(imgui: Commands, slide: Slide.Locked, slidewin: SlideCommands) = imgui.run {

		// apply initial settings if needed
		if (!isInitialized) {
			slidewin.renderSettings.set(initialSettings)
			isInitialized = true
		}

		if (pOpen.value) {

			window("Render Settings##${slide.name}", pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

				// make sliders for the render settings
				sliderFloat("Color", Ref.of(slidewin.renderSettings::colorWeight), 0f, 1f, "%.2f")
				sliderFloat("Shading", Ref.of(slidewin.renderSettings::shadingWeight), 0f, 2f, "%.2f")
				sliderFloat("Light Intensity", Ref.of(slidewin.renderSettings::lightWeight), 0f, 4f, "%.2f")
				sliderFloat("Depth Fade", Ref.of(slidewin.renderSettings::depthWeight), 0f, 1f, "%.2f")
				indent(20f) {
					checkbox("Auto", Ref.of(slidewin.renderSettings::autoDepth))
					disabledIf(slidewin.renderSettings.autoDepth) {
						itemWidth(140f) {
							sliderFloat("Min Depth", Ref.of(slidewin.renderSettings::depthZMin), 0f, 1f, "%0.2f")
							sliderFloat("Max Depth", Ref.of(slidewin.renderSettings::depthZMax), 0f, 1f, "%0.2f")
						}
					}
				}
				// TEMP: disable ambient occlusion, since it can hang the window for a bit
				//sliderFloat("Ambient Occlusion", Ref.of(slidewin.renderSettings::ambientOcclusionWeight), 0f, 4f, "%.2f")
				// TODO: optimize ambient occlusion so it's not using a brute force calculation and doesn't chug so hard
			}
		}
	}
}
