package edu.duke.cs.osprey.molscope.render

import cuchaz.kludge.tools.AutoCloser
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.toFlagsString
import cuchaz.kludge.vulkan.*
import edu.duke.cs.osprey.Osprey


class VulkanDevice(
	val vulkanExtensions: Set<String> = emptySet()
) : AutoCloseable {

	private val autoCloser = AutoCloser()
	private fun <R:AutoCloseable> R.autoClose() = also { autoCloser.add(this@autoClose) }
	override fun close() = autoCloser.close()

	private val canDebug = Osprey.dev && Vulkan.DebugExtension in Vulkan.supportedExtensions
	private val canReport = Osprey.dev && Vulkan.ReportExtension in Vulkan.supportedExtensions

	// make the main vulkan instance with the extensions we need
	val vulkan =
		Vulkan(
			extensionNames = vulkanExtensions.toMutableSet().apply {
				if (canDebug) {
					add(Vulkan.DebugExtension)
				} else if (canReport) {
					add(Vulkan.ReportExtension)
				}
			},
			layerNames = if (Osprey.dev) {
				mutableSetOf<String>().apply {
					if (Vulkan.StandardValidationLayer in Vulkan.supportedLayers) {
						add(Vulkan.StandardValidationLayer)
					}
					if (Vulkan.ValidationLayer in Vulkan.supportedLayers) {
						add(Vulkan.ValidationLayer)
					}
				}
			} else {
				emptySet()
			}
		)
		.autoClose()
		.apply {
			// listen to problems from vulkan, prefer the debug extension over the report extension
			if (canDebug) {
				debugMessenger(
					severities = IntFlags.of(DebugMessenger.Severity.Error, DebugMessenger.Severity.Warning)
				) { severity, type, msg ->
					println("VULKAN: ${severity.toFlagsString()} ${type.toFlagsString()} $msg")
					Exception("Stack Trace").printStackTrace(System.out)
				}.autoClose()
			} else if (canReport) {
				reportMessenger(
					flags = IntFlags.of(ReportMessenger.Flags.Error, ReportMessenger.Flags.Warning)
				) { flags, msg ->
					println("VULKAN: ${flags.toFlagsString()} $msg")
					Exception("Stack Trace").printStackTrace(System.out)
				}
				.autoClose()
			}
		}

	// pick a physical device: prefer discrete GPU
	// TODO: do we want to pick the physical device that is rendering the desktop?
	val physicalDevice = vulkan.physicalDevices
		.asSequence()
		.sortedBy { if (it.properties.type == PhysicalDevice.Type.DiscreteGpu) 0 else 1 }
		.first()

	// flag the features we need for when we create the logical device
	val deviceFeatures = PhysicalDevice.Features().apply {

		// required features

		// image formats
		if (!physicalDevice.features.shaderStorageImageExtendedFormats) {
			throw UnsupportedOperationException("Your GPU and/or video drivers don't support required image formats")
		}
		shaderStorageImageExtendedFormats = true

		// independent blend
		if (!physicalDevice.features.independentBlend) {
			throw UnsupportedOperationException("Your GPU and/or video drivers don't support independent blending, which is required for g-buffer rendering techniques")
		}
		independentBlend = true

		// optional features

		// geometry shaders
		if (physicalDevice.features.geometryShader) {
			geometryShader = true
		}
	}
}
