package edu.duke.cs.osprey.gui.features.win

import cuchaz.kludge.imgui.Commands
import edu.duke.cs.osprey.Osprey
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.WindowFeature
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.gui.forcefield.amber.OperatingSystem
import edu.duke.cs.osprey.gui.io.OspreyService
import edu.duke.cs.osprey.gui.io.UserSettings
import edu.duke.cs.osprey.service.OspreyService as Server
import java.nio.file.Paths


class LocalService : WindowFeature {

	override val id = FeatureId("service.local")

	override fun menu(imgui: Commands, win: WindowCommands) = imgui.run {

		// add a checkbox to toggle the local service
		checkbox("Local Service", LocalServiceRunner.isRunning)?.let { isChecked ->
			if (isChecked) {
				LocalServiceRunner.start()
			} else {
				LocalServiceRunner.close()
			}
		}

		Unit
	}
}


object LocalServiceRunner : AutoCloseable {

	private val serviceDir = Paths.get("")

	private var service: Server.Instance? = null

	val isRunning get() = service != null

	init {
		// if we're a linux developer, start a local service by default
		if (Osprey.dev && OperatingSystem.get() == OperatingSystem.Linux) {
			try {
				start()
			} catch (t: Throwable) {

				System.err.println("Failed to start local development service")
				t.printStackTrace(System.err)

				// weirdly enough, if these exeptions come from Netty, sometimes they have no stack traces!
				// super helpful, eh? So at least print something we can see in the logs
				if (t.stackTrace.isEmpty()) {
					Thread.dumpStack()
				}
			}
		}
	}

	fun start() {
		if (service == null) {
			service = Server.Instance(serviceDir, wait = false, useVersionPrefix=true)
			OspreyService.provider = UserSettings.ServiceProvider("localhost", https=false)
		}
	}

	override fun close() {
		service?.close()
		service = null
		OspreyService.provider = null
	}
}