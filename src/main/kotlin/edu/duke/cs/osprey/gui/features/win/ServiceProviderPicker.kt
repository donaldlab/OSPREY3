package edu.duke.cs.osprey.gui.features.win

import cuchaz.kludge.imgui.Commands
import cuchaz.kludge.tools.IntFlags
import cuchaz.kludge.tools.Ref
import edu.duke.cs.osprey.gui.io.UserSettings
import edu.duke.cs.osprey.molscope.gui.WindowCommands
import edu.duke.cs.osprey.molscope.gui.WindowFeature
import edu.duke.cs.osprey.molscope.gui.enabledIf
import edu.duke.cs.osprey.molscope.gui.features.FeatureId
import edu.duke.cs.osprey.molscope.gui.features.WindowState
import kotlinx.coroutines.runBlocking
import java.util.concurrent.atomic.AtomicReference
import kotlin.concurrent.thread
import edu.duke.cs.osprey.service.OspreyService as OspreyServiceServer
import edu.duke.cs.osprey.gui.io.OspreyService as OspreyServiceClient
import kotlin.math.min


class ServiceProviderPicker : WindowFeature {

	override val id = FeatureId("settings.serviceProviders")

	private val winState = WindowState()
	private val addTitle = "Add a new Service Provider"
	private val editTitle = "Edit a Service Provider"

	private var selected: UserSettings.ServiceProvider? = null

	private val textHostname = Commands.TextBuffer(1024)
	private val pPort = Ref.of(0)
	private val pHttps = Ref.of(false)

	private var testResults = AtomicReference("")

	override fun menu(imgui: Commands, win: WindowCommands) = imgui.run {
		if (menuItem("Service Providers")) {
			winState.isOpen = true
		}
	}

	override fun gui(imgui: Commands, win: WindowCommands) = imgui.run {

		// render the main window
		winState.render(
			whenOpen = {

				// draw the window
				window("Service Providers", winState.pOpen, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

					// show the active service provider
					text("Active service provider:")
					indent(20f)

					text(UserSettings.serviceProvider.toString())

					if (button("test")) {
						test()
					}
					sameLine(spacing = 20f)
					text(testResults.get())

					unindent(20f)

					spacing()

					// show the current service providers
					text("Available service providers:")
					child("providers", 500f, 300f, border = true) {
						for (provider in UserSettings.serviceProviders) {

							// convert the provider into a unique label
							val label = provider.toString()

							if (radioButton("##$provider", UserSettings.serviceProvider === provider)) {
								set(provider)
							}
							sameLine()
							if (selectable(label, selected === provider)) {
								selected = provider
							}
						}
					}

					spacing()

					if (button("Add")) {
						textHostname.text = ""
						pPort.value = OspreyServiceServer.defaultPort
						pHttps.value = true
						openPopup(addTitle)
					}

					val selected = selected

					sameLine(spacing = 20f)

					enabledIf(selected != null) {
						if (button("Edit") && selected != null) {
							textHostname.text = selected.hostname
							pPort.value = selected.port
							pHttps.value = selected.https
							openPopup(editTitle)
						}
					}

					sameLine(spacing = 20f)

					// allow removal if more than one left
					enabledIf(UserSettings.serviceProviders.size > 1 && selected != null) {
						if (button("Remove") && selected != null) {
							remove(selected)
						}
					}

					// render any popups, if they're open
					renderAdd(imgui)
					renderEdit(imgui)
				}
			}
		)
	}

	private fun renderAdd(imgui: Commands) = imgui.run {

		popupModal(addTitle, null, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

			spacing()
			providerForm(imgui)
			spacing()

			if (button("Add")) {
				add(UserSettings.ServiceProvider(
					textHostname.text,
					pPort.value,
					pHttps.value
				))
				closeCurrentPopup()
			}

			sameLine(pos = 300f)

			if (button("Cancel")) {
				closeCurrentPopup()
			}
		}
	}

	private fun renderEdit(imgui: Commands) = imgui.run {

		val selected = selected ?: return

		popupModal(editTitle, null, IntFlags.of(Commands.BeginFlags.AlwaysAutoResize)) {

			spacing()
			providerForm(imgui)
			spacing()

			if (button("Save")) {
				edit(selected, UserSettings.ServiceProvider(
					textHostname.text,
					pPort.value,
					pHttps.value
				))
				closeCurrentPopup()
			}

			sameLine(pos = 300f)

			if (button("Cancel")) {
				closeCurrentPopup()
			}
		}
	}

	private fun providerForm(imgui: Commands) = imgui.run {

		inputText("Hostname", textHostname)
		inputInt("Port", pPort)
		checkbox("Https", pHttps)
	}

	private fun set(provider: UserSettings.ServiceProvider) {

		UserSettings.serviceProvider = provider
		// NOTE: setting serviceProvider automatically saves the new settings

		// reset the test results
		testResults.set("")
	}

	private fun add(provider: UserSettings.ServiceProvider) {

		UserSettings.serviceProviders.add(provider)

		// changing the serviceProviders list doesn't automatically save, so explictly save now
		UserSettings.save()
	}

	private fun edit(old: UserSettings.ServiceProvider, new: UserSettings.ServiceProvider) {

		// replace the old provider with the new one
		UserSettings.serviceProviders.replaceAll {
			if (it === old) {
				new
			} else {
				old
			}
		}

		// apply the changes
		if (UserSettings.serviceProvider === old) {
			UserSettings.serviceProvider = new
			// NOTE: setting serviceProvider automatically saves the new settings
		} else {
			// changing the serviceProviders list doesn't automatically save, so explictly save now
			UserSettings.save()
		}

		// update the selection too
		selected = new

		// reset the test results
		testResults.set("")
	}

	private fun remove(provider: UserSettings.ServiceProvider) {

		// can't remove the last provider
		if (UserSettings.serviceProviders.size <= 1) {
			return
		}

		// remove the given provider
		val i = UserSettings.serviceProviders
			.indexOfFirst { it === provider }
			.takeIf { it >= 0 }
			?: return
		UserSettings.serviceProviders.removeAt(i)

		// pick the next provider arbitrarily
		UserSettings.serviceProvider = UserSettings.serviceProviders[min(i, UserSettings.serviceProviders.size - 1)]
		// NOTE: setting serviceProvider automatically saves the new settings

		// reset the test results
		testResults.set("")
	}

	private fun test() {
		thread(
			name = "Service Tester",
			isDaemon = true,
			start = true
		) {
			try {
				testResults.set("Testing ...")
				runBlocking {
					val response = OspreyServiceClient.about()
					testResults.set("Success: ${response.name} ${response.version}")
				}
			} catch (t: Throwable) {
				testResults.set("Failed!")
			}
		}
	}
}
