package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.service.read
import edu.duke.cs.osprey.service.OspreyService
import org.tomlj.Toml
import java.nio.file.Files
import java.nio.file.Paths


object UserSettings {

	private val file = Paths.get(System.getProperty("user.home"))
		.resolve(".osprey")
		.resolve("settings.toml")

	data class ServiceProvider(
		val hostname: String,
		val port: Int = OspreyService.defaultPort,
		val https: Boolean = true
	) {

		override fun toString(): String =
			"${if (https) "https" else "https"}://$hostname:$port"
	}

	private var isLoading = true

	// init settings with defaults

	val serviceProviders = mutableListOf(
		ServiceProvider("olympias.cs.duke.edu", https=true)
	)
	var serviceProvider: ServiceProvider = serviceProviders.first()
		set(value) {
			if (value !in serviceProviders) {
				serviceProviders.add(value)
			}
			field = value
			if (!isLoading) {
				save()
			}
		}

	var openSaveDir = Paths.get(System.getProperty("user.home"))
		set(value) {
			field = value
			if (!isLoading) {
				save()
			}
		}

	init {
		// try to load values on first use
		if (file.exists) {
			load()
		} else {
			// no file? just save the defaults
			save()
		}
	}

	private fun load() {

		// don't let errors here crash the whole program
		try {

			isLoading = true

			val doc = Toml.parse(file.read())

			// load the service providers
			serviceProviders.clear()
			val serviceProvidersArray = doc.getArray("serviceProviders")
			if (serviceProvidersArray != null) {
				for (i in 0 until serviceProvidersArray.size()) {
					val providerTable = serviceProvidersArray.getTable(i)
					if (providerTable != null) {

						serviceProviders.add(ServiceProvider(
							hostname = providerTable.getString("host") ?: continue,
							port = providerTable.getInt("port") ?: continue,
							https = providerTable.getBoolean("https") ?: continue
						))
					}
				}
			}
			doc.getInt("serviceProvider")?.let { index ->
				serviceProviders.getOrNull(index)?.let { serviceProvider = it }
			}

			// load the open/save dir
			doc.getString("openSaveDir")?.let { pathname ->
				val path = Paths.get(pathname)
				if (Files.isDirectory(path)) {
					openSaveDir = path
				}
			}

		} catch (t: Throwable) {

			// just report the exception and move on
			// worst case, we just use default settings for this user
			t.printStackTrace(System.err)

		} finally {

			isLoading = false
		}
	}

	fun save() {

		val buf = StringBuilder()
		fun write(str: String, vararg args: Any) = buf.append(String.format(str, *args))

		// write the service providers
		if (serviceProviders.isNotEmpty()) {
			write("serviceProviders = [\n")
			for ((index, endpoint) in serviceProviders.withIndex()) {
				write("\t{ host = %s, port = %d, https = %b }, # %d\n",
					endpoint.hostname.quote(),
					endpoint.port,
					endpoint.https,
					index
				)
			}
			write("]\n")
		}
		serviceProviders
			.indexOfFirst { it === serviceProvider }
			.takeIf { it >= 0 }
			?.let {
				write("serviceProvider = %d\n", it)
			}

		// write the open/save dir
		write("openSaveDir = %s\n", openSaveDir.toAbsolutePath().toString().quote())

		Files.createDirectories(file.parent)
		buf.toString().write(file)
	}

	private fun String.quote() =
		replace("\\", "\\\\")
			.replace("\"", "\\\"")
			.let { str -> "\"$str\"" }
}
