package edu.duke.cs.osprey.service

import edu.duke.cs.osprey.Osprey
import edu.duke.cs.osprey.service.services.*
import io.ktor.application.call
import io.ktor.application.install
import io.ktor.client.features.json.serializer.*
import io.ktor.features.ContentNegotiation
import io.ktor.http.ContentType
import io.ktor.request.receive
import io.ktor.response.respond
import io.ktor.response.respondText
import io.ktor.routing.Routing
import io.ktor.routing.get
import io.ktor.routing.post
import io.ktor.routing.routing
import io.ktor.serialization.*
import io.ktor.server.engine.embeddedServer
import io.ktor.server.engine.stop
import io.ktor.server.netty.Netty
import kotlinx.serialization.json.Json
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.SerializersModule
import kotlinx.serialization.modules.polymorphic
import kotlinx.serialization.modules.subclass
import org.slf4j.LoggerFactory
import java.nio.charset.Charset
import java.nio.file.Path
import java.util.concurrent.TimeUnit


object OspreyService {

	const val name = "Osprey Service"
	const val defaultPort = 44342

	val version = Osprey.versionService

	fun getResourceAsStream(path: String) = OspreyService.javaClass.getResourceAsStream(path)

	fun getResourceAsString(path: String, charset: Charset = Charsets.UTF_8) =
		getResourceAsStream(path)
			?.use { stream -> stream.reader(charset).readText() }
			?: throw NoSuchElementException("resource not found: $path")

	fun getResourceAsBytes(path: String) =
		getResourceAsStream(path)
			?.use { stream -> stream.readBytes() }
			?: throw NoSuchElementException("resource not found: $path")

	val log = LoggerFactory.getLogger(OspreyService::class.java)


	interface Provider {
		fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {}
		fun registerErrors(errors: PolymorphicModuleBuilder<ErrorInfo>) {}
		fun registerService(instance: Instance, routing: Routing, prefix: String)
	}

	val services: List<Provider> = listOf(
		AboutService,
		MissingAtomsService,
		BondsService,
		ProtonationService,
		ProtonateService,
		TypesService,
		MoleculeFFInfoService,
		ForcefieldParamsService,
		MinimizeService,
		ClashesService
	)


	class Instance(
		val dir: Path,
		wait: Boolean,
		port: Int = defaultPort,
		/**
		 * If true, prefix all endpoints with "v$version" to match the effect of the Caddy reverse proxy
		 * that a production deployment of the service would use.
		 */
		useVersionPrefix: Boolean = false
	) : AutoCloseable {

		private val service = embeddedServer(Netty, port) {

			install(ContentNegotiation) {
				register(ContentType.Application.Json, SerializationConverter(json))
			}

			routing {

				// serve a simple webpage at the root
				get("/") {
					val html = getResourceAsString("index.html")
						.replace("\$name", name)
						.replace("\$version", version)
					call.respondText(html, ContentType.Text.Html)
				}

				// determine the endpoint prefix
				val prefix = if (useVersionPrefix) {
					"/v$version"
				} else {
					""
				}

				// map each service to a URL
				for (service in services) {
					service.registerService(this@Instance, this@routing, prefix)
				}
			}
		}

		init {
			service.start(wait)
		}

		override fun close() {
			service.stop(0L, 5L, TimeUnit.SECONDS)
		}
	}

	// register types for each service
	private val serializationModule = SerializersModule {

		polymorphic(ResponseInfo::class) {
			for (service in services) {
				service.registerResponses(this)
			}
		}
		polymorphic(ErrorInfo::class) {

			// register built-in types
			subclass(InternalError::class)
			subclass(RequestError::class)

			for (service in services) {
				service.registerErrors(this)
			}
		}
	}

	val json = Json {
		serializersModule = serializationModule
	}

	val serializer = KotlinxSerializer(json)
}

inline fun <reified R:ResponseInfo> Routing.service(instance: OspreyService.Instance, path: String, crossinline func: (OspreyService.Instance) -> ServiceResponse<R>) {
	get(path) {
		try {
			call.respond(func(instance))
		} catch (t: Throwable) {
			call.respondError(t)
		}
	}
}

inline fun <reified T:Any, reified R:ResponseInfo> Routing.service(instance: OspreyService.Instance, path: String, crossinline func: (OspreyService.Instance, T) -> ServiceResponse<R>) {
	post(path) {
		try {
			call.respond(func(instance, call.receive()))
		} catch (t: Throwable) {
			call.respondError(t)
		}
	}
}
