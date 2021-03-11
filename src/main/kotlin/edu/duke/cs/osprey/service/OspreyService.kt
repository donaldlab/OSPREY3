package edu.duke.cs.osprey.service

import edu.duke.cs.osprey.Osprey
import edu.duke.cs.osprey.service.services.*
import io.ktor.application.call
import io.ktor.application.install
import io.ktor.features.ContentNegotiation
import io.ktor.http.ContentType
import io.ktor.request.receive
import io.ktor.response.respond
import io.ktor.response.respondText
import io.ktor.routing.Routing
import io.ktor.routing.get
import io.ktor.routing.post
import io.ktor.routing.routing
import io.ktor.server.engine.embeddedServer
import io.ktor.server.engine.stop
import io.ktor.server.netty.Netty
import kotlinx.serialization.ImplicitReflectionSerializer
import kotlinx.serialization.json.Json
import kotlinx.serialization.json.JsonConfiguration
import kotlinx.serialization.modules.SerializersModule
import kotlinx.serialization.serializer
import org.slf4j.LoggerFactory
import java.nio.charset.Charset
import java.nio.file.Path
import java.util.concurrent.TimeUnit
import kotlin.reflect.KClass


object OspreyService {

	const val name = "Osprey Service"

	val version = Osprey.versionService

	fun getResourceAsStream(path: String) = OspreyService.javaClass.getResourceAsStream(path)

	fun getResourceAsString(path: String, charset: Charset = Charsets.UTF_8) =
		getResourceAsStream(path).use { stream -> stream.reader(charset).readText() }

	fun getResourceAsBytes(path: String) =
		getResourceAsStream(path).use { stream -> stream.readBytes() }

	val log = LoggerFactory.getLogger(OspreyService::class.java)


	class Instance(val dir: Path, wait: Boolean, port: Int = 8080) : AutoCloseable {

		private val service = embeddedServer(Netty, port) {

			install(ContentNegotiation) {
				serializationForServiceResponse()
			}

			routing {

				// serve a simple webpage at the root
				get("/") {
					val html = getResourceAsString("index.html")
						.replace("\$name", name)
						.replace("\$version", version)
					call.respondText(html, ContentType.Text.Html)
				}

				// map each service to a URL
				service("/about", AboutService::run)
				service("/missingAtoms", MissingAtomsService::run)
				service("/bonds", BondsService::run)
				service("/protonation", ProtonationService::run)
				service("/protonate", ProtonateService::run)
				service("/types", TypesService::run)
				service("/moleculeFFInfo", MoleculeFFInfoService::run)
				service("/forcefieldParams", ForcefieldParamsService::run)
				service("/minimize", MinimizeService::run)
				service("/clashes", ClashesService::run)
			}
		}

		init {
			service.start(wait)
		}

		override fun close() {
			service.stop(0L, 5L, TimeUnit.SECONDS)
		}

		private inline fun <reified R:ResponseInfo> Routing.service(path: String, crossinline func: (Instance) -> ServiceResponse<R>) {
			get(path) {
				try {
					call.respond(func(this@Instance))
				} catch (t: Throwable) {
					call.respondError(t)
				}
			}
		}

		private inline fun <reified T:Any, reified R:ResponseInfo> Routing.service(path: String, crossinline func: (Instance, T) -> ServiceResponse<R>) {
			post(path) {
				try {
					call.respond(func(this@Instance, call.receive()))
				} catch (t: Throwable) {
					call.respondError(t)
				}
			}
		}
	}

	// register types for each service
	@UseExperimental(ImplicitReflectionSerializer::class)
	val serializationModule = SerializersModule {

		val registrar = ResponseRegistrar()

		// register built-in types
		registrar.addError<InternalError>()
		registrar.addError<RequestError>()

		// ask each service to register their responses and errors
		AboutService.registerResponses(registrar)
		MissingAtomsService.registerResponses(registrar)
		BondsService.registerResponses(registrar)
		ProtonationService.registerResponses(registrar)
		ProtonateService.registerResponses(registrar)
		TypesService.registerResponses(registrar)
		MoleculeFFInfoService.registerResponses(registrar)
		ForcefieldParamsService.registerResponses(registrar)
		MinimizeService.registerResponses(registrar)
		ClashesService.registerResponses(registrar)

		polymorphic<ResponseInfo> {
			for (response in registrar.responses) {
				@Suppress("UNCHECKED_CAST")
				val c = response as KClass<ResponseInfo>
				addSubclass(c, c.serializer())
			}
		}
		polymorphic<ErrorInfo> {
			for (error in registrar.errors) {
				@Suppress("UNCHECKED_CAST")
				val c = error as KClass<ErrorInfo>
				addSubclass(c, c.serializer())
			}
		}
	}

	val json = Json(
		configuration = JsonConfiguration.Stable.copy(
			encodeDefaults = true,
			strictMode = false,
			unquoted = false,
			prettyPrint = false,
			useArrayPolymorphism = true
		),
		context = serializationModule
	)
}
