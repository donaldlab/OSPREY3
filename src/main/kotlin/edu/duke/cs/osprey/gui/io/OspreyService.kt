package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.service.Point3d
import edu.duke.cs.osprey.service.ServiceResponse
import edu.duke.cs.osprey.service.services.*
import edu.duke.cs.osprey.service.OspreyService as Server
import io.ktor.client.HttpClient
import io.ktor.client.call.TypeInfo
import io.ktor.client.engine.cio.CIO
import io.ktor.client.features.json.JsonFeature
import io.ktor.client.features.json.JsonSerializer
import io.ktor.client.features.json.serializer.KotlinxSerializer
import io.ktor.client.request.HttpRequestBuilder
import io.ktor.client.request.request
import io.ktor.http.ContentType
import io.ktor.http.HttpMethod
import io.ktor.http.contentType
import io.ktor.util.KtorExperimentalAPI
import io.ktor.utils.io.core.Input
import io.ktor.utils.io.core.readText
import kotlinx.serialization.ImplicitReflectionSerializer
import kotlinx.serialization.serializer
import org.joml.Vector3d
import kotlin.reflect.jvm.jvmErasure


object OspreyService {

	@UseExperimental(KtorExperimentalAPI::class)
	private val client =
		HttpClient(CIO) {
			install(JsonFeature) {
				serializer = ServiceSerializer()
			}
			engine {
				// some requests can take a long time, so disable timeouts
				requestTimeout = 0L
			}
		}

	private fun HttpRequestBuilder.path(path: String) {
		url {
			val provider = UserSettings.serviceProvider
			host = provider.hostname
			port = provider.port
			path(path)
		}
	}

	private fun HttpRequestBuilder.get(path: String) {
		path(path)
		method = HttpMethod.Get
	}

	private fun HttpRequestBuilder.post(path: String, obj: Any) {
		path(path)
		method = HttpMethod.Post
		contentType(ContentType.Application.Json)
		body = obj
	}

	// NOTE: tragically, we can't remove some code duplication here due to a compiler bug
	// see: https://youtrack.jetbrains.com/issue/KT-34051
	// maybe someday we can?

	suspend fun about() =
		client.request<ServiceResponse<AboutResponse>> {
			get("about")
		}
		.responseOrThrow()

	suspend fun missingAtoms(request: MissingAtomsRequest) =
		client.request<ServiceResponse<MissingAtomsResponse>> {
			post("missingAtoms", request)
		}
		.responseOrThrow()

	suspend fun bonds(request: BondsRequest) =
		client.request<ServiceResponse<BondsResponse>> {
			post("bonds", request)
		}
		.responseOrThrow()

	suspend fun protonation(request: ProtonationRequest) =
		client.request<ServiceResponse<ProtonationResponse>> {
			post("protonation", request)
		}
		.responseOrThrow()

	suspend fun protonate(request: ProtonateRequest) =
		client.request<ServiceResponse<ProtonateResponse>> {
			post("protonate", request)
		}
		.responseOrThrow()

	suspend fun types(request: TypesRequest) =
		client.request<ServiceResponse<TypesResponse>> {
			post("types", request)
		}
		.responseOrThrow()

	suspend fun moleculeFFInfo(request: MoleculeFFInfoRequest) =
		client.request<ServiceResponse<MoleculeFFInfoResponse>> {
			post("moleculeFFInfo", request)
		}
		.responseOrThrow()

	suspend fun forcefieldParams(request: ForcefieldParamsRequest) =
		client.request<ServiceResponse<ForcefieldParamsResponse>> {
			post("forcefieldParams", request)
		}
		.responseOrThrow()

	suspend fun minimize(request: MinimizeRequest) =
		client.request<ServiceResponse<MinimizeResponse>> {
			post("minimize", request)
		}
		.responseOrThrow()

	suspend fun clashes(request: ClashesRequest) =
		client.request<ServiceResponse<ClashesResponse>> {
			post("clashes", request)
		}
		.responseOrThrow()
}

private class ServiceSerializer : JsonSerializer {

	// get ktor's usual serializer
	private val ktorSerializer = KotlinxSerializer()

	// for writes, just pass through to the usual serializer
	override fun write(data: Any, contentType: ContentType) =
		ktorSerializer.write(data, contentType)

	// for reads, look for our ServiceResponse type and handle it specially
	@UseExperimental(ImplicitReflectionSerializer::class)
	override fun read(type: TypeInfo, body: Input) =
		when (type.type) {

			ServiceResponse::class -> {

				// get the response type from the given type info
				val rtype = type.kotlinType!!.arguments[0].type!!.jvmErasure

				// read the response from the server
				val text = body.readText()

				// de-serialize the json
				try {
					Server.json.parse(ServiceResponse.serializer(rtype.serializer()), text)
				} catch (t: Throwable) {
					// super rare concurrency bug here... can't reproduce reliably
					// best I can do for now is add more info to the error log
					throw IllegalArgumentException("""
						|can't deserialize JSON
						|read ${text.length} chars, end of input? ${body.endOfInput}
						|Response from server:
						|
						|$text
					""".trimMargin(), t)
				}
			}

			// otherwise, pass on to the usual serializer
			else -> ktorSerializer.read(type, body)
		}
}

fun Point3d.toVector3d() = Vector3d(x, y, z)
