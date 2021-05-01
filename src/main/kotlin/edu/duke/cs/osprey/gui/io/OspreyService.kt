package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.service.Point3d
import edu.duke.cs.osprey.service.ResponseInfo
import edu.duke.cs.osprey.service.ServiceResponse
import edu.duke.cs.osprey.service.services.*
import io.ktor.client.HttpClient
import io.ktor.client.engine.cio.CIO
import io.ktor.client.features.json.JsonFeature
import io.ktor.client.request.HttpRequestBuilder
import io.ktor.client.request.request
import io.ktor.http.ContentType
import io.ktor.http.HttpMethod
import io.ktor.http.contentType
import org.joml.Vector3d


object OspreyService {

	private val client =
		HttpClient(CIO) {
			install(JsonFeature) {
				serializer = edu.duke.cs.osprey.service.OspreyService.serializer
			}
			engine {
				// some requests can take a long time, so disable timeouts
				requestTimeout = Long.MAX_VALUE
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

	private suspend inline fun <reified R:ResponseInfo> get(path: String) =
		client.request<ServiceResponse<R>> {
			path(path)
			method = HttpMethod.Get
		}
		.responseOrThrow()

	private suspend inline fun <reified R:ResponseInfo> post(path: String, arg: Any) =
		client.request<ServiceResponse<R>> {
			path(path)
			method = HttpMethod.Post
			contentType(ContentType.Application.Json)
			body = arg
		}
		.responseOrThrow()


	suspend fun about() =
		get<AboutResponse>("about")

	suspend fun missingAtoms(request: MissingAtomsRequest) =
		post<MissingAtomsResponse>("missingAtoms", request)

	suspend fun bonds(request: BondsRequest) =
		post<BondsResponse>("bonds", request)

	suspend fun protonation(request: ProtonationRequest) =
		post<ProtonationResponse>("protonation", request)

	suspend fun protonate(request: ProtonateRequest) =
		post<ProtonateResponse>("protonate", request)

	suspend fun types(request: TypesRequest) =
		post<TypesResponse>("types", request)

	suspend fun moleculeFFInfo(request: MoleculeFFInfoRequest) =
		post<MoleculeFFInfoResponse>("moleculeFFInfo", request)

	suspend fun forcefieldParams(request: ForcefieldParamsRequest) =
		post<ForcefieldParamsResponse>("forcefieldParams", request)

	suspend fun minimize(request: MinimizeRequest) =
		post<MinimizeResponse>("minimize", request)

	suspend fun clashes(request: ClashesRequest) =
		post<ClashesResponse>("clashes", request)
}


fun Point3d.toVector3d() = Vector3d(x, y, z)
