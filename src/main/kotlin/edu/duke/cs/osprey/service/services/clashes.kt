package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.clashes.Probe
import io.ktor.routing.*
import kotlinx.serialization.Serializable
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.subclass


object ClashesService : OspreyService.Provider {

	override fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {
		responses.subclass(ClashesResponse::class)
	}

	override fun registerErrors(errors: PolymorphicModuleBuilder<ErrorInfo>) {
		errors.subclass(ClashesError::class)
	}

	override fun registerService(instance: OspreyService.Instance, routing: Routing) {
		routing.service(instance, "/clashes", ::run)
	}

	fun run(instance: OspreyService.Instance, request: ClashesRequest): ServiceResponse<ClashesResponse> {

		val results = Probe.run(
			instance.dir,
			request.pdb
		)

		return ServiceResponse.Success(ClashesResponse(
			results.groups.mapValues { (id, group) ->
				ClashesResponse.Group(
					id,
					dots = group.dots,
					vectors = group.vectors
				)
			}
		))
	}
}

@Serializable
data class ClashesRequest(
	val pdb: String
)

@Serializable
data class ClashesResponse(
	val groups: Map<String,Group>
) : ResponseInfo {

	@Serializable
	data class Group(
		val id: String,
		val dots: Map<String,List<Point3d>>,
		val vectors: Map<String,List<Pair<Point3d,Point3d>>>
	)
}

@Serializable
data class ClashesError(
	override val msg: String,
	val probeLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nProbe:\n$probeLog"
}
