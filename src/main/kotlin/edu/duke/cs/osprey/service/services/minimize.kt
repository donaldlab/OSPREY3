package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Sander
import io.ktor.routing.*
import kotlinx.serialization.Serializable
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.subclass


object MinimizeService : OspreyService.Provider {

	override fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {
		responses.subclass(MinimizeResponse::class)
	}

	override fun registerErrors(errors: PolymorphicModuleBuilder<ErrorInfo>) {
		errors.subclass(MinimizeError::class)
	}

	override fun registerService(instance: OspreyService.Instance, routing: Routing, prefix: String) {
		routing.service(instance, "$prefix/minimize", ::run)
	}

	fun run(instance: OspreyService.Instance, request: MinimizeRequest): ServiceResponse<MinimizeResponse> {

		val commands = ArrayList<String>()

		// TODO: report progress info to the caller somehow
		val reportEveryCycles = 10

		commands += listOf(
			"imin=1", // do cartesian minimization
			"maxcyc=${request.numCycles}",
			"ntpr=$reportEveryCycles",
			"ntxo=1" // format the output coords in ASCII
		)

		// TODO: expose more options?
		// ntmin     minimization type

		if (request.restraintMask != null) {

			// alas, the restraint mask has a limited size
			if (request.restraintMask.length > Sander.maxRestraintMaskSize) {
				return ServiceResponse.Failure(RequestError(
					"Alas, the restraintmask for sander can only be ${Sander.maxRestraintMaskSize} characters." +
					" ${request.restraintMask.length} characters is too long."
				))
			}

			commands += listOf(
				"ntr=1", // turn on cartesian restraints
				"restraint_wt=${request.restraintWeight}",
				"restraintmask='${Sander.sanitizeRestraintMask(request.restraintMask)}'"
			)
		}

		commands += listOf(
			"igb=1" // use generalized borne solvation
		)
		// TODO: expose more solvation options?

		val results = Sander.run(
			instance.dir,
			request.top,
			request.crd,
			"""
				|Header
				|&cntrl
				|${commands.joinToString(",\n")}
				|/
			""".trimMargin()
		)

		val coords = results.coords
			?: return ServiceResponse.Failure(MinimizeError(
				"Sander didn't produce output coordinates",
				results.console.joinToString("\n")
			))

		return ServiceResponse.Success(MinimizeResponse(coords))
	}
}

@Serializable
data class MinimizeRequest(
	val top: String,
	val crd: String,
	val numCycles: Int,
	val restraintMask: String? = null,
	val restraintWeight: Double = 1.0
)

@Serializable
data class MinimizeResponse(
	val coords: List<Point3d>
) : ResponseInfo

@Serializable
data class MinimizeError(
	override val msg: String,
	val sanderLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nSander:\n$sanderLog"
}
