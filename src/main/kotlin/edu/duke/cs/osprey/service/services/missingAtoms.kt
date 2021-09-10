package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Leap
import io.ktor.routing.*
import kotlinx.serialization.Serializable
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.subclass


object MissingAtomsService : OspreyService.Provider {

	override fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {
		responses.subclass(MissingAtomsResponse::class)
	}

	override fun registerErrors(errors: PolymorphicModuleBuilder<ErrorInfo>) {
		errors.subclass(MissingAtomsError::class)
	}

	override fun registerService(instance: OspreyService.Instance, routing: Routing, prefix: String) {
		routing.service(instance, "$prefix/missingAtoms", ::run)
	}

	fun run(instance: OspreyService.Instance, request: MissingAtomsRequest): ServiceResponse<MissingAtomsResponse> {

		// run LEaP to infer all the missing atoms
		val results = Leap.run(
			instance.dir,
			filesToWrite = mapOf("in.pdb" to request.pdb),
			commands = """
				|verbosity 2
				|source leaprc.${Leap.sanitizeToken(request.ffname)}
				|mol = loadPDB in.pdb
				|saveMol2 mol out.mol2 0
			""".trimMargin(),
			filesToRead = listOf("out.mol2")
		)

		val mol2 = results.files["out.mol2"]
			?: return ServiceResponse.Failure(MissingAtomsError(
				"LEaP didn't produce an output molecule",
				results.console.joinToString("\n")
			))

		return ServiceResponse.Success(MissingAtomsResponse(mol2))
	}
}

@Serializable
data class MissingAtomsRequest(
	val pdb: String,
	val ffname: String
)

@Serializable
data class MissingAtomsResponse(
	val mol2: String
) : ResponseInfo

@Serializable
data class MissingAtomsError(
	override val msg: String,
	val leapLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nLEaP:\n$leapLog"
}
