package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Parmchk
import io.ktor.routing.*
import kotlinx.serialization.Serializable
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.subclass


/* Runs Amber's parmchk2 command to generate parameters that might be missing for molecules */
object MoleculeFFInfoService : OspreyService.Provider {

	override fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {
		responses.subclass(MoleculeFFInfoResponse::class)
	}

	override fun registerErrors(errors: PolymorphicModuleBuilder<ErrorInfo>) {
		errors.subclass(MoleculeFFInfoError::class)
	}

	override fun registerService(instance: OspreyService.Instance, routing: Routing, prefix: String) {
		routing.service(instance, "$prefix/moleculeFFInfo", ::run)
	}

	fun run(instance: OspreyService.Instance, request: MoleculeFFInfoRequest): ServiceResponse<MoleculeFFInfoResponse> {

		val atomTypes = Parmchk.AtomTypes.from(request.ffname)
			?: return ServiceResponse.Success(MoleculeFFInfoResponse(null))

		val results = Parmchk.run(
			instance.dir,
			request.mol2,
			atomTypes
		)

		if (results.frcmod == null) {
			return ServiceResponse.Failure(MoleculeFFInfoError(
				"No forcefield info generate for molecule",
				results.console.joinToString("\n")
			))
		}

		return ServiceResponse.Success(MoleculeFFInfoResponse(results.frcmod))
	}
}

@Serializable
data class MoleculeFFInfoRequest(
	val mol2: String,
	val ffname: String
)

@Serializable
data class MoleculeFFInfoResponse(
	val ffinfo: String?
) : ResponseInfo

@Serializable
data class MoleculeFFInfoError(
	override val msg: String,
	val parmchkLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nParmchk:\n$parmchkLog"
}
