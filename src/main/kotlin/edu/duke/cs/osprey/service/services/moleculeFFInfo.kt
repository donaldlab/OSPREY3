package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Parmchk
import kotlinx.serialization.Serializable


object MoleculeFFInfoService {

	fun registerResponses(registrar: ResponseRegistrar) {
		registrar.addResponse<MoleculeFFInfoResponse>()
		registrar.addError<MoleculeFFInfoError>()
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
				"No forecfield info generate for molecule",
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
