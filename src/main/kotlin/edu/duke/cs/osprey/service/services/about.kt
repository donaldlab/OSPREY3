package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import kotlinx.serialization.Serializable


object AboutService {

	fun registerResponses(registrar: ResponseRegistrar) {
		registrar.addResponse<AboutResponse>()
	}

	fun run(instance: OspreyService.Instance): ServiceResponse<AboutResponse> =
		ServiceResponse.Success(AboutResponse(OspreyService.name, OspreyService.version))
}

@Serializable
data class AboutResponse(
	val name: String,
	val version: String
) : ResponseInfo
