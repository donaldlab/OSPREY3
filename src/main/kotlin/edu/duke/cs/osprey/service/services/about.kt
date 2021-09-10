package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import io.ktor.routing.*
import kotlinx.serialization.Serializable
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.subclass


object AboutService : OspreyService.Provider {

	override fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {
		responses.subclass(AboutResponse::class)
	}

	override fun registerService(instance: OspreyService.Instance, routing: Routing, prefix: String) {
		routing.service(instance, "$prefix/about", ::run)
	}

	fun run(instance: OspreyService.Instance): ServiceResponse<AboutResponse> =
		ServiceResponse.Success(AboutResponse(OspreyService.name, OspreyService.version))
}

@Serializable
data class AboutResponse(
	val name: String,
	val version: String
) : ResponseInfo
