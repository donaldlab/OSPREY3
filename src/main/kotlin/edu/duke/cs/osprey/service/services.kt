package edu.duke.cs.osprey.service

import io.ktor.application.ApplicationCall
import io.ktor.http.HttpStatusCode
import io.ktor.response.respond
import kotlinx.serialization.Serializable


@Serializable
sealed class ServiceResponse<R:ResponseInfo> {

	@Serializable
	data class Success<R:ResponseInfo>(
		val response: R
	) : ServiceResponse<R>()

	@Serializable
	data class Failure<R:ResponseInfo>(
		val error: ErrorInfo
	) : ServiceResponse<R>()

	fun responseOrThrow(): R =
		when (this) {
			is Success -> response
			is Failure -> throw ServiceException(error)
		}
}

interface ResponseInfo

interface ErrorInfo {
	/** a short message quickly describing the error */
	val msg: String
	/** show all of the error information as a string */
	fun message(): String = msg
}

/** for when something unexpected happened in the service */
@Serializable
class InternalError(override val msg: String) : ErrorInfo

/** for when the request is bad */
@Serializable
class RequestError(override val msg: String) : ErrorInfo

class ServiceException(val error: ErrorInfo) : RuntimeException(error.message())


suspend fun ApplicationCall.respond(response: ServiceResponse<*>) {

	// figure out the status code based on the response
	val status = when (response) {
		is ServiceResponse.Success -> HttpStatusCode.OK
		is ServiceResponse.Failure ->
			when (response.error) {
				is RequestError -> HttpStatusCode.BadRequest
				else -> HttpStatusCode.InternalServerError
			}
	}

	respond(status, response)
}

suspend fun ApplicationCall.respondError(t: Throwable) {

	// log the full error information for debugging
	OspreyService.log.error("An error occurred while executing a service", t)

	// but don't leak exception info to the world
	// send back a sanitized error message instead
	respond(
		status = HttpStatusCode.InternalServerError,
		message = ServiceResponse.Failure<ResponseInfo>(InternalError("Internal server error"))
	)
}
