package edu.duke.cs.osprey.service

import io.ktor.application.ApplicationCall
import io.ktor.application.call
import io.ktor.features.ContentConverter
import io.ktor.features.ContentNegotiation
import io.ktor.features.suitableCharset
import io.ktor.http.ContentType
import io.ktor.http.HttpStatusCode
import io.ktor.http.content.TextContent
import io.ktor.http.withCharset
import io.ktor.request.ApplicationReceiveRequest
import io.ktor.response.respond
import io.ktor.serialization.SerializationConverter
import io.ktor.util.KtorExperimentalAPI
import io.ktor.util.pipeline.PipelineContext
import kotlinx.serialization.ImplicitReflectionSerializer
import kotlinx.serialization.KSerializer
import kotlinx.serialization.Serializable
import kotlinx.serialization.serializer
import java.util.ArrayList
import kotlin.reflect.KClass


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


class ResponseRegistrar {

	private val _responses = ArrayList<KClass<out ResponseInfo>>()
	private val _errors = ArrayList<KClass<out ErrorInfo>>()

	val responses: List<KClass<out ResponseInfo>> get() = _responses
	val errors: List<KClass<out ErrorInfo>> get() = _errors

	fun addResponse(c: KClass<out ResponseInfo>) = _responses.add(c)
	fun addError(c: KClass<out ErrorInfo>) = _errors.add(c)

	inline fun <reified R:ResponseInfo> addResponse() = addResponse(R::class)
	inline fun <reified E:ErrorInfo> addError() = addError(E::class)
}


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


@UseExperimental(ImplicitReflectionSerializer::class, KtorExperimentalAPI::class)
fun ContentNegotiation.Configuration.serializationForServiceResponse() {

	// get ktor's usual serialization converter
	val ktorConverter = SerializationConverter(OspreyService.json)

	// make a serialization converter that understands our ServiceResponse container type
	val converter = object : ContentConverter {

		override suspend fun convertForSend(context: PipelineContext<Any,ApplicationCall>, contentType: ContentType, value: Any) =
			when (value) {
				is ServiceResponse<*> -> {

					// get the nested type
					val type = when (value) {
						is ServiceResponse.Success<*> -> value.response::class
						is ServiceResponse.Failure<*> -> value.error::class
					}

					// make the serializer for the response
					val serializer = ServiceResponse.serializer(type.serializer())

					// send the usual serialized response
					@Suppress("UNCHECKED_CAST")
					val content = OspreyService.json.stringify(serializer as KSerializer<Any>, value)
					TextContent(content, contentType.withCharset(context.call.suitableCharset()))
				}
				else -> ktorConverter.convertForSend(context, contentType, value)
			}

		// for reads, just pass through to the usual serializer
		override suspend fun convertForReceive(context: PipelineContext<ApplicationReceiveRequest,ApplicationCall>) =
			ktorConverter.convertForReceive(context)
	}

	register(ContentType.Application.Json, converter)
}
