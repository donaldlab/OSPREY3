package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Antechamber
import edu.duke.cs.osprey.service.amber.Leap
import io.ktor.routing.*
import kotlinx.serialization.Serializable
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.subclass


/* Use LEaP/Antechamber to generate a MOL2 structure, which contains atom connectivity */
object BondsService : OspreyService.Provider {

	override fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {
		responses.subclass(BondsResponse::class)
	}

	override fun registerErrors(errors: PolymorphicModuleBuilder<ErrorInfo>) {
		errors.subclass(BondsLeapError::class)
		errors.subclass(BondsAntechamberError::class)
	}

	override fun registerService(instance: OspreyService.Instance, routing: Routing, prefix: String) {
		routing.service(instance, "$prefix/bonds", ::run)
	}

	fun run(instance: OspreyService.Instance, request: BondsRequest): ServiceResponse<BondsResponse> {

		// if we got a forcefield name, use LEaP
		val ffname = request.ffname
		if (ffname != null) {

			// run LEaP to infer all the missing atoms
			val results = Leap.run(
				instance.dir,
				filesToWrite = mapOf("in.pdb" to request.pdb),
				commands = """
					|verbosity 2
					|source leaprc.${Leap.sanitizeToken(ffname)}
					|mol = loadPDB in.pdb
					|saveMol2 mol out.mol2 0
				""".trimMargin(),
				filesToRead = listOf("out.mol2")
			)

			val mol2 = results.files["out.mol2"]
				?: return ServiceResponse.Failure(BondsLeapError(
					"LEaP didn't produce an output molecule",
					results.console.joinToString("\n")
				))

			return ServiceResponse.Success(BondsResponse(mol2))

		// we didn't get a forcefield, so use antechamber
		} else {

			val results = Antechamber.run(
				instance.dir,
				request.pdb,
				Antechamber.InType.Pdb,
				Antechamber.AtomTypes.SYBYL
			)

			val mol2 = results.mol2
				?: return ServiceResponse.Failure(BondsAntechamberError(
					"Antechamber didn't produce an output molecule",
					results.console.joinToString("\n")
				))

			return ServiceResponse.Success(BondsResponse(mol2))
		}
	}
}


@Serializable
data class BondsRequest(
	val pdb: String,
	val ffname: String?
)

@Serializable
data class BondsResponse(
	val mol2: String
) : ResponseInfo

@Serializable
data class BondsLeapError(
	override val msg: String,
	val leapLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nLEaP:\n$leapLog"
}

@Serializable
data class BondsAntechamberError(
	override val msg: String,
	val antechamberLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nAntechamber:\n$antechamberLog"
}
