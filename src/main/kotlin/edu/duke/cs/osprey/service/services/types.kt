package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Antechamber
import edu.duke.cs.osprey.service.amber.Leap
import edu.duke.cs.osprey.service.amber.SQM
import kotlinx.serialization.Serializable


object TypesService {

	fun registerResponses(registrar: ResponseRegistrar) {
		registrar.addResponse<TypesResponse>()
		registrar.addError<TypesLeapError>()
		registrar.addError<TypesAntechamberError>()
	}

	fun run(instance: OspreyService.Instance, request: TypesRequest): ServiceResponse<TypesResponse> {

		val molSettings = request.molSettings
		if (molSettings != null) {

			val results = Leap.run(
				instance.dir,
				filesToWrite = mapOf("in.pdb" to molSettings.pdb),
				commands = """
					|verbosity 2
					|source leaprc.${Leap.sanitizeToken(molSettings.ffname)}
					|mol = loadPDB in.pdb
					|saveMol2 mol out.mol2 1
				""".trimMargin(),
				filesToRead = listOf("out.mol2")
			)

			val mol2 = results.files["out.mol2"]
				?: return ServiceResponse.Failure(TypesLeapError(
					"LEaP didn't produce an output molecule",
					results.console.joinToString("\n")
				))

			return ServiceResponse.Success(TypesResponse(mol2))
		}

		val smallMolSettings = request.smallMolSettings
		if (smallMolSettings != null) {

			// get the atom types
			val atomTypes = Antechamber.AtomTypes.values()
				.find { it.id == smallMolSettings.atomTypes }
				?: return ServiceResponse.Failure(RequestError("unrecognized atom types: ${smallMolSettings.atomTypes}"))

			// get the charge settings, if any
			val chargeSettings = smallMolSettings.chargeSettings
			val results = if (chargeSettings != null) {

				// get the charge method
				val chargeMethod = Antechamber.ChargeMethod.values()
					.find { it.id == chargeSettings.chargeMethod }
					?: return ServiceResponse.Failure(RequestError("unrecognized charge method: ${chargeSettings.chargeMethod}"))

				Antechamber.run(
					instance.dir,
					smallMolSettings.mol2,
					Antechamber.InType.Mol2,
					atomTypes,
					useACDoctor = false,
					generateCharges = chargeMethod,
					netCharge = chargeSettings.netCharge,
					sqmOptions = SQM.Options(
						maxCyc = chargeSettings.numMinimizationSteps,
						ntpr = 1 // print out progress info as frequently as possible
					)
				)

			} else {

				// don't generate charges, just get atom types
				Antechamber.run(
					instance.dir,
					smallMolSettings.mol2,
					Antechamber.InType.Mol2,
					atomTypes,
					useACDoctor = false
				)
			}

			val mol2 = results.mol2
				?: return ServiceResponse.Failure(TypesAntechamberError(
					"Antechamber didn't produce an output molecule",
					results.console.joinToString("\n")
				))

			return ServiceResponse.Success(TypesResponse(mol2))
		}

		return ServiceResponse.Failure(RequestError("no settings given"))
	}
}

@Serializable
data class TypesRequest(
	val molSettings: MoleculeSettings? = null,
	val smallMolSettings: SmallMoleculeSettings? = null
) {

	@Serializable
	data class MoleculeSettings(
		val pdb: String,
		val ffname: String
	) {
		fun toRequest() =
			TypesRequest(molSettings = this)
	}

	@Serializable
	data class SmallMoleculeSettings(
		val mol2: String,
		val atomTypes: String,
		val chargeSettings: ChargeSettings? = null
	) {
		fun toRequest() =
			TypesRequest(smallMolSettings = this)
	}

	@Serializable
	data class ChargeSettings(
		val chargeMethod: String,
		val netCharge: Int,
		val numMinimizationSteps: Int
	)
}

@Serializable
data class TypesResponse(
	val mol2: String
) : ResponseInfo

@Serializable
data class TypesLeapError(
	override val msg: String,
	val leapLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nLEaP:\n$leapLog"
}

@Serializable
data class TypesAntechamberError(
	override val msg: String,
	val antechamberLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nAntechamber:\n$antechamberLog"
}
