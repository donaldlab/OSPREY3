package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Antechamber
import edu.duke.cs.osprey.service.amber.Leap
import kotlinx.serialization.Serializable


object ProtonationService {

	fun registerResponses(registrar: ResponseRegistrar) {
		registrar.addResponse<ProtonationResponse>()
		registrar.addError<ProtonationLeapError>()
		registrar.addError<ProtonationAntechamberError>()
	}

	fun run(instance: OspreyService.Instance, request: ProtonationRequest): ServiceResponse<ProtonationResponse> {

		// get the atom types, if any
		val atomTypes = request.atomTypes?.let { id ->
			Antechamber.AtomTypes.values().find { it.id == id }
				?: return ServiceResponse.Failure(RequestError("unrecognized atom types: $id"))
		}

		if (atomTypes == null) {

			// run LEaP to add the hydrogens
			val results = Leap.run(
				instance.dir,
				filesToWrite = mapOf("in.pdb" to request.pdb),
				commands = """
					|verbosity 2
					|source leaprc.${Leap.sanitizeToken(request.ffname)}
					|mol = loadPDB in.pdb
					|addH mol
					|saveMol2 mol out.mol2 0
				""".trimMargin(),
				filesToRead = listOf("out.mol2")
			)

			val mol2 = results.files["out.mol2"]
				?: return ServiceResponse.Failure(ProtonationLeapError(
					"LEaP didn't produce an output molecule",
					results.console.joinToString("\n")
				))

			return ServiceResponse.Success(ProtonationResponse(mol2))

		} else {

			// run antechamber to infer all the atom and bond types
			val antechamberResults = Antechamber.run(
				instance.dir,
				request.pdb,
				Antechamber.InType.Pdb,
				atomTypes
			)

			val atomTypesMol2 = antechamberResults.mol2
				?: return ServiceResponse.Failure(ProtonationAntechamberError(
					"Antechamber didn't produce an output molecule",
					antechamberResults.console.joinToString("\n")
				))

			// run LEaP to add the hydrogens
			val leapResults = Leap.run(
				instance.dir,
				filesToWrite = mapOf("in.mol2" to atomTypesMol2),
				commands = """
					|verbosity 2
					|source leaprc.${Leap.sanitizeToken(request.ffname)}
					|mol = loadMol2 in.mol2
					|addH mol
					|saveMol2 mol out.mol2 0
				""".trimMargin(),
				filesToRead = listOf("out.mol2")
			)

			val protonatedMol2 = leapResults.files["out.mol2"]
				?: return ServiceResponse.Failure(ProtonationLeapError(
					"LEaP didn't produce an output molecule",
					leapResults.console.joinToString("\n")
				))

			return ServiceResponse.Success(ProtonationResponse(protonatedMol2))
		}
	}
}


@Serializable
data class ProtonationRequest(
	val pdb: String,
	val ffname: String,
	val atomTypes: String?
)

@Serializable
data class ProtonationResponse(
	val mol2: String
) : ResponseInfo

@Serializable
data class ProtonationLeapError(
	override val msg: String,
	val leapLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nLEaP:\n$leapLog"
}

@Serializable
data class ProtonationAntechamberError(
	override val msg: String,
	val antechamberLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nAntechamber:\n$antechamberLog"
}
