package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Leap
import kotlinx.serialization.Serializable


object ProtonateService {

	fun registerResponses(registrar: ResponseRegistrar) {
		registrar.addResponse<ProtonateResponse>()
		registrar.addError<ProtonateError>()
	}

	fun run(instance: OspreyService.Instance, request: ProtonateRequest): ServiceResponse<ProtonateResponse> {

		fun String.sanitized() = Leap.sanitizeToken(this)

		// build the LEaP commands
		val commands = ArrayList<String>().apply {

			// show all the info in the console
			add("verbosity 2")

			// load the generalized amber forcefield (for small molecules)
			// TODO: let caller pick the forcefield? (will have to have AmberTypes for each possible choice)
			add("source leaprc.gaff2")

			// read our molecule fragment
			add("mol = loadMol2 in.mol2")

			// set the central atom type
			add("set mol.1.${request.atomName.sanitized()} type ${request.atomType.sanitized()}")

			// add the bonds
			for (bond in request.bonds) {
				add("set mol.1.${bond.atomName.sanitized()} type ${bond.atomType.sanitized()}")
				add("deleteBond mol.1.${request.atomName.sanitized()} mol.1.${bond.atomName.sanitized()}")
				add("bond mol.1.${request.atomName.sanitized()} mol.1.${bond.atomName.sanitized()} ${bond.bondType.sanitized()}")
			}

			// add the hydrogens
			for ((i, h) in request.hydrogens.withIndex()) {
				add("h$i = createAtom ${h.atomName.sanitized()} ${h.atomType.sanitized()} 0.0")
				add("set h$i element H")
				add("add mol.1 h$i")
				add("bond mol.1.${request.atomName.sanitized()} h$i")
				add("select h$i")
			}
			add("rebuildSelectedAtoms mol")

			// save the built molecule
			add("saveMol2 mol out.mol2 0")
		}

		// run LEaP
		val results = Leap.run(
			instance.dir,
			filesToWrite = mapOf(
				"in.mol2" to request.mol2
			),
			commands = commands.joinToString("\n"),
			filesToRead = listOf("out.mol2"),
			// use the debug build of teLeap, so we get more info in the console
			debugFiles = listOf("model.c")
		)

		val protonatedMol2 = results.files["out.mol2"]
			?: return ServiceResponse.Failure(ProtonateError(
				"LEaP didn't produce an output molecule",
				request.mol2,
				results.console.joinToString("\n")
			))

		return ServiceResponse.Success(ProtonateResponse(protonatedMol2))
	}
}

@Serializable
data class ProtonateRequest(
	val mol2: String,
	val atomName: String,
	val atomType: String,
	val bonds: List<Bond>,
	val hydrogens: List<Hydrogen>
) {

	@Serializable
	data class Bond(
		val atomName: String,
		val atomType: String,
		val bondType: String
	)

	@Serializable
	data class Hydrogen(
		val atomName: String,
		val atomType: String
	)
}

@Serializable
data class ProtonateResponse(
	val mol2: String
) : ResponseInfo

@Serializable
data class ProtonateError(
	override val msg: String,
	val mol2: String,
	val leapLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nMol2:\n$mol2\n\nLEaP:\n$leapLog"
}
