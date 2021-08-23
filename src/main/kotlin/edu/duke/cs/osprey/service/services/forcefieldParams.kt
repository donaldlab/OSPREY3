package edu.duke.cs.osprey.service.services

import edu.duke.cs.osprey.service.*
import edu.duke.cs.osprey.service.amber.Leap
import io.ktor.routing.*
import kotlinx.serialization.Serializable
import kotlinx.serialization.modules.PolymorphicModuleBuilder
import kotlinx.serialization.modules.subclass
import java.util.ArrayList
import java.util.HashMap


object ForcefieldParamsService : OspreyService.Provider {

	override fun registerResponses(responses: PolymorphicModuleBuilder<ResponseInfo>) {
		responses.subclass(ForcefieldParamsResponse::class)
	}

	override fun registerErrors(errors: PolymorphicModuleBuilder<ErrorInfo>) {
		errors.subclass(ForcefieldParamsError::class)
	}

	override fun registerService(instance: OspreyService.Instance, routing: Routing, prefix: String) {
		routing.service(instance, "$prefix/forcefieldParams", ::run)
	}

	fun run(instance: OspreyService.Instance, request: ForcefieldParamsRequest): ServiceResponse<ForcefieldParamsResponse> {

		fun String.sanitize() = Leap.sanitizeToken(this)

		val molname = "mol.%d.mol2"
		val frcname = "mol.%d.%d.frc"
		val topname = "mol.top"
		val crdname = "mol.crd"

		val commands = ArrayList<String>()
		commands += "verbosity 2"

		// write each unique forcefield name only once
		request.molecules
			.map { it.ffname }
			.toSet()
			.forEach { ffname ->
				commands += "source leaprc.${ffname.sanitize()}"
			}

		for ((moli, mol) in request.molecules.withIndex()) {
			for (ffinfoi in mol.ffinfos.indices) {
				commands += "loadAmberParams ${frcname.format(moli, ffinfoi)}"
			}
			commands += "mol$moli = loadMol2 ${molname.format(moli)}"
		}
		commands += if (request.molecules.size > 1) {
			"combined = combine { ${request.molecules.indices.joinToString(" ") { "mol$it" }} }"
		} else {
			"combined = mol0"
		}
		commands += "saveAmberParm combined $topname, $crdname"

		// run LEaP to get the params
		val results = Leap.run(
			instance.dir,
			filesToWrite = HashMap<String,String>().apply {
				for ((moli, mol) in request.molecules.withIndex()) {
					put(molname.format(moli), mol.mol2)
					for ((ffinfoi, ffinfo) in mol.ffinfos.withIndex()) {
						put(frcname.format(moli, ffinfoi), ffinfo)
					}
				}
			},
			commands = commands.joinToString("\n"),
			filesToRead = listOf(topname, crdname)
		)

		val top = results.files[topname]
			?: return ServiceResponse.Failure(ForcefieldParamsError(
				"no topology file",
				results.console.joinToString("\n")
			))

		val crd = results.files[crdname]
			?: return ServiceResponse.Failure(ForcefieldParamsError(
				"no coordinates file",
				results.console.joinToString("\n")
			))

		return ServiceResponse.Success(ForcefieldParamsResponse(top, crd))
	}
}

@Serializable
data class ForcefieldParamsRequest(
	val molecules: List<MolInfo>
) {

	constructor(vararg molecules: MolInfo) : this(molecules.toList())

	@Serializable
	data class MolInfo(
		val mol2: String,
		val ffname: String,
		val ffinfos: List<String>
	)
}

@Serializable
data class ForcefieldParamsResponse(
	val params: String,
	val coords: String
) : ResponseInfo

@Serializable
data class ForcefieldParamsError(
	override val msg: String,
	val leapLog: String
) : ErrorInfo {
	override fun message() = "$msg\n\nLEaP:\n$leapLog"
}
