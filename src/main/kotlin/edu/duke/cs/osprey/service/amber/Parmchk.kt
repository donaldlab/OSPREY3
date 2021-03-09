package edu.duke.cs.osprey.service.amber

import edu.duke.cs.osprey.service.read
import edu.duke.cs.osprey.service.stream
import edu.duke.cs.osprey.service.tempFolder
import edu.duke.cs.osprey.service.write
import java.io.IOException
import java.nio.file.Path


object Parmchk {

	data class Results(
		val errorCode: Int,
		val console: List<String>,
		val frcmod: String?
	)

	enum class AtomTypes(val id: String) {

		Gaff("gaff"),
		Gaff2("gaff2");

		companion object {

			fun from(ffname: String) =
				values().find { it.id == ffname }

			fun fromOrThrow(ffname: String) =
				from(ffname) ?: throw IllegalArgumentException("forcefield $ffname is not supported by Parmchk")
		}
	}

	fun run(serviceDir: Path, mol2: String, atomTypes: AtomTypes): Results {

		val homePath = serviceDir.resolve("progs/ambertools").toAbsolutePath()
		val parmchkPath = homePath.resolve("bin/parmchk2").toAbsolutePath()

		tempFolder("parmchk") { cwd ->

			// write the input files
			mol2.write(cwd.resolve("mol.mol2"))

			// start parmchk
			val process = ProcessBuilder()
				.command(
					parmchkPath.toString(),
					"-i", "mol.mol2",
					"-f", "mol2",
					"-s", atomTypes.id,
					"-o", "frcmod"
				)
				.apply {
					environment().apply {
						put("AMBERHOME", homePath.toString())
					}
				}
				.directory(cwd.toFile())
				.stream()
				.waitFor()

			// return the results
			return Results(
				process.exitCode,
				process.console.toList(),
				try {
					cwd.resolve("frcmod").read()
				} catch (ex: IOException) {
					null
				}
			)
		}
	}

	class Exception(val msg: String, val mol2: String, val results: Results) : RuntimeException(StringBuilder().apply {
		append(msg)
		append("\n\n")
		append("Input:\n")
		append(mol2)
		append("\n\n")
		append("console:\n")
		results.console.forEach {
			append(it)
			append("\n")
		}
	}.toString())
}
