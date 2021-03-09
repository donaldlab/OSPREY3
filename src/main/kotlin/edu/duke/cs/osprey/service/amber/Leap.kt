package edu.duke.cs.osprey.service.amber

import edu.duke.cs.osprey.service.*
import java.io.IOException
import java.nio.file.Path


object Leap {

	/** remove all characters that aren't in the allowed charaters list, to prevent any injection attacks on the LEAP scripts */
	fun sanitizeToken(ffname: String) =
		ffname.replace(Regex("[^\\w.']"), "")

	data class Results(
		val exitCode: Int,
		val console: List<String>,
		val files: Map<String,String?>
	) {

		/* TODO: scan the console outfor for errors/warnings/etc
		     apparently LEaP *does* use consistent tagging in the output.
		     See: Amber manual 13.5, Error Handling and Reporting
		*/
	}

	fun run(
		serviceDir: Path,
		filesToWrite: Map<String,String>,
		commands: String,
		filesToRead: List<String> = emptyList(),
		debugFiles: List<String> = emptyList()
	): Results {

		val leapPath = serviceDir.resolve("progs/ambertools/bin/teLeap").toAbsolutePath()
		val datPath = serviceDir.resolve("progs/ambertools/dat/leap").toAbsolutePath()

		tempFolder("leap") { cwd ->

			// write the files
			for ((filename, content) in filesToWrite) {
				content.write(cwd.resolve(filename))
			}

			// write the leap commands
			val commandsPath = cwd.resolve("commands")
			commands.write(commandsPath)

			// start leap
			val process = ProcessBuilder()
				.command(
					if (debugFiles.isEmpty()) {
						leapPath.toString()
					} else {
						"$leapPath.debug"
					},
					"-I", datPath.resolve("prep").toString(),
					"-I", datPath.resolve("prep/oldff").toString(),
					"-I", datPath.resolve("lib").toString(),
					"-I", datPath.resolve("lib/oldff").toString(),
					"-I", datPath.resolve("parm").toString(),
					"-I", datPath.resolve("cmd").toString(),
					"-I", datPath.resolve("cmd/oldff").toString(),
					"-f", commandsPath.toAbsolutePath().toString()
				)
				.apply {
					if (debugFiles.isNotEmpty()) {
						environment()["MESSAGEON"] = debugFiles.joinToString(" ")
					}
				}
				.directory(cwd.toFile())
				.stream()
				.waitFor()

			// return the results
			return Results(
				process.exitCode,
				process.console.toList(),
				filesToRead.associateWith {
					try {
						cwd.resolve(it).read()
					} catch (ex: IOException) {
						null
					}
				}
			)
		}
	}

	class Exception(val msg: String, val inMol: String, val results: Results) : RuntimeException(StringBuilder().apply {
		append(msg)
		append("\n\n")
		append("input PDB or Mol2:\n")
		append(inMol)
		append("\n\n")
		append("console:\n")
		results.console.forEach {
			append(it)
			append("\n")
		}
	}.toString())
}
