package edu.duke.cs.osprey.service.amber

import edu.duke.cs.osprey.service.*
import java.nio.file.Path


object Sander {

	// we customized sander to give us bigger limits here
	const val maxCommandLineSize = 65664
	const val maxRestraintMaskSize = 65536

	fun sanitizeRestraintMask(mask: String) =
		mask.replace(Regex("[^\\w',@]"), "")

	data class Results(
		val exitCode: Int,
		val console: List<String>,
		val coords: List<Point3d>?
	)

	fun run(
		serviceDir: Path,
		top: String,
		crd: String,
		commands: String
	): Results {

		val sanderPath = serviceDir.resolve("progs/ambertools/bin/sander").toAbsolutePath()

		tempFolder("sander") { cwd ->

			// inputs:
			// mdin input control data for the min/md run
			// mdout output user readable state info and diagnostics
			//    or "-o stdout" will send output to stdout (to the terminal) instead of to a file.
			// mdinfo output latest mdout-format energy info
			//    (is this really required?)
			// prmtop input molecular topology, force field, periodic box type, atom and residue names
			// inpcrd input initial coordinates and (optionally) velocities and periodic box size

			val topname = "mol.top"
			val crdname = "mol.crd"
			val cmdname = "mol.cmd"
			val outname = "mol.restrt"

			// write the input files
			top.write(cwd.resolve(topname))
			crd.write(cwd.resolve(crdname))

			// check the command line sizes
			commands
				.split("\n")
				.filter { it.length > maxCommandLineSize }
				.takeIf { it.isNotEmpty() }
				?.let { lines ->
					throw IllegalArgumentException(
						"Sander commands lines size are over the limit of $maxCommandLineSize:\n" +
						lines.joinToString("\n")
					)
				}

			// sander crashes if the commands don't end with a newline, so just add one to be safe
			"$commands\n".write(cwd.resolve(cmdname))

			// start sander
			val process = ProcessBuilder()
				.command(
					sanderPath.toString(),
					"-i", cmdname,
					"-o", "stdout",
					"-p", topname,
					"-c", crdname,
					"-r", outname,
					"-ref", crdname
				)
				.directory(cwd.toFile())
				.stream()
				.waitFor()

			// TODO: track progress info

			// parse the output coords
			val coords = cwd.resolve(outname)
				.takeIf { it.exists }
				?.read()
				?.let { text ->
					ArrayList<Point3d>().apply {
						var lines = text.lines()

						// skip the first two lines
						lines = lines.subList(2, lines.size)

						for (line in lines) {

							// lines look like e.g.:
							//  14.7619439  27.0623578  24.0946254  13.9092238  25.8758637  24.2319541
							//   3.7678268  22.1883445   9.1170323
							// It turns out lines can also look like this:
							// 207.4266747-244.7747434 207.2027093 -49.6475693  79.7673583-122.2447495
							// So we also want to split using a lookahead for the symbol '-'

							fun String.toDoubleOrThrow() =
								toDoubleOrNull()
								?: throw IllegalArgumentException("$this doesn't appear to be a number\nin line\n$line")

							val parts = line
								.split(Regex("(\\s|(?=-))"))
								.filter { it.isNotBlank() }

							if (parts.size >= 3) {
								add(Point3d(
									parts[0].toDoubleOrThrow(),
									parts[1].toDoubleOrThrow(),
									parts[2].toDoubleOrThrow()
								))
							}
							if (parts.size >= 6) {
								add(Point3d(
									parts[3].toDoubleOrThrow(),
									parts[4].toDoubleOrThrow(),
									parts[5].toDoubleOrThrow()
								))
							}
						}
					}
				}

			return Results(
				process.exitCode,
				process.console.toList(),
				coords
			)
		}
	}

	class Exception(val msg: String, val console: Collection<String>) : RuntimeException(StringBuilder().apply {
		append(msg)
		append("\n\n")
		append("console:\n")
		console.forEach {
			append(it)
			append("\n")
		}
	}.toString())
}
