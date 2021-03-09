package edu.duke.cs.osprey.service.clashes

import edu.duke.cs.osprey.service.Point3d
import edu.duke.cs.osprey.service.stream
import edu.duke.cs.osprey.service.tempFolder
import edu.duke.cs.osprey.service.write
import java.nio.file.Path
import java.util.regex.Pattern


object Probe {

	// yeah, probe isn't technically part of AmberTools,
	// but it's convenient to use the multi-platform support we created for AmberTools

	data class Results(
		val exitCode: Int,
		val console: List<String>,
		val groups: Map<String,Group>
	)

	data class Group(
		val id: String,
		val dots: Map<String,List<Point3d>>,
		val vectors: Map<String,List<Pair<Point3d,Point3d>>>
	)

	fun run(
		serviceDir: Path,
		pdb: String
	): Results {

		val probePath = serviceDir.resolve("progs/probe/probe.2.16r72.130520").toAbsolutePath()

		tempFolder("probe") { cwd ->

			// write the input files
			pdb.write(cwd.resolve("mol.pdb"))

			// start antechamber
			val process = ProcessBuilder()
				.command(
					probePath.toString(),
					"mol.pdb"
					// TOOD: specialized args?
				)
				.directory(cwd.toFile())
				.stream()
				.waitFor()

			class TempGroup {
				val dots = HashMap<String,ArrayList<Point3d>>()
				val vectors = HashMap<String,ArrayList<Pair<Point3d,Point3d>>>()
			}
			val groups = HashMap<String,TempGroup>()
			var group = null as TempGroup?

			// e.g.
			// @dotlist {x} color=white master={vdw contact}
			// @vectorlist {x} color=red master={small overlap}
			val masterPattern = Pattern.compile(".*\\s+master=\\{?([\\w\\s]+)}?")

			fun getOrMakeGroup(line: String): TempGroup? =
				masterPattern.matcher(line)
					.takeIf { it.matches() && it.groupCount() >= 1 }
					?.let {
						val id = it.group(1)
						groups.getOrPut(id) { TempGroup() }
					}

			val prefixPattern = "\\{[^}]+}(\\w+)\\s+.+"
			val numberPattern = "[0-9\\.\\-\\+]+"
			val pointPattern = "($numberPattern),($numberPattern),($numberPattern)"

			// e.g.
			// 	{ C   ALA   2  A}blue  'M' 12.753,24.762,22.059
			//	{"}blue  'M' 12.501,24.791,22.060
			val dotPattern = Pattern.compile("$prefixPattern\\s+$pointPattern")

			// e.g.
			// 	{ CG  GLU   3  A}yellowtint P  'P' 10.611,25.725,19.578 {"}yellowtint   'P' 10.612,25.727,19.577
			//	{"}yellowtint P  'P' 10.614,25.868,19.794 {"}yellowtint   'P' 10.616,25.873,19.790
			//val vectorPattern = Pattern.compile("$prefixPattern\\s+$pointPattern\\s+.*\\s+$pointPattern")
			val vectorPattern = Pattern.compile("$prefixPattern\\s+$pointPattern\\s+.+\\s+$pointPattern")

			var isDot = false

			// read the dots and vectors from the console
			for (line in process.console) {

				if (line.startsWith("@dotlist")) {

					group = getOrMakeGroup(line)
					isDot = true

				} else if (line.startsWith("@vectorlist")) {

					group = getOrMakeGroup(line)
					isDot = false

				} else if (group != null && line.startsWith("{")) {

					if (isDot) {

						dotPattern.matcher(line)
							.takeIf { it.matches() && it.groupCount() >= 4 }
							?.let {
								val color = it.group(1)
								group.dots
									.getOrPut(color) { ArrayList() }
									.add(Point3d(
										it.group(2).toDouble(),
										it.group(3).toDouble(),
										it.group(4).toDouble()
									))
							}

					} else {

						vectorPattern.matcher(line)
							.takeIf { it.matches() && it.groupCount() >= 7 }
							?.let {
								val color = it.group(1)
								group.vectors
									.getOrPut(color) { ArrayList() }
									.add(Point3d(
										it.group(2).toDouble(),
										it.group(3).toDouble(),
										it.group(4).toDouble()
									) to Point3d(
										it.group(5).toDouble(),
										it.group(6).toDouble(),
										it.group(7).toDouble()
									))
							}
					}
				}
			}

			return Results(
				process.exitCode,
				process.console.toList(),
				groups.mapValues { (name, group) ->
					Group(name, group.dots, group.vectors)
				}
			)
		}
	}
}
