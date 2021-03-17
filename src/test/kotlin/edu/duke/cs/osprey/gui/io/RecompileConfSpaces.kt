package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.tools.LZMA2
import edu.duke.cs.osprey.gui.compiler.ConfSpaceCompiler
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.prep.ConfSpace
import java.nio.file.Files
import java.nio.file.Paths


/**
 * Just a simple tool to easily recompile all the conf spaces
 * in Osprey's test suite and examples to the new version of the compiled conf space format.
 */
fun main() = withService {

	val extension = ".confspace"
	val ospreyDir = Paths.get(".")
	val dirs = listOf(
		"src/test/resources/confSpaces",
		"examples/python.ccs/F98Y"
	).map { ospreyDir.resolve(it) }

	val netChargesByMolName = mapOf(
		"NDP" to -3,
		"06W" to -1,
		"EPE" to 0
	)

	for (dir in dirs) {
		Files.list(dir)
			.filter { it.fileName.toString().endsWith(extension) }
			.forEach { inPath ->

				val filename = inPath.fileName.toString()
				val basename = filename.substring(0, filename.length - extension.length)

				// load the conf space
				println("loading $basename ...")
				val confSpace = ConfSpace.fromToml(inPath.read())

				// compile it
				ConfSpaceCompiler(confSpace).run {

					// use default setings, to match classic osprey
					forcefields.add(Forcefield.Amber96)
					forcefields.add(Forcefield.EEF1)

					// add necessary net charges
					for ((type, mol) in confSpace.mols) {
						netCharges.get(mol, type)?.netCharge = netChargesByMolName.getValue(mol.name)
					}

					println("compiling $basename ...")
					val report = compile().run {
						printUntilFinish(5000)
						report!!
					}

					// if there was an error, throw it
					report.error?.let { throw Error("can't compile $basename", it) }

					// otherwise, yay it worked!
					val compiledConfSpace = report.compiled!!

					// save the compressed conf space
					val outPath = inPath.parent.resolve("$basename.ccsx")
					println("saving $basename ...")
					LZMA2.compress(compiledConfSpace.toBytes()).write(outPath)
				}
			}
	}
}
