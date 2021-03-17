package edu.duke.cs.osprey.gui.io

import edu.duke.cs.osprey.gui.compiler.ConfSpaceCompiler
import edu.duke.cs.osprey.gui.forcefield.Forcefield
import edu.duke.cs.osprey.gui.prep.ConfSpace
import java.nio.file.Paths


fun main() {

	val toml = Paths.get("").read()
	val confSpace = ConfSpace.fromToml(toml)

	// compile it
	withService {

		ConfSpaceCompiler(confSpace).run {

			// use default setings
			forcefields.add(Forcefield.Amber96)
			forcefields.add(Forcefield.EEF1)

			// add necessary net charges
			netCharges[confSpace.findMol("ANP")]?.netCharge = -1

			println("compiling ...")
			val report = compile().run {
				printUntilFinish(5000)
				report!!
			}

			// if there was an error, throw it
			report.error?.let { throw Error("can't compile", it) }
		}
	}
}
