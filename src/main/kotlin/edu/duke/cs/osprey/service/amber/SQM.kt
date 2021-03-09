package edu.duke.cs.osprey.service.amber


object SQM {

	enum class ErrorType {
		NoConvergence,
		BadNetCharge
	}

	/**
	 * See AMBER 19 manual pg 148
	 */
	data class Options(

		/** level of theory to use for the QM region of the simulation (Hamiltonian) */
		var qmTheory: String = "AM1",

		/** gradient RMS tolerence */
		var grmsTol: Double = 0.0005,

		/** self-consistent field convergence */
		var scfConv: Double = 1e-10,

		/** number of iterations for direct inversion of the iterative subspace */
		var ndiisAttempts: Int = 700,

		/** maximum number of minimization cycles to allow, using the xmin minimizer */
		var maxCyc: Int = 9999,

		/** print the progress of the minimization every ntpr steps */
		var ntpr: Int = 10
	) {
		/**
			Converts the double value to fortran format, eg:
				1e-10 -> 1.d-10
				1.2435 -> 1.2345d0
		*/
		private fun Double.toFortran() =
			"%e".format(this).let { str ->
				if (str.contains('.')) {
					str.replace("e", "d")
				} else {
					str.replace("e", ".d")
				}
			}

		override fun toString() =
			mutableListOf(
				"qm_theory" to "\"$qmTheory\"",
				"grms_tol" to grmsTol.toFortran(),
				"scfconv" to scfConv.toFortran(),
				"ndiis_attempts" to "$ndiisAttempts",
				"maxcyc" to "$maxCyc",
				"ntpr" to "$ntpr"
			)
			.joinToString(", ") { (key, value) ->
				"$key=$value"
			}
			.let { "$it," } // need a trailing comma, since this list gets prepended to another list by antechamber
	}

	data class Results(
		val exitCode: Int?,
		var console: List<String>?,
		val inFile: String,
		val outFile: String?
	) {

		/**
		 * Sometimes, other programs call SQM for us, so we don't
		 * get the exitCode and console. In that case, just handle
		 * the in/out files.
		 */
		constructor (inFile: String, outFile: String) : this(
			exitCode = null,
			console = null,
			inFile = inFile,
			outFile = outFile
		)

		/**
		 * Any errors we could parse from the logs.
		 */
		val errors: List<Pair<ErrorType,String>> = ArrayList<Pair<ErrorType,String>>().apply {

			if (outFile == null) {
				return@apply
			}

			val lines = outFile.lines().iterator()
			while (lines.hasNext()) {

				val line = lines.next().trim()

				/* eg:
					QMMM: ERROR!
					QMMM: Unable to achieve self consistency to the tolerances specified
					QMMM: No convergence in SCF after   1000 steps.
					QMMM: E =  -0.1882E+06 DeltaE =  -0.1728E+02 DeltaP =   0.1252E+00
					QMMM: Smallest DeltaE =   0.1630E-08 DeltaP =   0.3595E-05 Step =      3
				*/
				if (line == "QMMM: ERROR!") {

					add(ErrorType.NoConvergence to """
						|${lines.next()}
						|${lines.next()}
					""".trimMargin())

				/* or eg:
					QMMM: System specified with odd number of electrons (   53)
					QMMM: but odd spin (  1). You most likely have the charge of
					QMMM: QM region (qmcharge) set incorrectly. Correct error and re-run calculation.
				*/
				} else if (line.startsWith("QMMM: System specified with odd number of electrons")) {

					add(ErrorType.BadNetCharge to """
						|$line
						|${lines.next()}
						|${lines.next()}
					""".trimMargin())
				}
			}
		}
	}
}