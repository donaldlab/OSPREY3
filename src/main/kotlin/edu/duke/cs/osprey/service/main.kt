package edu.duke.cs.osprey.service

import java.nio.file.Paths


fun main(args: Array<String>) {

	// read command-line arguments
	val port = args
		.indexOfFirst { it == "--port" }
		.takeIf { it >= 0 }
		?.let { i ->
			// read the following argument
			val str = args.getOrNull(i + 1)
				?: throw IllegalArgumentException("missing port")
			str.toIntOrNull()
				?: throw IllegalArgumentException("port not recognized: $str")
		}
		?: OspreyService.defaultPort

	// get the current working directory
	val cwd = Paths.get(System.getProperty("user.dir"))

	// start the server
	println("Starting ${OspreyService.name} v ${OspreyService.version} on port $port ...")
	OspreyService.Instance(
		cwd,
		wait = true,
		port = port
	)
}
