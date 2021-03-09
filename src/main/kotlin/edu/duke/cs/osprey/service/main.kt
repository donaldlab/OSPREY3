package edu.duke.cs.osprey.service

import java.nio.file.Paths


fun main() {

	val cwd = Paths.get(System.getProperty("user.dir"))

	// start the server
	println("Starting ${OspreyService.name} v ${OspreyService.version} ...")
	OspreyService.Instance(cwd, wait = true)
}
