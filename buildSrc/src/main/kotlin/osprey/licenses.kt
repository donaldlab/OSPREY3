package osprey

import org.gradle.api.Project
import org.gradle.kotlin.dsl.creating
import org.gradle.kotlin.dsl.getValue
import java.nio.file.Path


fun Project.makeLicenseTasks() {

	// TODO: replace this with the nifty licenser plugin?
	//   https://github.com/Minecrell/licenser
	@Suppress("UNUSED_VARIABLE")
	val updateLicenseHeaders by tasks.creating {
		group = "build"
		description = "updates license headers in all source files"
		doLast {
			updateLicenseHeaders()
		}
	}
}


enum class HeaderResult {
	Updated,
	Ignored
}


// TODO: replace this custom license header code with the licenser plugin
// https://github.com/Minecrell/licenser

fun Project.updateLicenseHeaders() {

	val header = """
		|This file is part of OSPREY 3.0
		|
		|OSPREY Protein Redesign Software Version 3.0
		|Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
		|
		|OSPREY is free software: you can redistribute it and/or modify
		|it under the terms of the GNU General Public License version 2
		|as published by the Free Software Foundation.
		|
		|You should have received a copy of the GNU General Public License
		|along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
		|
		|OSPREY relies on grants for its development, and since visibility
		|in the scientific literature is essential for our success, we
		|ask that users of OSPREY cite our papers. See the CITING_OSPREY
		|document in this distribution for more information.
		|
		|Contact Info:
		|   Bruce Donald
		|   Duke University
		|   Department of Computer Science
		|   Levine Science Research Center (LSRC)
		|   Durham
		|   NC 27708-0129
		|   USA
		|   e-mail: www.cs.duke.edu/brd/
		|
		|<signature of Bruce Donald>, Mar 1, 2018
		|Bruce Donald, Professor of Computer Science
	""".trimMargin().lines()

	val deleteTheseAutoHeaders = listOf(
		""" |/*
			| * To change this template, choose Tools | Templates
			| * and open the template in the editor.
			| */
		""".trimMargin(),
		""" |/*
			| * To change this license header, choose License Headers in Project Properties.
			| * To change this template file, choose Tools | Templates
			| * and open the template in the editor.
			| */
		""".trimMargin()
	)

	fun applyCHeader(lines: MutableList<String>): HeaderResult {

		// extract the existing license header, if any
		var readMode = 0
		var existingHeader = lines
			.takeWhile {
				val line = it.trim()
				when (readMode) {
					0 -> {
						if (line.startsWith("/*")) {
							readMode = 1
							return@takeWhile true
						}
					}
					1 -> {
						if (line.startsWith("**")) {
							return@takeWhile true
						} else if (line.startsWith("*/")) {
							readMode = 2
							return@takeWhile true
						}
					}
				}
				return@takeWhile false
			}
			.map { it.substring(2).trim() }
		if (existingHeader.size >= 3) {
			for (i in 0 until existingHeader.size) {
				lines.removeAt(0)
			}
			existingHeader = existingHeader.subList(1, existingHeader.size - 1)
		}

		// if it matches the desired header, then we're done
		if (existingHeader == header) {
			return HeaderResult.Ignored
		}

		// trim blank lines
		while (lines.firstOrNull()?.isBlank() == true) {
			lines.removeAt(0)
		}

		// add the new header
		lines.add(0, "")
		lines.add(0, "*/")
		for (i in 0 until header.size) {
			lines.add(0, "** " + header[header.size - i - 1])
		}
		lines.add(0, "/*")

		return HeaderResult.Updated
	}

	fun applyPythonHeader(lines: MutableList<String>): HeaderResult {

		// extract the existing license header, if any
		val existingHeader = lines
			.takeWhile { it.startsWith("##") }
			.map { it.substring(2).trim() }
		for (i in 0 until existingHeader.size) {
			lines.removeAt(0)
		}

		// if it matches the desired header, then we're done
		if (existingHeader == header) {
			return HeaderResult.Ignored
		}

		// trim blank lines
		while (lines.firstOrNull()?.isBlank() == true) {
			lines.removeAt(0)
		}

		// add the new header
		lines.add(0, "")
		for (i in 0 until header.size) {
			lines.add(0, "## " + header[header.size - i - 1])
		}

		return HeaderResult.Updated
	}

	fun applyHeader(path: Path, applier: (MutableList<String>) -> HeaderResult) {

		var text = path.read()

		// remove any headers automatically added by IDEs or other tools
		var removedAutoHeaders = false
		for (autoHeader in deleteTheseAutoHeaders) {
			val newtext = text.replace(autoHeader, "")
			if (newtext != text) {
				removedAutoHeaders = true
				text = newtext
			}
		}

		val lines = text.lines().toMutableList()

		// trim blank lines from the top
		while (lines.firstOrNull()?.isBlank() == true) {
			lines.removeAt(0)
		}

		// keep one blank line on the bottom
		// NOTE: a trailing newline creates a blank line at the end of the file,
		// so it's sufficient to remove all blank entries in the lines list
		while (lines.lastOrNull()?.isBlank() == true) {
			lines.removeAt(lines.size - 1)
		}

		// anything left?
		if (lines.isEmpty()) {
			return
		}

		if (removedAutoHeaders || applier(lines) == HeaderResult.Updated) {
			path.writeLines(lines)
			println("updated: $path")
		}
	}

	fun applyHeaders(dirname: String, filter: (String) -> Boolean, applier: (MutableList<String>) -> HeaderResult) {

		// for each matched file in the folder (and subfolders)
		val dir = projectPath / dirname
		dir.walk { stream ->
			stream
				.filter { filter(it.fileName.toString()) }
				.forEach { applyHeader(it, applier) }
		}
	}

	// apply header to java files
	for (dirname in listOf("src", "test")) {
		applyHeaders(
			dirname,
			filter = { it.endsWith(".java") },
			applier = ::applyCHeader
		)
	}

	// apply header to this file
	applyHeader(projectPath / "build.gradle.kts", ::applyCHeader)

	// apply header to kernel files
	for (dirname in listOf("src/main/resources/gpuKernels")) {
		applyHeaders(
			dirname,
			filter = { it.endsWith(".cu") || it.endsWith(".cl") },
			applier = ::applyCHeader
		)
	}

	// apply header to python files
	for (dirname in listOf("python")) {
		applyHeaders(
			dirname,
			filter = { it.endsWith(".py") },
			applier = ::applyPythonHeader
		)
	}

	// NOTE: don't apply the header to the python example scripts.
	// there's no need to scare osprey users with legalese in the tutorials
}
