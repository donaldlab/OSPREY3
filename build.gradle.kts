
import java.io.BufferedWriter
import java.io.OutputStreamWriter
import java.net.URL
import java.nio.charset.StandardCharsets
import java.nio.file.attribute.PosixFilePermission
import java.nio.file.StandardCopyOption
import java.nio.file.Files
import java.nio.file.Path

plugins {
	application
	idea
}

val projectDir = project.projectDir.toPath().toAbsolutePath()
val pythonSrcDir = projectDir.resolve("python")
val pythonBuildDir = projectDir.resolve("build/python")
val pythonWheelDir = pythonBuildDir.resolve("wheel")
val pythonWheelhouseDir = pythonSrcDir.resolve("wheelhouse")
val docSrcDir = pythonSrcDir.resolve("doc")
val docBuildDir = pythonBuildDir.resolve("doc")


group = "edu.duke.cs"
version = Files.readAllLines(projectDir.resolve("resources/config/version"))[0]


repositories {
	jcenter()
}

java {
	sourceCompatibility = JavaVersion.VERSION_1_8
	targetCompatibility = JavaVersion.VERSION_1_8
	sourceSets {
		get("main").apply {
			java.srcDir("src")
			resources.srcDir("resources")
		}
		get("test").apply {
			java.srcDir("test")
			resources.srcDir("test-resources")
		}
	}
}

idea {
	module {
		// use the same output folders as gradle, so the pythonDevelop task works correctly
		outputDir = java.sourceSets["main"].output.classesDirs.singleFile
		testOutputDir = java.sourceSets["test"].output.classesDirs.singleFile
		inheritOutputDirs = false
	}
}

application {
	mainClassName = "edu.duke.cs.osprey.control.Main"
}

dependencies {

	// test dependencies
	testCompile("org.hamcrest:hamcrest-all:1.3")
	testCompile("junit:junit:4.12")

	// compile dependencies
	compile("colt:colt:1.2.0")
	compile("org.apache.commons:commons-math3:3.6.1")
	compile("org.apache.commons:commons-collections4:4.1")
	compile("commons-io:commons-io:2.5")
	compile("com.joptimizer:joptimizer:3.5.1")
	compile("org.ojalgo:ojalgo:41.0.0")
	compile("org.jogamp.gluegen:gluegen-rt:2.3.2")
	compile("org.jogamp.jocl:jocl:2.3.2")
	compile("org.mapdb:mapdb:3.0.5")
	compile("org.apache.xmlgraphics:batik-svggen:1.9.1")
	compile("org.apache.xmlgraphics:batik-svg-dom:1.9.1")

	// for JCuda, gradle tries (and fails) download the natives jars automatically,
	// so turn off transitive dependencies. we'll deal with natives manually
	compile("org.jcuda:jcuda:0.8.0") {
		isTransitive = false
	}

	// build systems never seem to like downloading from arbitrary URLs for some reason...
	// so make a helper func to do it for us
	fun url(urlString: String): ConfigurableFileCollection {

		// parse the URL
		val url = URL(urlString)
		val filename = urlString.split("/").last()
		val libDir = projectDir.resolve("lib")
		val filePath = libDir.resolve(filename)

		// create a gradle task to download the file
		val downloadTask by project.tasks.creating {
			Files.createDirectories(libDir)
			url.openStream().use {
				Files.copy(it, filePath, StandardCopyOption.REPLACE_EXISTING)
			}
		}

		// return a files collection that depends on the task
		return files(filePath) {
			builtBy(downloadTask)
		}
	}

	// TPIE-Java isn't in the maven/jcenter repos yet, download directly from Github
	compile(url("https://github.com/donaldlab/TPIE-Java/releases/download/v1.1/edu.duke.cs.tpie-1.1.jar"))

	// native libs for GPU stuff
	listOf("natives-linux-amd64", "natives-macosx-universal", "natives-windows-amd64").forEach {
		runtime("org.jogamp.gluegen:gluegen-rt:2.3.2:$it")
		runtime("org.jogamp.jocl:jocl:2.3.2:$it")
	}
	listOf("linux-x86_64", "apple-x86_64", "windows-x86_64").forEach {
		runtime("org.jcuda:jcuda-natives:0.8.0:$it")
	}
}

distributions {

	get("main").apply {
		baseName = "osprey-cli"
		contents {
			into("") { // project root
				from("README.rst")
				from("LICENSE.txt")
				from("CONTRIBUTING.rst")
			}
			listOf("1CC8", "1FSV", "2O9S_L2", "2RL0.kstar", "3K75.3LQC", "4HEM", "4NPD").forEach {
				into("examples/$it") {
					from("examples/$it")
				}
			}
		}
	}

	create("python").apply {
		baseName = "osprey-python"
		contents {
			into("") { // project root
				from("README.rst")
				from("LICENSE.txt")
				from("CONTRIBUTING.rst")
				from(pythonBuildDir) {
					include("install.sh")
					include("install.bat")
					include("uninstall.sh")
					include("uninstall.bat")
				}
			}
			into("doc") {
				from(docBuildDir) {
					exclude(".doctrees")
				}
			}
			listOf("python.GMEC", "python.KStar", "gpu").forEach {
				into("examples/$it") {
					from("examples/$it")
				}
			}
			into("wheelhouse") {
				from(pythonWheelhouseDir)
				from(pythonBuildDir) {
					include("*.whl")
				}
			}
		}
	}
}

val pythonCmd = "python2"
val pipCmd = "pip2"

tasks {

	// get vals for tasks we care about that were defined by plugins
	val jar = tasks.getByName("jar")

	// turn off tar distributions
	"distTar" {
		enabled = false
	}
	"pythonDistTar" {
		enabled = false
	}

	val compileCuda_residueForcefield by creating(Exec::class) {
		nvcc(this, "residueForcefield")
	}

	val compileCuda_residueCcd by creating(Exec::class) {
		nvcc(this, "residueCcd", maxRegisters=64)
	}

	val compileCuda by creating {
		description = "Compile cuda kernels"
		dependsOn(
			compileCuda_residueForcefield,
			compileCuda_residueCcd
		)
	}

	val cleanDoc by creating(Delete::class) {
		group = "documentation"
		description = "Cleans python documentation"
		delete(docBuildDir)
	}
	val makeDoc by creating(Exec::class) {
		group = "documentation"
		description = "Build python documentation"
		workingDir = docSrcDir.toFile()
		commandLine("sphinx-build", "-b", "html", "-j", "2", ".", "$docBuildDir")
	}
	val remakeDoc by creating {
		group = "documentation"
		description = "runs cleanDoc, then makeDoc to refresh all changes to documentation"
		dependsOn(cleanDoc, makeDoc)
	}

	val pythonDevelop by creating(Exec::class) {
		group = "develop"
		description = "Install python package in development mode"
		workingDir = pythonSrcDir.toFile()
		commandLine(pipCmd, "install",
			"--user", "--editable",
			".", // path to package to install, ie osprey
			"--no-index", "--use-wheel", "--find-links=$pythonWheelhouseDir" // only use wheelhouse to resolve dependencies
		)
		doLast {
			Files.createDirectories(pythonBuildDir)
			val classpathPath = pythonBuildDir.resolve("classpath.txt")
			Files.write(classpathPath, java.sourceSets["main"].runtimeClasspath.files.map { it.toString() })
		}
	}

	val pythonUndevelop by creating(Exec::class) {
		group = "develop"
		description = "Uninstall development mode python package"
		workingDir = pythonSrcDir.toFile()
		commandLine(pipCmd, "uninstall", "--yes", "osprey")
	}

	val pythonWheel by creating(Exec::class) {
		group = "build"
		description = "Build python wheel"
		inputs.file(tasks.getByName("jar"))
		outputs.dir(pythonWheelDir)
		doFirst {

			// delete old cruft
			delete {
				delete(pythonWheelDir) // TODO: this apparently does not delete the folder at all??!
				delete(fileTree(pythonBuildDir) {
					include("*.whl")
				})
			}

			// copy python sources
			copy {
				from(pythonSrcDir) {
					includeEmptyDirs = false
					include("osprey/*.py")
				}
				from(".") {
					include("*.rst")
				}
				into(pythonWheelDir.toFile())
			}

			// copy setup.py, but change the rootDir
			copy {
				from(pythonSrcDir) {
					include("setup.py")
				}
				into(pythonWheelDir.toFile())
				filter { line ->
					if (line.startsWith("rootDir = ")) {
						"rootDir = \"$projectDir\""
					} else {
						line
					}
				}
			}

			val libDir = pythonWheelDir.resolve("osprey/lib")

			// copy osprey jar
			copy {
				from(jar)
				into(libDir.toFile())
			}

			// copy java libs
			copy {
				from(java.sourceSets["main"].runtimeClasspath.files
					.filter { it.extension == "jar" }
				)
				into(libDir.toFile())
			}
		}
		workingDir = pythonWheelDir.toFile()
		commandLine(pythonCmd, "setup.py", "bdist_wheel", "--dist-dir", pythonBuildDir.toString())
	}

	val pythonInstallScripts by creating {
		group = "build"
		description = "Make install scripts for python distribution"
		doLast {
			writeScripts(
				"install",
				"""
				|$pipCmd uninstall -y osprey JPype-py2
				|$pipCmd install --user numpy
				|$pipCmd install --user osprey --no-index --use-wheel --find-link=wheelhouse
				""".trimMargin()
			)
		}
	}

	val pythonUninstallScripts by creating {
		group = "build"
		description = "Make uninstall scripts for python distribution"
		doLast {
			writeScripts(
				"uninstall",
				"$pipCmd uninstall -y osprey JPype-py2"
			)
		}
	}

	// insert some build steps before we build the python dist
	"pythonDistZip" {
		dependsOn(pythonWheel, makeDoc, pythonInstallScripts, pythonUninstallScripts)
	}
}

fun nvcc(exec: Exec, kernelName: String, maxRegisters: Int? = null, profile: Boolean = false) {

	val args = mutableListOf("nvcc")

	if (profile) {
		// if profiling, compile for one arch with profiling/debug info
		// NOTE: change this to your GPU's arch
		args.addAll(listOf("-cubin", "-gencode=arch=compute_61,code=sm_61", "-lineinfo", "--ptxas-options=-v"))
	} else {
		// otherwise, compile for all archs
		// see Maxwell compatibility guide:
		// http://docs.nvidia.com/cuda/maxwell-compatibility-guide/index.html#building-maxwell-compatible-apps-using-cuda-6-0
		args.addAll(listOf("-fatbin",
			"-gencode=arch=compute_20,code=sm_20",
			"-gencode=arch=compute_30,code=sm_30",
			"-gencode=arch=compute_35,code=sm_35",
			"-gencode=arch=compute_50,code=sm_50",
			"-gencode=arch=compute_52,code=sm_52",
			"-gencode=arch=compute_60,code=sm_60",
			"-gencode=arch=compute_61,code=sm_61",
			"-gencode=arch=compute_62,code=sm_62",
			"-gencode=arch=compute_62,code=compute_62"
		))
	}

	if (maxRegisters != null) {
		args.addAll(listOf("-maxrregcount", "$maxRegisters"))
	}

	args.addAll(listOf("$kernelName.cu", "-o", "$kernelName.bin"))

	exec.workingDir = file("resources/gpuKernels/cuda")
	exec.commandLine(args)
}

fun List<String>.writeToFile(path: Path, newline: String) {
	BufferedWriter(OutputStreamWriter(Files.newOutputStream(path), StandardCharsets.UTF_8.newEncoder())).use { out ->
		for (line in this) {
			out.append(line)
			out.append(newline)
		}
	}
}

fun writeShellScript(filename: String, cmd: String) {

	val file = pythonBuildDir.resolve("$filename.sh")

	"""
		|#! /bin/sh
		|$cmd
	""".trimMargin()
		.split("\n")
		.writeToFile(file, "\n")

	// set the shell script executable
	Files.setPosixFilePermissions(file, Files.getPosixFilePermissions(file).apply {
		add(PosixFilePermission.OWNER_EXECUTE)
	})
}

fun writeBatchScript(filename: String, cmd: String) {

	val file = pythonBuildDir.resolve("$filename.bat")

	"""
		|@echo off
		|$cmd
	""".trimMargin()
		.split("\n")
		.writeToFile(file, "\r\n")
}

fun writeScripts(filename: String, cmd: String) {
	writeShellScript(filename, cmd)
	writeBatchScript(filename, cmd)
}
