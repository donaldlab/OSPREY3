/*
** This file is part of OSPREY 3.0
** 
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
** 
** OSPREY relies on grants for its development, and since visibility
** in the scientific literature is essential for our success, we
** ask that users of OSPREY cite our papers. See the CITING_OSPREY
** document in this distribution for more information.
** 
** Contact Info:
**    Bruce Donald
**    Duke University
**    Department of Computer Science
**    Levine Science Research Center (LSRC)
**    Durham
**    NC 27708-0129
**    USA
**    e-mail: www.cs.duke.edu/brd/
** 
** <signature of Bruce Donald>, Mar 1, 2018
** Bruce Donald, Professor of Computer Science
*/

import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

import osprey.*
import osprey.build.*


plugins {
	`java-library`
	idea
	kotlin("jvm") // no version here, already specified in buildSrc
	kotlin("plugin.serialization") // no version here, already specified in buildSrc
	id("org.openjfx.javafxplugin") version("0.0.7")
	id("org.beryx.runtime") // no version here, already specified in buildSrc
}

javafx {
	version = "11"
	modules("javafx.controls")
}


group = "edu.duke.cs"

/**
 * Version number for Osprey itself
 *
 * This version number is largely cosmetic, compared to the versioning scheme for the Osprey Service
 * But it does help users provide some information to developers when reporting issues.
 */
version = "3.2"

repositories {
	mavenCentral()
}


java {
	toolchain {
		languageVersion.set(JavaLanguageVersion.of(Jvm.javaLangVersion))
	}
}

idea {
	module {
		// use the same output folders as gradle, so the pythonDevelop task works correctly
		outputDir = sourceSets.main.output.classesDirs.files.first()
		testOutputDir = sourceSets.test.output.classesDirs.files.first()
		inheritOutputDirs = false
	}
}

dependencies {

	// kotlin runtime
	implementation(kotlin("stdlib-jdk8"))
	implementation(kotlin("reflect"))

	// test dependencies
	testImplementation("org.hamcrest:hamcrest-all:1.3")
	testImplementation("junit:junit:4.12")
	testImplementation("io.kotlintest:kotlintest-runner-junit5:3.4.0")

	// handle logging
	implementation("ch.qos.logback:logback-classic:1.2.3")
	implementation("org.slf4j:jul-to-slf4j:1.7.30")

	// internal osprey libs
	implementation("colt:colt:1.2.0")
	implementation("org.apache.commons:commons-math3:3.6.1")
	implementation("org.apache.commons:commons-collections4:4.1")
	implementation("com.joptimizer:joptimizer:3.5.1")
	implementation("org.ojalgo:ojalgo:41.0.0")
	implementation("org.jogamp.gluegen:gluegen-rt:2.3.2")
	implementation("org.jogamp.jocl:jocl:2.3.2")
	implementation("org.mapdb:mapdb:3.0.5")
	implementation("org.apache.xmlgraphics:batik-svggen:1.9.1")
	implementation("org.apache.xmlgraphics:batik-svg-dom:1.9.1")
	implementation("com.github.haifengl:smile-core:1.5.1")
	implementation("com.github.haifengl:smile-netlib:1.5.1")
	implementation("de.lmu.ifi.dbs.elki:elki:0.7.1")
	implementation("ch.obermuhlner:big-math:2.3.0")
	implementation("org.joml:joml:1.9.19")
	implementation("org.tukaani:xz:1.8")
	implementation("com.hazelcast:hazelcast:4.0")
	implementation("net.java.dev.jna:jna:5.5.0")
	implementation("com.google.guava:guava:29.0-jre")
	implementation("org.apache.commons:commons-lang3:3.4")
	implementation("commons-io:commons-io:2.5")
	implementation("org.tomlj:tomlj:1.0.0")
	implementation(fileTree("lib"))

	val autoValueVersion = 1.7
	implementation("ch.obermuhlner:big-math:2.3.0")
	implementation("com.fasterxml.jackson.core:jackson-databind:2.9.9")
	implementation("com.fasterxml.jackson.dataformat:jackson-dataformat-yaml:2.9.9")
	implementation("com.beust:jcommander:1.72")
	implementation("one.util:streamex:0.7.3")
	implementation(platform("software.amazon.awssdk:bom:2.15.48"))
	implementation("software.amazon.awssdk:s3")
	implementation("org.postgresql:postgresql:42.2.16")
	implementation("org.sql2o:sql2o:1.6.0")
	implementation("com.google.auto.value:auto-value-annotations:$autoValueVersion")
	annotationProcessor("com.google.auto.value:auto-value:$autoValueVersion")
	testImplementation("org.junit.jupiter:junit-jupiter:5.4.2")
	testImplementation("org.assertj:assertj-core:3.18.1")


	val ktorVersion = "1.5.4"

	// used by the gui
	implementation("com.cuchazinteractive:kludge:0.2")
	implementation("io.ktor:ktor-client-apache:$ktorVersion")
	implementation("io.ktor:ktor-client-serialization-jvm:$ktorVersion")

	// used by the service
	implementation("io.ktor:ktor-server-netty:$ktorVersion")
	implementation("io.ktor:ktor-serialization:$ktorVersion")

	// for JCuda, gradle tries (and fails) download the natives jars automatically,
	// so turn off transitive dependencies. we'll deal with natives manually
	implementation("org.jcuda:jcuda:0.8.0") {
		isTransitive = false
	}
	implementation("one.util:streamex:0.7.3")

	// native libs for GPU stuff
	listOf("natives-linux-amd64", "natives-macosx-universal", "natives-windows-amd64").forEach {
		runtimeOnly("org.jogamp.gluegen:gluegen-rt:2.3.2:$it")
		runtimeOnly("org.jogamp.jocl:jocl:2.3.2:$it")
	}
	listOf("linux-x86_64", "apple-x86_64", "windows-x86_64").forEach {
		runtimeOnly("org.jcuda:jcuda-natives:0.8.0:$it")
	}

	// used by the build system
	testImplementation("org.json:json:20210307")
	val dokkaVersion = "1.6.0"
	testImplementation("org.jetbrains.dokka:dokka-cli:$dokkaVersion")
	testImplementation("org.jetbrains.dokka:dokka-base:$dokkaVersion")
	testImplementation("org.jetbrains.dokka:dokka-core:$dokkaVersion")
}


/* NOTE: the IDE thinks `jvmArgs` and `args` are not nullable and shows warnings
	(and the Kotlin language rules agree with that, as far as I can tell),
	but for some reason, the Kotlin compiler thinks they are nullable
	so we need the not null assertions (ie !!). Maybe a compiler bug?
*/
@Suppress("UNNECESSARY_NOT_NULL_ASSERTION")
tasks.withType<JavaExec> {
	Jvm.addModuleArgs(jvmArgs)
	doFirst {

		// make Gradle dump the java command to the console when running
		// since IntelliJ doesn't do it when using Gradle as the runner
		val jvmArgsStr = jvmArgs!!.joinToString(" ")
		val classpathStr = classpath.joinToClasspath()
		val argsStr = args!!.joinToString(" ")
		println("""
			|Java:
			|	$workingDir
			|	$executable
			|	$jvmArgsStr -cp "$classpathStr" $main $argsStr
		""".trimMargin())
	}
}


tasks.withType<JavaCompile> {
	options.compilerArgs.addAll(Jvm.moduleArgs)
}

tasks.withType<KotlinCompile> {

	kotlinOptions {

		jvmTarget = "1.8"

		// enable experimental features so we can use the fancy ktor stuff
		freeCompilerArgs += "-Xuse-experimental=kotlin.Experimental"
		freeCompilerArgs += "-XXLanguage:+InlineClasses"
	}
}

tasks.withType<Test> {
	// the default 512m is too little memory to run test designs
	maxHeapSize = "2g"
	useJUnit()
    failFast = true
	Jvm.addModuleArgs(jvmArgs)

	testLogging {
		setExceptionFormat("full")
        events("passed", "skipped", "failed")
	}
}


// assume we're doing a dev build if the top task is "classes"
if (gradle.startParameter.taskNames.any { it.endsWith(":classes") }) {
	System.setProperty("isDev", true.toString())
}
val isDev = object {
	override fun toString() = System.getProperty("isDev") ?: false.toString()
}


tasks {

	processResources {

		// always update the build properties, versions, etc
		outputs.upToDateWhen { false }

		var version = version.toString()

		// append the CI build ID, if available
		version += if (rootProject.hasProperty("AZURE_BUILD_ID")) {
			val versionId = rootProject.property("AZURE_BUILD_ID")
			".$versionId"
		} else {
			"-dev"
		}

		// write the build properties so the app can find them
		from(sourceSets.main.resources.srcDirs) {
			include("${Jvm.packagePath}/build.properties")
			expand(
				"dev" to isDev,
				"version" to version,
				"versionService" to Builds.service.version
			)
		}

		// I have no idea why Gradle thinks the build properties are duplicated
		// (there's only one instance of build.properties anywhere in the project),
		// but make sure we include the expanded version of the file rather than ignoring it
		duplicatesStrategy = DuplicatesStrategy.INCLUDE
	}
}


// add a whole bunch of other tasks to do interesting things
makeCudaTasks()
makePythonTasks()
makeLicenseTasks()
makeDocsTasks()
makeBuildTasks()
makeBuildServerTasks()
makeBuildDesktopTasks()
makeBuildServiceTasks()
makeBuildServiceDockerTasks()
