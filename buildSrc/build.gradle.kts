
plugins {
	`kotlin-dsl`
	// NOTE: putting plugins here only affects this script file
	// it has no effect on the buildSrc/src/** files
	// to do that, you'll need to add a dependency instead
	// which means somehow translating the plugin ids to maven coordinates
	// this might help, good luck:
	// https://plugins.gradle.org/
}

repositories {
	mavenCentral()
	maven("https://plugins.gradle.org/m2/")
}

dependencies {

	val kotlinVersion = "1.5.31"

	// kotlin("jvm")
	implementation("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")

	// kotlin("plugin.serialization")
	implementation("org.jetbrains.kotlin:kotlin-serialization:$kotlinVersion")

	// id("org.beryx.runtime") version "1.12.5"
	implementation("org.beryx:badass-runtime-plugin:1.12.5")

	// SSH client, BSD license: https://github.com/mwiede/jsch
	implementation("com.github.mwiede:jsch:0.1.66")
}
