
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

	val kotlinVersion = "1.6.10"

	// kotlin("jvm")
	implementation("org.jetbrains.kotlin:kotlin-gradle-plugin:$kotlinVersion")

	// kotlin("plugin.serialization")
	implementation("org.jetbrains.kotlin:kotlin-serialization:$kotlinVersion")

	// id("org.jetbrains.dokka")
	implementation("org.jetbrains.dokka:dokka-gradle-plugin:$kotlinVersion")

	// id("org.beryx.runtime") version "1.12.5"
	// https://badass-runtime-plugin.beryx.org
	implementation("org.beryx:badass-runtime-plugin:1.12.5")

	// SSH client, BSD license: https://github.com/mwiede/jsch
	implementation("com.github.mwiede:jsch:0.1.66")

	// JSON library, The JSON License: https://json.org/license.html
	implementation("org.json:json:20210307")

	// used by the kotlin dokka plugin
	implementation("org.jetbrains.dokka:dokka-base:$kotlinVersion")
	implementation("org.jetbrains.dokka:dokka-core:$kotlinVersion")

	// test dependencies
	testImplementation("org.hamcrest:hamcrest-all:1.3")
	testImplementation("junit:junit:4.12")
}

tasks.withType<Test>().configureEach {
	useJUnit()
}
