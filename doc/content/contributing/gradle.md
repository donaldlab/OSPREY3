+++
title = "Gradle Build System"
weight = 1
+++


OSPREY uses [Gradle][gradle] to manage builds and automate other development tasks too.

[gradle]: https://gradle.org/


## Build definition files

 * `/build.gradle.kts`:
   This is the main definition file for the build system.
   This file is a [Kotlin][kotlin] script using the [Gradle Kotlin DSL][gradle-kotlin-dsl]
   that configures the build system to compile OSPREY, among other development tasks.
 * `/settings.gradle.kts`:
   This is a configuration file for the Gradle build.
 * `/buildSrc`:
   This is where most of the build system definitions are for OSPREY.
   The [buildSrc][build-src] folder is like an extension of the `/build.gradle.kts` file,
   but can be organized into packages, multiple files, and be tested just like regular Kotlin code.
 * `/buildSrc/build.gradle.kts`:
   This is the main definition file for the `buildSrc` part of the build definitions.
   It mostly just declares the dependencies for the code inside `buildSrc`.
 * `buildSrc/src/main/kotlin/**/*.kt`:
   These files are the main collection of build definitions for OSPREY.
   All of the custom build tasks are defined here. If you find yourself needing to change
   OSPREY's build behavior, or add new build behavior, this is probably where you'll want to look first.
 * `buildSrc/src/main/java/**/*.java`:
   Here lies the implemention of OSPREY's Javadoc tools, to help build the documentation for the Java code.
 * `buildSrc/src/main/python/**/*.py`:
   Here lies the code for OSPREY's tools for generating documentation for the Python code.
 * `buildSrc/src/test/**/*`:
   These are automated tests for the build code.
 
[kotlin]: https://kotlinlang.org/
[gradle-kotlin-dsl]: https://docs.gradle.org/current/userguide/kotlin_dsl.html
[build-src]: https://docs.gradle.org/current/userguide/organizing_gradle_projects.html#sec:build_sources


## Running Gradle Tasks

Your IDE should give you tools to run Gradle Tasks easily.

For example, in IntelliJ IDEA, there is a Gradle
panel in the upper right of your window. Clicking on that shows all the tasks available. Double-clicking
one of those will run the selected task.

If for some reason, your IDE can't run the Gradle task, you can also use the command line.
In Linux and OSX, open a terminal in the OSPREY project folder and run:
```shell
./gradlew $TASK
```
where `$TASK` is the name of the Gradle task you want to run.

In windows, the command looks slightly different:
```shell
gradlew.bat $TASK
```


## Stopping long-running Gradle tasks {#stop-task}

For Gradle tasks that keep running until explicitly stopped (like servers), how you stop the task depends on how
you started it.

For Gradle tasks started in IntelliJ IDEA, look for the square red Stop button in the Run window.

For Gradle tasks started in a console, simply press Ctrl-C to ask the process to exit.
