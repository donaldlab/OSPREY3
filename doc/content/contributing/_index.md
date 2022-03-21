+++
title = "Contributing"
weight = 4
+++

# Contributing to Osprey

Welcome to OSPREY development! OSPREY is open-source software and contributions are welcome.
See the instructions below to get started with the code and compiling it.


## Set up development environment

OSPREY develoment depends on a few general software packages:

* [git](https://git-scm.com/)
* [Java JDK](https://adoptium.net/) (v17 or newer)
* [Python](https://www.python.org/) (v3 or newer)
* [pip](https://pip.pypa.io/en/stable/)

Many developer machines have these tools installed already.
If you don't have them yet, install them before proceeding to the next sections.


### Clone the git repostory

The main OSPREY code repository is hosted at [Github](https://github.com):

https://github.com/donaldlab/OSPREY3

To develop OSPREY, first clone the git repository to a folder of your choosing:
```shell
mkdir osprey
cd osprey
git clone git@github.com:donaldlab/OSPREY3.git .
```


### Overview

OSPREY is primarily written in the [Java][java] programming language, with a small [Python][python]
API for users. Some newer parts of the code are written in the [Kotlin] programming language as well.

[java]: https://en.wikipedia.org/wiki/Java_(software_platform)
[python]: https://www.python.org/
[kotlin]: https://kotlinlang.org/

OSPREY uses the [Gradle][gradle] build system with the [Kotlin DSL][gradle-kotlin-dsl]
to automate certain development tasks like:

* Download Java libraries (e.g. jar files)
* Maintain the Java classpath needed to compile or run OSPREY
* Compile code and build the OSPREY executables
* Compile Osprey documentation
* Build the Osprey distribution files
* Configure Java IDEs
* Configure Python development environment

[gradle]: https://gradle.org/
[gradle-kotlin-dsl]: https://blog.gradle.org/kotlin-meets-gradle

You don't need to install Gradle, since it's already included in the OSPREY repo.


### Compile the code

To begin developing OSPREY, just use your favorite Java IDE to import the Gradle project.
Gradle has built-in or plugin support in many popular IDEs, including [IntelliJ IDEA][idea],
[Eclipse][eclipse], and [Netbeans][netbeans]. We usually use IntelliJ IDEA in the Donald Lab,
so you're likely to have the best experience with that one. IDEA also has the best support for the
[Kotlin][kotlin] language, which is used in the Gradle build scripts and the newer Desktop and Service code.

The import process in your IDE should walk you through the steps needed
to setup OSPREY for compilation and testing.

[eclipse]: https://www.eclipse.org/
[idea]: https://www.jetbrains.com/idea/
[netbeans]: https://netbeans.org/

Once your IDE seems happy with the OSPREY project, tell it to run the main method in the
`edu.duke.cs.osprey.CompileIt` class. This class is in the `src/test/java` folder.
This entry point to the code does nothing on its own,
so it's a nice and simple test of the compile-and-run process of your IDE.

Running the main method there starts the Gradle build process and then
lanches the Java program. When it's done, you should see a simple message on the console:
```
OSPREY compiled and ran!
```

Congratulations! You can now write code, compile, and run OSPREY!


### Set up the Python development environment for OSPREY Server

OSPREY Server's main user interface is Python scripting. To enable developers to run Python scripts
that use the code from the local git clone of OSPREY, we'll need to install OSPREY into the Python
environment. To do that, run the Gradle task `pythonDevelop` (in the `develop` group).

{{% notice tip %}}
To see how to run tasks in Gradle, go to [Gradle Tasks]({{< ref "gradle-tasks" >}})
{{% /notice %}}

Once installed into the Python environment, calls to e.g., `import osprey` in Python scripts will run
the Python module in `src/main/python/osprey` directly from the source code.

However, the Java and Kotlin codes can't be run directly from source like the Python code can.
When developing your Python scripts, you'll still need to compile the Java and Kotlin code before
your changes will be available to the Python scripts. The `CompileIt.main()` entry point is a useful
tool for compiling the rest of the code for the Python scripts. Run `CompileIt.main()` after editing
the Java/Kotlin code but before running your Python scripts to make sure the Python scripts can see
your latest changes.

To uninstall OSPREY from your Python environment, run the `pythonUndevelop` Gradle task (also in
the `develop` group).

This will remove the OSPREY package from the Python environment. Removing the development verison
of OSPREY is necessary if you want to install any release version of OSPREY Server.
