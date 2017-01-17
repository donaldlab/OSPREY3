
# OSPREY

Open-Source Protein REdesign for You!

OSPREY is developed and maintained by the [Donald Lab](http://www.cs.duke.edu/donaldlab/home.php)
in the [Department of Computer Science](http://www.cs.duke.edu/)
at [Duke University](https://www.duke.edu/)


## Installation

TODO


## Run Osprey

### using Python scripts

TODO: explain python interface, should be preferred way to using Osprey


### using the command-line interface

Run Osprey at the command line by entering the following command into a shell
```
java -jar /path/to/osprey.jar [commands]
```

For `commands`, try `help` to show a list of commands.

To show the version of your Osprey installation, try:
```
java -jar osprey.jar version
```

To try a GMEC-based design, try:
```
java -jar osprey.jar FindGMEC /path/to/config1 /path/to/config2 ...
```

## Contributing

### Setup Eclipse build path:

Osprey uses the [Jerkar](http://project.jerkar.org/) build system.
* [installed Jerkar if needed](http://project.jerkar.org/documentation/latest/getting_started.html)
* run this in a shell at the project root:
```
jerkar eclipse#generateFiles
```

This will download the needed Java libraries generate the `.classpath` file for Eclipse.
Then Eclipse can build the project and you'll be able to launch Osprey from within Eclipse for debugging.


### Build Osprey Jar

`jerkar doPack`
