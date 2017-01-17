
# Setup Eclipse build path:

* make sure Jerkar is installed (http://project.jerkar.org/)
* run this in a shell at the project root:
```
jerkar eclipse#generateFiles
```

This will download the needed Java libraries generate the `.classpath` file for Eclipse.


# Build Osprey Jar

`jerkar doPack`


# Run Osprey

Run Osprey by entering the following command into a shell
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

# TODO: write more of this