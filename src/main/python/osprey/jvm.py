## This file is part of OSPREY 3.0
## 
## OSPREY Protein Redesign Software Version 3.0
## Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
## 
## OSPREY is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License version 2
## as published by the Free Software Foundation.
## 
## You should have received a copy of the GNU General Public License
## along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
## 
## OSPREY relies on grants for its development, and since visibility
## in the scientific literature is essential for our success, we
## ask that users of OSPREY cite our papers. See the CITING_OSPREY
## document in this distribution for more information.
## 
## Contact Info:
##    Bruce Donald
##    Duke University
##    Department of Computer Science
##    Levine Science Research Center (LSRC)
##    Durham
##    NC 27708-0129
##    USA
##    e-mail: www.cs.duke.edu/brd/
## 
## <signature of Bruce Donald>, Mar 1, 2018
## Bruce Donald, Professor of Computer Science


'''
Functions to manange the Java Virtual Machine (JVM)
'''


import os, glob
import jpype


c = None
'''
Class factory for JVM runtime classes.

Osprey exposes raw Java classes to Python code using a bridge library called
[JPype](https://github.com/jpype-project/jpype). To access the raw Java classes from Python code,
call properties and functions on this factory object.

For example, to construct an instance of a Java class called `ArrayList` in the package `java.util` call its constructor:
```python
osprey.jvm.c.java.util.ArrayList()
```

Currently, the `java.*` and `javax.*` packages are supported.
'''

_classpath = []
_nativesDir = None


class Packages(object):
	pass


def addClasspath(path):
	'''
	Adds the given path to the class path.

	The next call to #start will use this classpath to launch the JVM.

	# Arguments
	path `str`: the classpath entry, usually a path to a `.jar` file, of a folder of `.class` files.
	'''

	global _classpath

	for path in glob.glob(path):
		_classpath.append(path)


def makeClasspath():
	'''
	Returns the classpath that has been collected by #addClasspath.

	# Returns
	`str` where all classpath entries have been joined by the system classpath separator character.
	'''

	global _classpath

	# on windows, the classpath separator is ; instead of :
	if os.name == 'nt':
		separator = ';'
	else:
		separator = ':'
	return separator.join(_classpath)


def setNativesDir(path):
	'''
	Sets the directory for native libraries.

	The next call to #start will use this directory to launch the JVM.
	'''

	global _nativesDir
	_nativesDir = path


def start_with_args(jrePath, jvmArgs):
	'''
	Start the JVM with the given arguments.

	Settings accumulated by #addClasspath and #setNativesDir are ignored.

	# Arguments
	jrePath `str`: Path to the JVM installation folder
	jvmArgs `[str]`: Command-line arguments to pass to the JVM
	'''

	# if no path to a JRE was given, assume Java is installed somewhere,
	# and try to determine the path automatically
	if jrePath is None:
		jrePath = jpype.getDefaultJVMPath()

	# build JVM launch args
	args = [
		jrePath,
		'-Djava.class.path=%s' % makeClasspath()
	]

	args += jvmArgs

	# start the JVM
	try:
		jpype.startJVM(*args, convertStrings=True)
	except TypeError:
		# JPype-py2 doesn't support the convertStrings kwarg
		jpype.startJVM(*args)

	# set up class factories
	global c
	c = Packages()
	c.java = jpype.JPackage('java')
	c.javax = jpype.JPackage('javax')


def start(jrePath, heapSizeMiB=1024, enableAssertions=False, stackSizeMiB=None, garbageSizeMiB=None, allowRemoteManagement=False, attachJvmDebugger=False):
	'''
	Starts the JVM using the settings accumulated by #addClasspath and #setNativesDir,
	as well as settings from the keyword arguments.

	# Arguments
	jrePath `str`: Path to the JVM installation folder
	heapSizeMiB `int`: The size of the JVM heap, in [Mebibytes](https://simple.wikipedia.org/wiki/Mebibyte), ie the `-Xmx` flag.
	enableAssertions `bool`: True to enable Java assertions, ie the `-ea` flag.
	stackSizeMiB `int`: The size of the JVM stack, in Mebibytes, ie the `-Xss` flag.
	{{% notice note %}}
	You can leave this at the default value, unless the JVM crashes with an error saying the stack size is too small.
	Then try running your program again with a larger stack size.
	{{% /notice %}}
	garbageSizeMiB `int`: The size of the JVM heap reserved for garbage collection, in Mebibytes, ie the `-XX:MaxNewSize` flag.
	                      By default, the JVM will reserve roughly 1/3 of the heap for garbage collection, but this may
	                      be inefficient for Osprey runs that need a lot of non-garbage memory.
	                      To reserve more non-garbage memory for Osprey calculations and less memory for garbage collection,
	                      supply a value here to limit the garbage reservation to a fixed size.
	                      Often, a gibibyte or two is plenty of memory for garbage collection, but the total heap size
	                      may be much larger than a few gibibytes.
	allowRemoteManagement `bool`: `True` to turn on [JMX remote management](https://www.oracle.com/technical-resources/articles/javase/jmx.html).
	attachJvmDebugger `bool`: `True` to turn on the JVM debugger service, see [Attaching a Debugger](https://jpype.readthedocs.io/en/latest/userguide.html#attaching-a-debugger).
	'''

	args = ['-Xmx%dM' % heapSizeMiB]

	if enableAssertions:
		args.append("-ea")
	if stackSizeMiB is not None:
		args.append('-Xss%sM' % stackSizeMiB)
	if garbageSizeMiB is not None:
		args.append('-XX:MaxNewSize=%dM' % garbageSizeMiB)
	if _nativesDir is not None:
		args.append('-Djava.library.path=%s' % _nativesDir)
	if allowRemoteManagement:
		args.append('-Dcom.sun.management.jmxremote.ssl=false')
		args.append('-Dcom.sun.management.jmxremote.authenticate=false')
		args.append('-Dcom.sun.management.jmxremote.port=9010')
		args.append('-Dcom.sun.management.jmxremote.rmi.port=9011')
		args.append('-Djava.rmi.server.hostname=localhost')
		args.append('-Dcom.sun.management.jmxremote.local.only=false')
	if attachJvmDebugger:
		# See https://jpype.readthedocs.io/en/latest/userguide.html#attaching-a-debugger for how to attach
		# a java debugger to the running process
		args.append("-Xdebug")
		args.append("-Xnoagent")
		args.append("-Xrunjdwp:transport=dt_socket,server=y,address=12999,suspend=n")

	# suppress warnings about illegal reflection operations in newer JVMs
	# yes, the code that generates these warnings can handle the exceptions that future JVMs
    # will eventually throw when the reflection operations get forcibly denied
	# there's no need to scare our users with warning messages
	args.append('--add-opens=java.base/java.text=ALL-UNNAMED')

	# enable FFI libraries in newer JVMs
	args.append('--add-modules=jdk.incubator.foreign')

	start_with_args(jrePath, args)

	if attachJvmDebugger:
		input("Attach the JVM debugger now, set your breakpoints, and hit [enter] to continue:")


def shutdown():
	'''
	Shuts down the JVM and releases all its resources from the process.
	The Python VM will keep running as usual though.
	'''

	jpype.shutdownJVM()


def toArrayList(items):
	'''
	Convert the python list into a java `ArrayList` object.

	# Arguments
	items `[]`: the python list

	# Returns
	${type_jre(java.util.ArrayList)}
	'''

	jlist = c.java.util.ArrayList()
	for item in items:
		jlist.add(item)
	return jlist


def toFile(path):
	'''
	Converts a path into a java `File` object.

	# Arguments
	path `str`: A path

	# Returns
	${type_jre(java.io.File)}
	'''
	return c.java.io.File(path)


def getJavaClass(classname):
	'''
	Returns the class object for the given class name.

	# Arguments
	classname `str`: The fully-qualified name of a Java class.
	                 Inner classes should be separated from outer classes with a `$` character,
	                 rather than a `.` character.

	# Returns
	${type_jre(java.lang.Class)}
	'''
	jclass = c.java.lang.Class
	classloader = c.java.lang.ClassLoader.getSystemClassLoader()
	return jclass.forName(classname, True, classloader)
	

def getInnerClass(jclass, inner_class_name):
	'''
	Gets the inner class from the Java outer class.

	# Arguments
	jclass: The JPype object representing the Java class, provided by a class factory like #c.
	inner_class_name `str`: The simple name of the inner class.

	# Returns
	${type_jre(java.lang.Class)}
	'''

	# get the class name, if this is even a class
	try:
		classname = jclass.__javaclass__.getName()
		return getJavaClass('%s$%s' % (classname, inner_class_name))
	except AttributeError:
		# must be a new version of jpype, try the newer API
		classname = jclass.class_.getName()
		return jpype.JClass(getJavaClass('%s$%s' % (classname, inner_class_name)))
	except TypeError:
		raise ValueError('%s is not a recognized Java class' % jclass)


def boxDouble(val):
	'''
	Converts the python floating-point number into a Java `Double` object.

	# Arguments
	val `float`: the number

	# Returns
	${type_jre(java.lang.Double)}
	'''

	if val is None:
		return None
	return c.java.lang.Double(val)


def boxInt(val):
	'''
	Converts the python integer number into a Java `Integer` object.

	# Arguments
	val `int`: the number

	# Returns
	${type_jre(java.lang.Integer)}
	'''
	if val is None:
		return None
	return c.java.lang.Integer(val)


def boxLong(val):
	'''
	Converts the python integer number into a Java `Long` object.

	# Arguments
	val `int`: the number

	# Returns
	${type_jre(java.lang.Long)}
	'''
	if val is None:
		return None
	return c.java.lang.Long(val)
