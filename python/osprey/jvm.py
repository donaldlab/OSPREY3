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

import os, glob
import jpype


c = None
_classpath = []
_nativesDir = None


class Packages(object):
	pass


def addClasspath(path):

	global _classpath

	for path in glob.glob(path):
		_classpath.append(path)


def makeClasspath():

	global _classpath

	# on windows, the classpath separator is ; instead of :
	if os.name == 'nt':
		separator = ';'
	else:
		separator = ':'
	return separator.join(_classpath)


def setNativesDir(path):
	global _nativesDir
	_nativesDir = path


def start_with_args(jrePath, jvmArgs):
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

	start_with_args(jrePath, args)

	if attachJvmDebugger:
		input("Attach the JVM debugger now, set your breakpoints, and hit [enter] to continue:")


def shutdown():
	# NOTE: jpype says this function doesn't even work for most JVMs
	# I guess we'll just include it anyway, for completeness' sake
	jpype.shutdownJVM()


def attachThread():
	jpype.attachThreadToJVM()


def toArrayList(items):
	jlist = c.java.util.ArrayList()
	for item in items:
		jlist.add(item)
	return jlist


def toFile(path):
	return c.java.io.File(path)


def getJavaClass(classname):
	jclass = c.java.lang.Class
	classloader = c.java.lang.ClassLoader.getSystemClassLoader()
	return jclass.forName(classname, True, classloader)
	

def getInnerClass(jclass, inner_class_name):
	'''
	Gets the inner class from the Java outer class.
	
	:param jclass: The Java class, provided by a JPype JPackage object.
	:param str inner_class_name: The simple name of the inner class.
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
	if val is None:
		return None
	return c.java.lang.Double(val)

def boxInt(val):
	if val is None:
		return None
	return c.java.lang.Integer(val)

def boxLong(val):
	if val is None:
		return None
	return c.java.lang.Long(val)

def makeIntArray(val):
    array = jpype.JArray(c.java.lang.Integer)(len(val) or val)
    for i in range(0, len(val)):
        array[i] = (c.java.lang.Integer(val[i]))
    return array

def makeStringArray(val):
    array = jpype.JArray(c.java.lang.String)(val)
    return array

