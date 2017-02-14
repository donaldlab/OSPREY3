
import os, glob
import jpype


c = None
_classpath = []


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


def start(heapSizeMB=1024, enableAssertions=False):

	# build JVM launch args
	args = [
		jpype.getDefaultJVMPath(),
		'-xmx%dM' % heapSizeMB,
		'-Djava.class.path=%s' % makeClasspath()
	]
	if enableAssertions:
		args.append("-ea")

	# start the JVM
	jpype.startJVM(*args)
	
	# set up class factories
	global c
	c = Packages()
	c.java = jpype.JPackage('java')
	c.javax = jpype.JPackage('javax')


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

