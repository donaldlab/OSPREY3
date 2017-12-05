
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


def start(heapSizeMiB=1024, enableAssertions=False, stackSizeMiB=None, garbageSizeMiB=None):

	# build JVM launch args
	args = [
		jpype.getDefaultJVMPath(),
		'-xmx%dM' % heapSizeMiB,
		'-Djava.class.path=%s' % makeClasspath()
	]
	if enableAssertions:
		args.append("-ea")
	if stackSizeMiB is not None:
		args.append('-Xss%sM' % stackSizeMiB)
	if garbageSizeMiB is not None:
		args.append('-XX:MaxNewSize=%dM' % garbageSizeMiB)
	if _nativesDir is not None:
		args.append('-Djava.library.path=%s' % _nativesDir)

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
	except TypeError:
		raise ValueError('%s is not a recognized Java class' % jclass)

	# get the inner class
	return getJavaClass('%s$%s' % (classname, inner_class_name))


def boxDouble(val):
	if val is None:
		return None
	return c.java.lang.Double(val)
