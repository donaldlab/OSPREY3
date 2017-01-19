
import sys
import jvm, jpype


# NOTE: this var gets set by the build system during packaging
# so the release version of this script will point to the final jar file for osprey
# instead of the development classes folder
_ospreyPaths = ['../../build/output/*.jar', '../../bin']
 
c = None

def _javaAwareExcepthook(exctype, value, traceback):

	# show original python exception info
	sys.__excepthook__(exctype, value, traceback)

	# if this is a java exception, show java info too
	if value.message is not None and value.stacktrace is not None:
		if value.message() is not None:
			print(value.message())
		print(value.stacktrace())


def start(heapSizeMB=1024, enableAssertions=False):

	# setup a global exception handler to show java exception info
	sys.excepthook = _javaAwareExcepthook
	
	# start the jvm
	for path in _ospreyPaths:
		jvm.addClasspath(path)
	jvm.start(heapSizeMB, enableAssertions)

	# set up class factories
	global c
	c = jvm.Packages()
	c = jpype.JPackage('edu.duke.cs.osprey')

	# print the preamble
	print("OSPREY %s" % c.control.Main.Version)


def loadPdb(path):
	return c.structure.PDBIO.readFile(path)

