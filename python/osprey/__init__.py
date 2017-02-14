
import sys, jpype
import jvm, augmentation


# NOTE: this var gets set by the build system during packaging
# so the release version of this script will point to the final jar file for osprey
# instead of the development classes folder
_ospreyPaths = ['../../build/output/*.jar', '../../bin']

c = None

WILD_TYPE = None

def _javaAwareExcepthook(exctype, value, traceback):

	# show original python exception info
	sys.__excepthook__(exctype, value, traceback)

	# if this is a java exception, show java info too
	if value.message is not None and 'stacktrace' in value and value.stacktrace is not None:
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
	c = jpype.JPackage('edu.duke.cs.osprey')

	augmentation.init()

	# init other globals
	global WILD_TYPE
	WILD_TYPE = c.confspace.Strand.WildType

	# print the preamble
	print("OSPREY %s" % c.control.Main.Version)


def loadPdb(path):
	return c.structure.PDBIO.readFile(path)


def Strand(path):
	mol = loadPdb(path)
	# TODO: expose builder args
	return c.confspace.Strand.builder(mol).build()


def ConfSpace(strands, shellDist=None):

	builder = c.confspace.SimpleConfSpace.builder()

	# get a list of strands, even if we were passed just one strand
	try:
		# try iterating over strands
		strandsList = [strand for strand in strands]
	except TypeError:
		# otherwise, it's just one strand
		strandsList = [strands]

	# now add all the strands to the conf space
	for strandInfo in strandsList:

		try:
			# check for args
			strand = strandInfo[0]
			flex = strandInfo[1:]
		except TypeError:
			# nope, just a strand
			strand = strandInfo
			flex = []

		builder.addStrand(strand, jvm.toArrayList(flex))

	# add the shell distance if needed
	if shellDist is not None:
		builder.setShellDistance(shellDist)
	
	return builder.build()


def StrandFlex():
	# TODO: implement me
	pass


def EnergyMatrix(confSpace, cacheFile=None):
	
	# TODO: expose builder options
	builder = c.ematrix.SimplerEnergyMatrixCalculator.builder(confSpace).build()

	if cacheFile is not None:
		return builder.calcEnergyMatrix(jvm.toFile(cacheFile))
	else:
		return builder.calcEnergyMatrix()


def GMECFinder(confSpace, emat, confLog=None, printIntermediateConfs=None):
	
	builder = c.gmec.SimpleGMECFinder.builder(confSpace, emat)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not None:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	return builder.build()

