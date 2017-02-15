
import sys, os, jpype
import jvm, augmentation

getJavaClass = augmentation.getJavaClass


# NOTE: this var gets set by the build system during packaging
# so the release version of this script will point to the final jar file for osprey
# instead of the development classes folder
_ospreyPaths = ['../../build/output/*.jar', '../../bin']

c = None

WILD_TYPE = None
Forcefield = None
SolvationForcefield = None
LovellRotamers = 0 # arbitrary value, doesn't matter

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
	global Forcefield
	Forcefield = getJavaClass('energy.forcefield.ForcefieldParams$Forcefield')
	global SolvationForcefield
	SolvationForcefield = getJavaClass('energy.forcefield.ForcefieldParams$SolvationForcefield')

	# print the preamble
	print("OSPREY %s" % c.control.Main.Version)


def readTextFile(path):
	return c.tools.FileTools.readFile(path)


def _ensureText(textOrPath):
	if os.path.exists(textOrPath):
		# it's a path, read the file
		return readTextFile(textOrPath)
	else:
		# it's text, just return it
		return textOrPath


def Parallelism(cpuCores=None, gpus=None, streamsPerGpu=None):
	if gpus is not None:
		return c.parallelism.Parallelism.makeGpu(gpus, streamsPerGpu)
	elif cpuCores is not None:
		return c.parallelism.Parallelism.makeCpu(cpuCores)
	else:
		return c.parallelism.Parallelism.makeDefault()


def TemplateLibrary(forcefield=None, templateCoords=None, rotamers=None, backboneDependentRotamers=None):
	builder = c.restypes.GenericResidueTemplateLibrary.builder()

	if forcefield is not None:
		builder.setForcefield(forcefield)

	if templateCoords is not None:
		builder.setTemplateCoords(_ensureText(templateCoords))

	if rotamers is not None:
		if rotamers == LovellRotamers:
			builder.setLovellRotamers()
		else:
			builder.setRotamers(_ensureText(rotamers))

	if backboneDependentRotamers is not None:
		builder.setBackboneDependentRotamers(_ensureText(backboneDependentRotamers))

	return builder.build()
	

def readPdb(path):
	mol = c.structure.PDBIO.readFile(path)
	print('read PDB file from file: %s' % path)
	return mol


def writePdb(mol, path):

	# unbox the mol if we need to
	if isinstance(mol, c.confspace.ParametricMolecule):
		mol = mol.mol

	c.structure.PDBIO.writeFile(mol, path)

	print('write PDB file to file: %s' % path)


def Strand(pathOrMol, residues=None):

	# did we get a path or a molecule?
	if isinstance(pathOrMol, c.structure.Molecule):
		mol = pathOrMol
	else:
		mol = readPdb(pathOrMol)

	builder = c.confspace.Strand.builder(mol)

	if residues is not None:
		builder.setResidues(residues[0], residues[1])
	
	return builder.build()


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


def ForcefieldParams():
	return c.energy.forcefield.ForcefieldParams()

def EnergyMatrix(confSpace, ffparams=None, parallelism=None, cacheFile=None):
	
	# TODO: expose builder options
	builder = c.ematrix.SimplerEnergyMatrixCalculator.builder(confSpace)

	if ffparams is not None:
		builder.setForcefieldParams(ffparams)

	if parallelism is not None:
		builder.setParallelism(parallelism)

	ematcalc = builder.build()

	if cacheFile is not None:
		return ematcalc.calcEnergyMatrix(jvm.toFile(cacheFile))
	else:
		return ematcalc.calcEnergyMatrix()


def AStarTraditional(emat, confSpace):
	builder = c.astar.conf.ConfAStarTree.builder(emat, confSpace)
	builder.setTraditional()
	return builder.build()


def AStarMPLP(emat, confSpace, numIterations=None, convergenceThreshold=None):
	mplpBuilder = c.astar.conf.ConfAStarTree.MPLPBuilder()

	if numIterations is not None:
		mplpBuilder.setNumIterations(numIterations)

	if convergenceThreshold is not None:
		mplpBuilder.setConvergenceThreshold(convergenceThreshold)

	builder = c.astar.conf.ConfAStarTree.builder(emat, confSpace)
	builder.setMPLP(mplpBuilder)
	return builder.build()


def MinimizingEnergyCalculator(confSpace, ffparams=None, parallelism=None, streaming=None):
	builder = c.minimization.SimpleConfMinimizer.builder(confSpace)

	if ffparams is not None:
		builder.setForcefieldParams(ffparams)

	if parallelism is not None:
		builder.setParallelism(parallelism)

	if streaming is not None:
		builder.setStreaming(streaming)

	return builder.build()

def GMECFinder(confSpace, emat=None, astar=None, energyCalculator=None, confLog=None, printIntermediateConfs=None):
	
	if emat is not None:
		builder = c.gmec.SimpleGMECFinder.builder(confSpace, emat)
	elif astar is not None:
		builder = c.gmec.SimpleGMECFinder.builder(confSpace, astar)
	else:
		raise TypeError('either emat or astar must be specified')

	if energyCalculator is not None:
		builder.setEnergyCalculator(energyCalculator)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not None:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	return builder.build()

