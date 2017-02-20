
import sys, os, jpype
import jvm, wraps

getJavaClass = wraps.getJavaClass


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
	'''
	Starts the Java Virtual Machine (JVM) that runs Osprey's computation libraries.

	Call :meth:`start` before using any of Osprey's other functions.

	:param int heapSizeMB: Size of the JVM heap in megabytes. This is essentially the amount of memory
		Osprey will have to do computations. 1024 MB is 1 GB, but for larger designs,
		you may want to use 2048 MB (2 GB), 4096 MB (4 GB), or even more memory.
	
	:param bool enableAssertions: pass ``True`` to enable JVM assertions. Only useful for debugging.
	'''

	# setup a global exception handler to show java exception info
	sys.excepthook = _javaAwareExcepthook
	
	# start the jvm
	for path in _ospreyPaths:
		jvm.addClasspath(path)
	jvm.start(heapSizeMB, enableAssertions)

	# set up class factories
	global c
	c = jpype.JPackage('edu.duke.cs.osprey')

	wraps.init()

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
	'''
	Useful for passing file contents to functions that expect strings rather than file paths

	:returns: the text of the file at ``path``
	:rtype: str
	'''
	return c.tools.FileTools.readFile(path)


def _ensureText(textOrPath):
	if os.path.exists(textOrPath):
		# it's a path, read the file
		return readTextFile(textOrPath)
	else:
		# it's text, just return it
		return textOrPath


def readPdb(path):
	'''
	loads a PDB file into a molecule object

	.. note:: Raw molecules cannot be used directly in designs.
		Create a :py:meth:`Strand` with the molecule to use it in a design.

	:param str path: path to PDB file
	:rtype: :java:ref:`.structure.Molecule`
	'''
	mol = c.structure.PDBIO.readFile(path)
	print('read PDB file from file: %s' % path)
	return mol


def writePdb(mol, path):
	'''
	save a molecule to a PDB file

	:param molecule: the molecule to save
	:type molecule: :java:ref:`.structure.Molecule`

	:param str path: path of the PDB file
	'''

	# unbox the mol if we need to
	if isinstance(mol, c.confspace.ParametricMolecule):
		mol = mol.mol

	c.structure.PDBIO.writeFile(mol, path)

	print('write PDB file to file: %s' % path)



#-------------------------------------#
# pythonic wrappers for Java builders #
#-------------------------------------#

def Parallelism(cpuCores=None, gpus=None, streamsPerGpu=None):
	'''
	:java:classdoc:`.parallelism.Parallelism`

	:builder_option cpuCores .parallelism.Parallelism#numThreads:
	:builder_option gpus .parallelism.Parallelism#numGpus:
	:builder_option streamsPerGpu .parallelism.Parallelism#numStreamsPerGpu:
	:rtype: :java:ref:`.parallelism.Parallelism`
	'''
	if gpus is not None:
		return c.parallelism.Parallelism.makeGpu(gpus, streamsPerGpu)
	elif cpuCores is not None:
		return c.parallelism.Parallelism.makeCpu(cpuCores)
	else:
		return c.parallelism.Parallelism.makeDefault()


def TemplateLibrary(forcefield=None, templateCoords=None, rotamers=None, backboneDependentRotamers=None):
	'''
	:java:classdoc:`.restypes.GenericResidueTemplateLibrary`

	:builder_option forcefield .restypes.GenericResidueTemplateLibrary$Builder#forcefield:
	:default forcefield: osprey.Forcefield.AMBER

	:builder_option templateCoords .restypes.GenericResidueTemplateLibrary$Builder#templateCoordsText:
	:type templateCoords: coords str or file path

	:builder_option rotamers .restypes.GenericResidueTemplateLibrary$Builder#rotamersText:
	:type rotamers: coords str or file path

	:builder_option backboneDependentRotamers .restypes.GenericResidueTemplateLibrary$Builder#backboneDependentRotamersText:
	:type backboneDependentRotamers: coords str or file path

	:builder_return .restypes.GenericResidueTemplateLibrary$Builder:
	'''

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
	

def Strand(pathOrMol, residues=None):
	'''
	:java:classdoc:`.confspace.Strand`

	:param pathOrMol: path to a PDB file, or a molecule instance
	:type pathOrMol: str or :java:ref:`.structure.Molecule`

	:param residues: range of residue numbers, inclusive. `None` to include all residues.
	:type residues: [int, int]

	:return type: :java:ref:`.confspace.Strand`
	'''

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
	'''
	:java:classdoc:`.confspace.SimpleConfSpace`

	:param strands: the strands to use
	:type strands: :java:ref:`.confspace.Strand` or list of Strands
	:builder_option shellDist .confspace.SimpleConfSpace$Builder#shellDist:
	:builder_return .confspace.SimpleConfSpace$Builder:
	'''

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
	'''
	:java:classdoc:`.energy.forcefield.ForcefieldParams`
	
	Configure the forcefield parameters by setting the properties of the :java:ref:`.energy.forcefield.ForcefieldParams` object.

	:rtype: :java:ref:`.energy.forcefield.ForcefieldParams`
	'''
	return c.energy.forcefield.ForcefieldParams()


def EnergyMatrix(confSpace, ffparams=None, parallelism=None, cacheFile=None):
	'''
	:java:classdoc:`.ematrix.SimplerEnergyMatrixCalculator`

	:builder_option confSpace .ematrix.SimplerEnergyMatrixCalculator$Builder#confSpace:
	:builder_option ffparams .ematrix.SimplerEnergyMatrixCalculator$Builder#ffparams:
	:builder_option parallelism .ematrix.SimplerEnergyMatrixCalculator$Builder#parallelism:
	:builder_option cacheFile .ematrix.SimplerEnergyMatrixCalculator$Builder#cacheFile:
	'''
	
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
	'''
	Creates an A* search using the traditional estimation function.

	:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:
	:param confSpace: The conformation space containing the residue conformations to search.
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:builder_return .astar.conf.ConfAStarTree$Builder:
	'''
	builder = c.astar.conf.ConfAStarTree.builder(emat, confSpace)
	builder.setTraditional()
	return builder.build()


def AStarMPLP(emat, confSpace, numIterations=None, convergenceThreshold=None):
	'''
	Creates an A* search using the newer estimation function based on Max Product Linear Programming (MPLP)

	:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:
	:param confSpace: The conformation space containing the residue conformations to search.
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:builder_option numIterations .astar.conf.ConfAStarTree$MPLPBuilder#numIterations:
	:builder_option convergenceThreshold .astar.conf.ConfAStarTree$MPLPBuilder#convergenceThreshold:
	:builder_return .astar.conf.ConfAStarTree$Builder:
	'''
	mplpBuilder = c.astar.conf.ConfAStarTree.MPLPBuilder()

	if numIterations is not None:
		mplpBuilder.setNumIterations(numIterations)

	if convergenceThreshold is not None:
		mplpBuilder.setConvergenceThreshold(convergenceThreshold)

	builder = c.astar.conf.ConfAStarTree.builder(emat, confSpace)
	builder.setMPLP(mplpBuilder)
	return builder.build()


def MinimizingEnergyCalculator(confSpace, ffparams=None, parallelism=None, streaming=None):
	'''
	:java:classdoc:`.minimization.SimpleConfMinimizer`

	:builder_option confSpace .minimization.SimpleConfMinimizer$Builder#confSpace:
	:builder_option ffparams .minimization.SimpleConfMinimizer$Builder#ffparams:
	:builder_option parallelism .minimization.SimpleConfMinimizer$Builder#parallelism:
	:builder_option streaming .minimization.SimpleConfMinimizer$Builder#isStreaming:
	:builder_return .minimization.SimpleConfMinimizer$Builder:
	'''
	builder = c.minimization.SimpleConfMinimizer.builder(confSpace)

	if ffparams is not None:
		builder.setForcefieldParams(ffparams)

	if parallelism is not None:
		builder.setParallelism(parallelism)

	if streaming is not None:
		builder.setStreaming(streaming)

	return builder.build()


def GMECFinder(confSpace, emat=None, astar=None, energyCalculator=None, confLog=None, printIntermediateConfs=None):
	'''
	:java:classdoc:`.gmec.SimpleGMECFinder`

	Exactly one of ``emat`` or ``astar`` should be specified.
	If ``emat`` is specified, a default A* implementation will be used and ``astar`` will be ingored.

	:builder_option confSpace .gmec.SimpleGMECFinder$Builder#space:
	:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:
	:builder_option astar .gmec.SimpleGMECFinder$Builder#search:
	:builder_option energyCalculator .gmec.SimpleGMECFinder$Builder#ecalc:
	:param str confLog: Path to file where conformations found during conformation space search should be logged.
	:builder_option printIntermediateConfs .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole:
	:builder_return .gmec.SimpleGMECFinder$Builder:
	'''

	if emat is not None:
		builder = c.gmec.SimpleGMECFinder.builder(confSpace, emat)
	elif astar is not None:
		builder = c.gmec.SimpleGMECFinder.builder(confSpace, astar)
	else:
		raise ValueError('either emat or astar must be specified')

	if energyCalculator is not None:
		builder.setEnergyCalculator(energyCalculator)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not None:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	return builder.build()

