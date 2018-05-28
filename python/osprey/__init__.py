
import sys, os, jpype
import jvm, wraps


# NOTE: These classes exists only to talk to the docstring processors.
# Variables that get assigned instances of these classes will
# be re-assigned with other values once the JVM is started
class DocstringVal:
	def __init__(self, val):
		self.val = val
	def __repr__(self):
		return self.val

class DocstringJavaDefault(DocstringVal):
	def __init__(self, val):
		DocstringVal.__init__(self, val)

c = None

WILD_TYPE = DocstringJavaDefault('.confspace.Strand#WildType')
'''Magic constant that refers to the wild-type template of the residue'''

Forcefield = None
SolvationForcefield = None
EnergyPartition = None
ExternalMemory = None
ConfSpaceType = None
BreakdownType = None

# make a special type to use in function signatures to explicitly
# signal that values should rely on defaults in the java code
class UseJavaDefault:
	pass
useJavaDefault = UseJavaDefault()


def _get_builder(jclass, builder_name='Builder'):
	return jvm.getInnerClass(jclass, builder_name)


def _java_aware_excepthook(exctype, value, traceback):

	# show original python exception info
	sys.__excepthook__(exctype, value, traceback)

	# try to print java exception info
	try:
		print('\n%s' % value.stacktrace())
	except (AttributeError, TypeError):
		# must not be a java exception
		pass


def start(heapSizeMiB=1024, enableAssertions=False, stackSizeMiB=16, garbageSizeMiB=128):
	'''
	Starts the Java Virtual Machine (JVM) that runs Osprey's computation libraries.

	Call :meth:`start` before using any of Osprey's other functions.

	:param int heapSizeMiB: Size of the JVM heap in megabytes. This is essentially the amount of memory
		Osprey will have to do computations. 1024 MiB is 1 GiB, but for larger designs,
		you may want to use 2048 MiB (2 GiB), 4096 MiB (4 GiB), or even more memory.
	
	:param bool enableAssertions: pass ``True`` to enable JVM assertions. Only useful for debugging.

	:param int stackSizeMiB: Size of the JVM stack portion of the heap in megabytes.
		Generally leave this at the default value, unless Osprey crashes because it's too small. Then try increasing it.

	:param int garbageSizeMiB: Size of the garbage portion of the JVM heap that is reserved for temporary objects.
		This default value is appropriate for the default heap size, but if using larger heap sizes, then increasing
		the garbage size to 256, 512, or even 1024 MiB can give a modest improvement in performance.
	'''

	# disable buffered output on stdout, so python log messages line up with java log messages
	sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)

	# setup a global exception handler to show java exception info
	sys.excepthook = _java_aware_excepthook

	# is this the real life, or is this just fantasy?
	lib_dir = os.path.join(os.path.dirname(__file__), 'lib')
	no_escape_from_reality = os.path.exists(lib_dir)

	# build jvm classpath
	if no_escape_from_reality:

		# release environment: use jars in lib folder
		for file in os.listdir(lib_dir):
			jvm.addClasspath(os.path.join(lib_dir, file))

	else:

		# development environment: use the gradle-defined classpath
		osprey_dir = os.path.join(os.path.dirname(__file__), '../../')
		classpath_path = os.path.join(osprey_dir, 'build/python/classpath.txt')
		if not os.path.isfile(classpath_path):
			raise Exception('dev classpath for python not generated yet. run ./gradlew pythonDevelop')

		for path in open(classpath_path, 'r').readlines():
			jvm.addClasspath(path.strip())

	# start the jvm
	jvm.start(heapSizeMiB, enableAssertions, stackSizeMiB, garbageSizeMiB)

	# set up class factories
	global c
	c = jpype.JPackage('edu.duke.cs.osprey')

	wraps.init(c)

	# init other globals
	global WILD_TYPE
	WILD_TYPE = c.confspace.Strand.WildType
	global Forcefield
	Forcefield = jvm.getInnerClass(c.energy.forcefield.ForcefieldParams, 'Forcefield')
	global SolvationForcefield
	SolvationForcefield = jvm.getInnerClass(c.energy.forcefield.ForcefieldParams, 'SolvationForcefield')
	global EnergyPartition
	EnergyPartition = c.energy.EnergyPartition
	global ExternalMemory
	ExternalMemory = c.externalMemory.ExternalMemory
	global ConfSpaceType
	ConfSpaceType = jvm.getInnerClass(c.kstar.KStar, 'ConfSpaceType')
	global BreakdownType
	BreakdownType = jvm.getInnerClass(c.energy.ResidueForcefieldBreakdown, 'Type')

	# expose static builder methods too
	Parallelism.makeCpu = c.parallelism.Parallelism.makeCpu
	Parallelism.make = c.parallelism.Parallelism.make

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


def writePdb(mol, path, comment=None):
	'''
	save a molecule to a PDB file

	:param mol: the molecule to save
	:type mol: :java:ref:`.structure.Molecule`
	:param str path: path of the PDB file
	:param str comment: Optional comment to add to the PDB file headers
	'''

	# deal with different molecule types
	if isinstance(mol, jvm.getInnerClass(c.energy.EnergyCalculator, 'EnergiedParametricMolecule')):
		c.structure.PDBIO.writeFile(mol, comment, path)
	elif isinstance(mol, c.confspace.ParametricMolecule):
		c.structure.PDBIO.writeFile(mol.mol, comment, None, path)
	else:
		c.structure.PDBIO.writeFile(mol, path)

	print('write PDB file to file: %s' % path)


def printGpuInfo():
	'''
	Prints information about GPU hardware present in the system
	and which GPUs are usable by Osprey.

	Both Cuda and OpenCL APIs are queried.
	'''
	c.gpu.cuda.Diagnostics.main(None)
	c.gpu.opencl.Diagnostics.main(None)


def initExternalMemory(internalSizeMiB, tempDir=None, tempSubdir=None):
	'''
	Initializes external memory for calculations.

	:param int internalSizeMiB: :java:methoddoc:`.externalMemory.ExternalMemory#setInternalLimit`
	:param str tempDir: Path to temporary directory to host external memory
	:default tempDir: <system temp dir>
	:param str tempSubdir: name of subdirectory within tempDir
	:default tempSubdir: <automatically generated>
	'''

	ExternalMemory.setInternalLimit(internalSizeMiB)
	if tempDir is not None:
		if tempSubdir is not None:
			ExternalMemory.setTempDir(tempDir, tempSubdir)
		else:
			ExternalMemory.setTempDir(tempDir)

	# make a function we can call to cleanup the external memory
	def cleanupExternalMemory():
		ExternalMemory.cleanup()
		pass

	# automatically cleanup the external memory at exit
	import atexit
	atexit.register(cleanupExternalMemory)

	# and cleanup at sigint (e.g., ctrl-c) too
	import signal
	def sigint_handler(signal, frame):
		cleanupExternalMemory()
	signal.signal(signal.SIGINT, sigint_handler)


#-------------------------------------#
# pythonic wrappers for Java builders #
#-------------------------------------#

def Parallelism(cpuCores=None, gpus=None, streamsPerGpu=None):
	'''
	:java:classdoc:`.parallelism.Parallelism`

	:builder_option cpuCores .parallelism.Parallelism$Builder#numCpus:
	:builder_option gpus .parallelism.Parallelism$Builder#numGpus:
	:builder_option streamsPerGpu .parallelism.Parallelism$Builder#numStreamsPerGpu:
	:builder_return .parallelism.Parallelism$Builder:
	'''
	builder = _get_builder(c.parallelism.Parallelism)()
	if cpuCores is not None:
		builder.setNumCpus(cpuCores)
	if gpus is not None:
		builder.setNumGpus(gpus)
	if streamsPerGpu is not None:
		builder.setNumStreamsPerGpu(streamsPerGpu)

	return builder.build()


def TemplateLibrary(
		forcefield=None,
		defaultTemplates=True, extraTemplates=[],
		defaultTemplateCoords=True, extraTemplateCoords=[],
		defaultRotamers=True, extraRotamers=[],
		extraBackboneDependentRotamers=[],
		defaultResidueEntropies=True, extraResidueEntropies=[],
		makeDAminoAcids=None,
		moleculesForWildTypeRotamers=[]
	):
	'''
	:java:classdoc:`.restypes.ResidueTemplateLibrary`

	:builder_option forcefield .restypes.ResidueTemplateLibrary$Builder#forcefield:

	:param bool defaultTemplates: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#clearTemplates`
	:param extraTemplates: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addTemplates`
	:type extraTemplates: [template string or file path]

	:param bool defaultTemplateCoords: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#clearTemplateCoords`
	:param extraTemplateCoords: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addTemplateCoords`
	:type extraTemplateCoords: [coords string or file path]

	:param bool defaultRotamers: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#clearRotamers`
	:param extraRotamers: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addRotamers`
	:type extraRotamers: [rotamers string or file path]

	:param extraBackboneDependentRotamers: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addBackboneDependentRotamers`
	:type extraBackboneDependentRotamers: [backbone-dependent rotamers string or file path]

	:builder_option makeDAminoAcids .restypes.ResidueTemplateLibrary$Builder#makeDAminoAcidTemplates:

	:param moleculesForWildTypeRotamers: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addMoleculeForWildTypeRotamers`
	:type moleculesForWildTypeRotamers: [:java:ref:`.structure.Molecule`]

	:builder_return .restypes.ResidueTemplateLibrary$Builder:
	'''

	if forcefield is None:
		builder = _get_builder(c.restypes.ResidueTemplateLibrary)()
	else:
		builder = _get_builder(c.restypes.ResidueTemplateLibrary)(forcefield)

	if not defaultTemplates:
		builder.clearTemplates()
	for text in [_ensureText(pathOrText) for pathOrText in extraTemplates]:
		builder.addTemplates(text)

	if not defaultTemplateCoords:
		builder.clearTemplateCoords()
	for text in [_ensureText(pathOrText) for pathOrText in extraTemplateCoords]:
		builder.addTemplateCoords(text)

	if not defaultRotamers:
		builder.clearRotamers()
	for text in [_ensureText(pathOrText) for pathOrText in extraRotamers]:
		builder.addRotamers(text)

	for text in [_ensureText(pathOrText) for pathOrText in extraBackboneDependentRotamers]:
		builder.addBackboneDependentRotamers(text)

	if makeDAminoAcids is not None:
		builder.setMakeDAminoAcidTemplates(makeDAminoAcids)

	for mol in moleculesForWildTypeRotamers:
		builder.addMoleculeForWildTypeRotamers(mol)

	return builder.build()
	

def Strand(pathOrMol, residues=None, templateLib=None):
	'''
	:java:classdoc:`.confspace.Strand`

	:param pathOrMol: path to a PDB file, or a molecule instance
	:type pathOrMol: str or :java:ref:`.structure.Molecule`

	:param residues: range of residue numbers, inclusive. `None` to include all residues.
	:type residues: [str, str]

	:builder_option templateLib .confspace.Strand$Builder#templateLib:

	:builder_return .confspace.Strand$Builder:
	'''

	# did we get a path or a molecule?
	if isinstance(pathOrMol, c.structure.Molecule):
		mol = pathOrMol
	else:
		mol = readPdb(pathOrMol)

	builder = _get_builder(c.confspace.Strand)(mol)

	if residues is not None:
		builder.setResidues(residues[0], residues[1])

	if templateLib is not None:
		builder.setTemplateLibrary(templateLib)
	
	return builder.build()


def ConfSpace(strands, shellDist=None):
	'''
	:java:classdoc:`.confspace.SimpleConfSpace`

	:param strands: the strands to use
	:type strands: :java:ref:`.confspace.Strand` or list of Strands
	:builder_option shellDist .confspace.SimpleConfSpace$Builder#shellDist:
	:builder_return .confspace.SimpleConfSpace$Builder:
	'''

	builder = _get_builder(c.confspace.SimpleConfSpace)()

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

		builder.addStrand(strand, flex)

	# add the shell distance if needed
	if shellDist is not None:
		builder.setShellDistance(shellDist)
	
	return builder.build()


def StrandFlex():
	# TODO: implement me
	pass


def ForcefieldParams(forcefield=None):
	'''
	:java:classdoc:`.energy.forcefield.ForcefieldParams`
	
	Configure the forcefield parameters by setting the properties of the :java:ref:`.energy.forcefield.ForcefieldParams` object.

	:builder_option forcefield .energy.forcefield.ForcefieldParams#forcefld:
	:rtype: :java:ref:`.energy.forcefield.ForcefieldParams`
	'''
	return c.energy.forcefield.ForcefieldParams()


def EnergyCalculator(confSpace, ffparams, parallelism=None, type=None, isMinimizing=None, infiniteWellEnergy=None):
	'''
	:java:classdoc:`.energy.EnergyCalculator`

	:param confSpace: The conformation space containing the residue templates to use for atom connectivities.
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:builder_option ffparams .energy.EnergyCalculator$Builder#ffparams:
	:builder_option parallelism .energy.EnergyCalculator$Builder#parallelism:
	:builder_option type .energy.EnergyCalculator$Builder#type:
	:builder_option isMinimizing .energy.EnergyCalculator$Builder#isMinimizing:
	:builder_option infiniteWellEnergy .energy.EnergyCalculator$Builder#infiniteWellEnergy:

	:builder_return .energy.EnergyCalculator$Builder:
	'''
	builder = _get_builder(c.energy.EnergyCalculator)(confSpace, ffparams)

	if parallelism is not None:
		builder.setParallelism(parallelism)

	if type is not None:
		builder.setType(type)

	if isMinimizing is not None:
		builder.setIsMinimizing(isMinimizing)

	if infiniteWellEnergy is not None:
		builder.setInfiniteWellEnergy(jvm.boxDouble(infiniteWellEnergy))

	return builder.build()


def SharedEnergyCalculator(ecalc, isMinimizing=None):
	'''
	:java:classdoc:`.energy.EnergyCalculator$SharedBuilder`

	:param ecalc: The existing energy calculator with which to share resources
	:type ecalc: :java:ref:`.energy.EnergyCalculator`
	:builder_option isMinimizing .energy.EnergyCalculator$Builder#isMinimizing:
	:builder_return .energy.EnergyCalculator$SharedBuilder:
	'''
	builder = jvm.getInnerClass(c.energy.EnergyCalculator, 'SharedBuilder')(ecalc)

	if isMinimizing is not None:
		builder.setIsMinimizing(isMinimizing)

	return builder.build()


def ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=None, addResEntropy=None, energyPartition=None):
	'''
	:java:classdoc:`.energy.ConfEnergyCalculator`

	:builder_option confSpace .energy.ConfEnergyCalculator$Builder#confSpace:
	:builder_option ecalc .energy.ConfEnergyCalculator$Builder#ecalc:
	:builder_option referenceEnergies .energy.ConfEnergyCalculator$Builder#eref:
	:builder_option addResEntropy .energy.ConfEnergyCalculator$Builder#addResEntropy:
	:builder_option energyPartition .energy.ConfEnergyCalculator$Builder#epart:
	:builder_return .energy.ConfEnergyCalculator$Builder:
	'''
	builder = _get_builder(c.energy.ConfEnergyCalculator)(confSpace, ecalc)

	if referenceEnergies is not None:
		builder.setReferenceEnergies(referenceEnergies)

	if energyPartition is not None:
		builder.setEnergyPartition(energyPartition)

	return builder.build()


def EnergyMatrix(confEcalc, cacheFile=None):
	'''
	:java:methoddoc:`.ematrix.SimplerEnergyMatrixCalculator#calcEnergyMatrix`

	:builder_option confEcalc .ematrix.SimplerEnergyMatrixCalculator$Builder#confEcalc:
	:builder_option cacheFile .ematrix.SimplerEnergyMatrixCalculator$Builder#cacheFile:
	'''
	
	builder = _get_builder(c.ematrix.SimplerEnergyMatrixCalculator)(confEcalc)

	if cacheFile is not None:
		builder.setCacheFile(jvm.toFile(cacheFile))

	return builder.build().calcEnergyMatrix()


def ReferenceEnergies(confSpace, ecalc, addResEntropy=None):
	'''
	:java:methoddoc:`.ematrix.SimplerEnergyMatrixCalculator#calcReferenceEnergies`

	:builder_option confSpace .ematrix.SimpleReferenceEnergies$Builder#confSpace:
	:builder_option ecalc .ematrix.SimpleReferenceEnergies$Builder#ecalc:
	:builder_option addResEntropy .ematrix.SimpleReferenceEnergies$Builder#addResEntropy:
	:builder_return .ematrix.SimpleReferenceEnergies$Builder:
	'''

	builder = _get_builder(c.ematrix.SimpleReferenceEnergies)(confSpace, ecalc)

	if addResEntropy is not None:
		builder.addResEntropy(addResEntropy)

	return builder.build()


def DEE(confSpace, emat, singlesThreshold=None, pairsThreshold=None, singlesGoldsteinDiffThreshold=None, pairsGoldsteinDiffThreshold=None, typeDependent=None, numIterations=None, showProgress=None):
	'''
	:java:classdoc:`.pruning.SimpleDEE$Runner`

	:param confSpace: The design conformation space
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param emat: An energy matrix computed for the conformation space
	:type emat: :java:ref:`.ematrix.EnergyMatrix`

	:builder_option singlesThreshold .pruning.SimpleDEE$Runner#singlesThreshold:
	:builder_option pairsThreshold .pruning.SimpleDEE$Runner#pairsThreshold:
	:builder_option singlesGoldsteinDiffThreshold .pruning.SimpleDEE$Runner#singlesGoldsteinDiffThreshold:
	:builder_option pairsGoldsteinDiffThreshold .pruning.SimpleDEE$Runner#pairsGoldsteinDiffThreshold:
	:builder_option typeDependent .pruning.SimpleDEE$Runner#typeDependent:
	:builder_option numIterations .pruning.SimpleDEE$Runner#numIterations:
	:builder_option showProgress .pruning.SimpleDEE$Runner#showProgress:
	'''

	runner = _get_builder(c.pruning.SimpleDEE, 'Runner')()

	def boxDouble(val):
		return jvm.c.java.lang.Double(val)

	if singlesThreshold is not None:
		runner.setSinglesThreshold(boxDouble(singlesThreshold))
	if pairsThreshold is not None:
		runner.setPairsThreshold(boxDouble(pairsThreshold))
	if singlesGoldsteinDiffThreshold is not None:
		runner.setSinglesGoldsteinDiffThreshold(boxDouble(singlesGoldsteinDiffThreshold))
	if pairsGoldsteinDiffThreshold is not None:
		runner.setPairsGoldsteinDiffThreshold(boxDouble(pairsGoldsteinDiffThreshold))
	if typeDependent is not None:
		runner.setTypeDependent(typeDependent)
	if numIterations is not None:
		runner.setNumIterations(numIterations)
	if showProgress is not None:
		runner.setShowProgress(showProgress)

	return runner.run(confSpace, emat)


def AStarTraditional(emat, confSpaceOrPmat, showProgress=True, useExternalMemory=False):
	'''
	:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#setTraditional`

	:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:
	:param confSpaceOrPmat: The conformation space containing the residue conformations to search.
	:type confSpaceOrPmat: :java:ref:`.confspace.SimpleConfSpace` or :java:ref:`.pruning.PruningMatrix`
	:param useExternalMemory: set to True to use external memory.

		:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#useExternalMemory`

	:type useExternalMemory: boolean
	:builder_return .astar.conf.ConfAStarTree$Builder:
	'''
	builder = _get_builder(c.astar.conf.ConfAStarTree)(emat, confSpaceOrPmat)
	builder.setShowProgress(showProgress)
	builder.setTraditional()

	if useExternalMemory == True:
		builder.useExternalMemory()

	return builder.build()


def EdgeUpdater():
	return c.astar.conf.scoring.mplp.EdgeUpdater()

def NodeUpdater():
	return c.astar.conf.scoring.mplp.NodeUpdater()

def AStarMPLP(emat, confSpaceOrPmat, updater=None, numIterations=None, convergenceThreshold=None, useExternalMemory=False):
	'''
	:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#setMPLP`

	:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:
	:param confSpaceOrPmat: The conformation space containing the residue conformations to search.
	:type confSpaceOrPmat: :java:ref:`.confspace.SimpleConfSpace` or :java:ref:`.pruning.PruningMatrix`
	:builder_option updater .astar.conf.ConfAStarTree$MPLPBuilder#updater:
	:builder_option numIterations .astar.conf.ConfAStarTree$MPLPBuilder#numIterations:
	:builder_option convergenceThreshold .astar.conf.ConfAStarTree$MPLPBuilder#convergenceThreshold:
	:param useExternalMemory: set to True to use external memory.

		:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#useExternalMemory`

	:type useExternalMemory: boolean
	:builder_return .astar.conf.ConfAStarTree$Builder:
	'''
	mplpBuilder = _get_builder(c.astar.conf.ConfAStarTree, 'MPLPBuilder')()

	if updater is not None:
		mplpBuilder.setUpdater(updater)

	if numIterations is not None:
		mplpBuilder.setNumIterations(numIterations)

	if convergenceThreshold is not None:
		mplpBuilder.setConvergenceThreshold(convergenceThreshold)

	builder = _get_builder(c.astar.conf.ConfAStarTree)(emat, confSpaceOrPmat)
	builder.setShowProgress(True)
	builder.setMPLP(mplpBuilder)

	if useExternalMemory == True:
		builder.useExternalMemory()

	return builder.build()


def GMECFinder(astar, confEcalc, confLog=None, printIntermediateConfs=None, useExternalMemory=None, resumeLog=None, confDBFile=None):
	'''
	:java:classdoc:`.gmec.SimpleGMECFinder`

	:builder_option astar .gmec.SimpleGMECFinder$Builder#search:

		Use one of :py:func:`AStarTraditional` or :py:func:`AStarMPLP` to get an A* implementation.

	:builder_option confEcalc .gmec.SimpleGMECFinder$Builder#confEcalc:

		Use :py:func:`ConfEnergyCalculator` to get a conformation energy calculator.

	:param str confLog: Path to file where conformations found during conformation space search should be logged.
	:builder_option printIntermediateConfs .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole:
	:builder_option useExternalMemory .gmec.SimpleGMECFinder$Builder#useExternalMemory:
	:param str resumeLog: Path to log file where resume info will be written or read, so designs can be resumed.
	:builder_return .gmec.SimpleGMECFinder$Builder:
	'''

	builder = _get_builder(c.gmec.SimpleGMECFinder)(astar, confEcalc)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not None:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	if useExternalMemory == True:
		builder.useExternalMemory()

	if resumeLog is not None:
		builder.setResumeLog(jvm.toFile(resumeLog))

	if confDBFile is not None:
		builder.setConfDB(jvm.toFile(confDBFile))

	return builder.build()


def DEEGMECFinder(emat, confSpace, ecalc, confEcalc, name, use_epic, use_lute, confLog=None, printIntermediateConfs=None, useExternalMemory=None):
	'''
	:java:classdoc:`.gmec.SimpleGMECFinder`

	:builder_option astar .gmec.SimpleGMECFinder$Builder#search:

			Use one of :py:func:`AStarTraditional` or :py:func:`AStarMPLP` to get an A* implementation.

	:builder_option confEcalc .gmec.SimpleGMECFinder$Builder#confEcalc:

			Use :py:func:`ConfEnergyCalculator` to get a conformation energy calculator.

	:param str confLog: Path to file where conformations found during conformation space search should be logged.
	:builder_option printIntermediateConfs .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole:
	:builder_option useExternalMemory .gmec.SimpleGMECFinder$Builder#useExternalMemory:
	:builder_return .gmec.SimpleGMECFinder$Builder:
	'''

	builder = _get_builder(c.gmec.DEEGMECFinder)(emat, confSpace, ecalc, confEcalc, name)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not None:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	if useExternalMemory == True:
		builder.useExternalMemory()

	gf = builder.build()

	if use_epic:
		gf.epicSettings = c.ematrix.epic.EPICSettings.defaultEPIC()
		gf.pruningSettings.algOption = 3
	if use_lute:
		gf.luteSettings = c.tupexp.LUTESettings.defaultLUTE()
		gf.pruningSettings.algOption = 3
		gf.pruningSettings.useTriples = True

	return gf


def DEEPerStrandFlex(strand, pert_file_name, flex_res_list, pdb_file):
	deeper_settings = c.dof.deeper.DEEPerSettings(True, pert_file_name, True, 'None', False, 2.5, 2.5, False, jvm.toArrayList(flex_res_list), pdb_file, False, strand.templateLib)
	bbflex = c.confspace.DEEPerStrandFlex(strand,deeper_settings)
	return bbflex


def KStar(proteinConfSpace, ligandConfSpace, complexConfSpace, ecalc, confEcalcFactory, astarFactory, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, energyMatrixCachePattern=useJavaDefault, confDBPattern=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None):
	'''
	:java:classdoc:`.kstar.KStar`

	For examples using K*, see the examples/python.KStar directory in your Osprey distribution.

	:param proteinConfSpace: :java:fielddoc:`.kstar.KStar#protein`
	:type proteinConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param ligandConfSpace: :java:fielddoc:`.kstar.KStar#ligand`
	:type ligandConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param complexConfSpace: :java:fielddoc:`.kstar.KStar#complex`
	:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param ecalc: :java:fielddoc:`.kstar.KStar#ecalc`
	:type ecalc: :java:ref:`.energy.EnergyCalculator`
	:param confEcalcFactory: :java:fielddoc:`.kstar.KStar#confEcalcFactory`
	:type confEcalcFactory: :java:ref:`.kstar.KStar$ConfEnergyCalculatorFactory`
	:param astarFactory: :java:fielddoc:`.kstar.KStar#confSearchFactory`
	:type astarFactory: :java:ref:`.kstar.KStar$ConfSearchFactory`
	:builder_option epsilon .kstar.KStar$Settings$Builder#epsilon:
	:builder_option stabilityThreshold .kstar.KStar$Settings$Builder#stabilityThreshold:
	:builder_option maxSimultaneousMutations .kstar.KStar$Settings$Builder#maxSimultaneousMutations:
	:builder_option energyMatrixCachePattern .kstar.KStar$Settings$Builder#energyMatrixCachePattern:
	:builder_option confDBPattern .kstar.KStar$Settings$Builder#confDBPattern:
	:param bool writeSequencesToConsole: True to write sequences and scores to the console
	:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

	:rtype: :java:ref:`.kstar.KStar`
	'''

	# convert functions from python to java
	confEcalcFactory = jpype.JProxy(jvm.getInnerClass(c.kstar.KStar, 'ConfEnergyCalculatorFactory'), dict={ 'make': confEcalcFactory })
	astarFactory = jpype.JProxy(jvm.getInnerClass(c.kstar.KStar, 'ConfSearchFactory'), dict={ 'make': astarFactory })

	# build settings
	settingsBuilder = _get_builder(jvm.getInnerClass(c.kstar.KStar, 'Settings'))()
	if epsilon is not useJavaDefault:
		settingsBuilder.setEpsilon(epsilon)
	if stabilityThreshold is not useJavaDefault:
		settingsBuilder.setStabilityThreshold(jvm.boxDouble(stabilityThreshold))
	if maxSimultaneousMutations is not useJavaDefault:
		settingsBuilder.setMaxSimultaneousMutations(maxSimultaneousMutations)
	if writeSequencesToConsole:
		settingsBuilder.addScoreConsoleWriter()
	if writeSequencesToFile is not None:
		settingsBuilder.addScoreFileWriter(jvm.toFile(writeSequencesToFile))
	if energyMatrixCachePattern is not useJavaDefault:
		settingsBuilder.setEnergyMatrixCachePattern(energyMatrixCachePattern)
	if confDBPattern is not useJavaDefault:
		settingsBuilder.setConfDBPattern(confDBPattern)
	settings = settingsBuilder.build()

	return c.kstar.KStar(proteinConfSpace, ligandConfSpace, complexConfSpace, ecalc, confEcalcFactory, astarFactory, settings)


def BBKStar(proteinConfSpace, ligandConfSpace, complexConfSpace, rigidEcalc, minimizingEcalc, confEcalcFactory, astarFactory, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, energyMatrixCachePattern=useJavaDefault, confDBPattern=useJavaDefault, numBestSequences=useJavaDefault, numConfsPerBatch=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None):
	'''
	:java:classdoc:`.kstar.BBKStar`

	For examples using BBK*, see the examples/python.KStar directory in your Osprey distribution.

	:param proteinConfSpace: :java:fielddoc:`.kstar.BBKStar#protein`
	:type proteinConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param ligandConfSpace: :java:fielddoc:`.kstar.BBKStar#ligand`
	:type ligandConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param complexConfSpace: :java:fielddoc:`.kstar.BBKStar#complex`
	:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param rigidEcalc: :java:fielddoc:`.kstar.BBKStar#rigidEcalc`
	:type rigidEcalc: :java:ref:`.energy.EnergyCalculator`
	:param minimizingEcalc: :java:fielddoc:`.kstar.BBKStar#minimizingEcalc`
	:type minimizingEcalc: :java:ref:`.energy.EnergyCalculator`
	:param confEcalcFactory: :java:fielddoc:`.kstar.BBKStar#confEcalcFactory`
	:type confEcalcFactory: :java:ref:`.kstar.KStar$ConfEnergyCalculatorFactory`
	:param astarFactory: :java:fielddoc:`.kstar.BBKStar#confSearchFactory`
	:type astarFactory: :java:ref:`.kstar.KStar$ConfSearchFactory`
	:builder_option epsilon .kstar.KStar$Settings$Builder#epsilon:
	:builder_option stabilityThreshold .kstar.KStar$Settings$Builder#stabilityThreshold:
	:builder_option maxSimultaneousMutations .kstar.KStar$Settings$Builder#maxSimultaneousMutations:
	:builder_option energyMatrixCachePattern .kstar.KStar$Settings$Builder#energyMatrixCachePattern:
	:builder_option confDBPattern .kstar.KStar$Settings$Builder#confDBPattern:
	:builder_option numBestSequences .kstar.BBKStar$Settings$Builder#numBestSequences:
	:builder_option numConfsPerBatch .kstar.BBKStar$Settings$Builder#numConfsPerBatch:
	:param bool writeSequencesToConsole: True to write sequences and scores to the console
	:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

	:rtype: :java:ref:`.kstar.BBKStar`
	'''

	# convert functions from python to java
	confEcalcFactory = jpype.JProxy(jvm.getInnerClass(c.kstar.KStar, 'ConfEnergyCalculatorFactory'), dict={ 'make': confEcalcFactory })
	astarFactory = jpype.JProxy(jvm.getInnerClass(c.kstar.KStar, 'ConfSearchFactory'), dict={ 'make': astarFactory })

	# build settings
	kstarSettingsBuilder = _get_builder(jvm.getInnerClass(c.kstar.KStar, 'Settings'))()
	if epsilon is not useJavaDefault:
		kstarSettingsBuilder.setEpsilon(epsilon)
	if stabilityThreshold is not useJavaDefault:
		kstarSettingsBuilder.setStabilityThreshold(jvm.boxDouble(stabilityThreshold))
	if maxSimultaneousMutations is not useJavaDefault:
		kstarSettingsBuilder.setMaxSimultaneousMutations(maxSimultaneousMutations)
	if writeSequencesToConsole:
		kstarSettingsBuilder.addScoreConsoleWriter()
	if writeSequencesToFile is not None:
		kstarSettingsBuilder.addScoreFileWriter(jvm.toFile(writeSequencesToFile))
	if energyMatrixCachePattern is not useJavaDefault:
		kstarSettingsBuilder.setEnergyMatrixCachePattern(energyMatrixCachePattern)
	if confDBPattern is not useJavaDefault:
		kstarSettingsBuilder.setConfDBPattern(confDBPattern)
	kstarSettings = kstarSettingsBuilder.build()

	bbkstarSettingsBuilder = _get_builder(jvm.getInnerClass(c.kstar.BBKStar, 'Settings'))()
	if numBestSequences is not useJavaDefault:
		bbkstarSettingsBuilder.setNumBestSequences(numBestSequences)
	if numConfsPerBatch is not useJavaDefault:
		bbkstarSettingsBuilder.setNumConfsPerBatch(numConfsPerBatch)
	bbkstarSettings = bbkstarSettingsBuilder.build()

	return c.kstar.BBKStar(proteinConfSpace, ligandConfSpace, complexConfSpace, rigidEcalc, minimizingEcalc, confEcalcFactory, astarFactory, kstarSettings, bbkstarSettings)


def ConfAnalyzer(confEcalc, emat):
	'''
	:java:classdoc:`.gmec.ConfAnalyzer`

	For examples using the conf analyzer, see examples/python.GMEC/analyzeConf.py in your Osprey distribution.

	:param confEcalc: :java:fielddoc:`.gmec.SimpleGMECFinder$Builder#confEcalc`
	:type confEcalc: :java:ref:`.energy.ConfEnergyCalculator`
	:param emat: The energy matrix that defines conformation scores
	:type emat: :java:ref:`.ematrix.EnergyMatrix`

	:rtype: :java:ref:`.gmec.ConfAnalyzer`
	'''

	return c.gmec.ConfAnalyzer(confEcalc, emat)


def SequenceAnalyzer(proteinConfSpace, ligandConfSpace, complexConfSpace, ecalc, confEcalcFactory, astarFactory, energyMatrixCachePattern=useJavaDefault, confDBPattern=useJavaDefault):
	'''
	:java:classdoc:`.kstar.SequenceAnalyzer`

	For examples using the sequence analyzer, see examples/python.KStar/analyzeSequence.py in your Osprey distribution.

	:param proteinConfSpace: :java:fielddoc:`.kstar.KStar#protein`
	:type proteinConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param ligandConfSpace: :java:fielddoc:`.kstar.KStar#ligand`
	:type ligandConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param complexConfSpace: :java:fielddoc:`.kstar.KStar#complex`
	:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param ecalc: :java:fielddoc:`.kstar.KStar#ecalc`
	:type ecalc: :java:ref:`.energy.EnergyCalculator`
	:param confEcalcFactory: :java:fielddoc:`.kstar.KStar#confEcalcFactory`
	:type confEcalcFactory: :java:ref:`.kstar.KStar$ConfEnergyCalculatorFactory`
	:param astarFactory: :java:fielddoc:`.kstar.KStar#confSearchFactory`
	:type astarFactory: :java:ref:`.kstar.KStar$ConfSearchFactory`
	:builder_option energyMatrixCachePattern .kstar.KStar$Settings$Builder#energyMatrixCachePattern:

	:rtype: :java:ref:`.kstar.SequenceAnalyzer`
	'''

	# convert functions from python to java
	confEcalcFactory = jpype.JProxy(jvm.getInnerClass(c.kstar.KStar, 'ConfEnergyCalculatorFactory'), dict={ 'make': confEcalcFactory })
	astarFactory = jpype.JProxy(jvm.getInnerClass(c.kstar.KStar, 'ConfSearchFactory'), dict={ 'make': astarFactory })

	# build settings
	settingsBuilder = _get_builder(jvm.getInnerClass(c.kstar.KStar, 'Settings'))()
	if energyMatrixCachePattern is not useJavaDefault:
		settingsBuilder.setEnergyMatrixCachePattern(energyMatrixCachePattern)
	if confDBPattern is not useJavaDefault:
		settingsBuilder.setConfDBPattern(confDBPattern)
	settings = settingsBuilder.build()

	return c.kstar.SequenceAnalyzer(proteinConfSpace, ligandConfSpace, complexConfSpace, ecalc, confEcalcFactory, astarFactory, settings)
