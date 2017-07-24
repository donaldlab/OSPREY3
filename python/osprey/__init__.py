
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


def start(heapSizeMB=1024, enableAssertions=False, stackSizeMB=8):
	'''
	Starts the Java Virtual Machine (JVM) that runs Osprey's computation libraries.

	Call :meth:`start` before using any of Osprey's other functions.

	:param int heapSizeMB: Size of the JVM heap in megabytes. This is essentially the amount of memory
		Osprey will have to do computations. 1024 MB is 1 GB, but for larger designs,
		you may want to use 2048 MB (2 GB), 4096 MB (4 GB), or even more memory.
	
	:param bool enableAssertions: pass ``True`` to enable JVM assertions. Only useful for debugging.

	:param int stackSizeMB: Size of the JVM stack in megabytes. Generally leave this at the default value,
		unless Osprey crashes because it's too small. Then try increasing it.
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
	jvm.start(heapSizeMB, enableAssertions, stackSizeMB)

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


def writePdb(mol, path):
	'''
	save a molecule to a PDB file

	:param mol: the molecule to save
	:type mol: :java:ref:`.structure.Molecule`

	:param str path: path of the PDB file
	'''

	# unbox the mol if we need to
	if isinstance(mol, c.confspace.ParametricMolecule):
		mol = mol.mol

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
	:type residues: [int, int]

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


def EnergyCalculator(confSpace, ffparams, parallelism=None, type=None):
	'''
	:java:classdoc:`.energy.EnergyCalculator`

	:param confSpace: The conformation space containing the residue templates to use for atom connectivities.
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:builder_option ffparams .energy.EnergyCalculator$Builder#ffparams:
	:builder_option parallelism .energy.EnergyCalculator$Builder#parallelism:
	:builder_option type .energy.EnergyCalculator$Builder#type:
	:builder_return .energy.EnergyCalculator$Builder:
	'''
	builder = _get_builder(c.energy.EnergyCalculator)(confSpace, ffparams)

	if parallelism is not None:
		builder.setParallelism(parallelism)

	if type is not None:
		builder.setType(type)

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


def AStarTraditional(emat, confSpaceOrPmat, useExternalMemory=False):
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
	builder.setMPLP(mplpBuilder)

	if useExternalMemory == True:
		builder.useExternalMemory()

	return builder.build()


def GMECFinder(astar, confEcalc, confLog=None, printIntermediateConfs=None, useExternalMemory=None):
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

	builder = _get_builder(c.gmec.SimpleGMECFinder)(astar, confEcalc)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not None:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	if useExternalMemory == True:
		builder.useExternalMemory()

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


def KStar(protein, ligand, complexConfSpace, ecalc, confEcalcFactory, astarFactory, epsilon=None, maxSimultaneousMutations=None):
	# TODO: docstring

	# convert functions from python to java
	confEcalcFactory = jpype.JProxy(jvm.getInnerClass(c.kstar.KStar, 'ConfEnergyCalculatorFactory'), dict={ 'make': confEcalcFactory })
	astarFactory = jpype.JProxy(c.gmec.ConfSearchFactory, dict={ 'make': astarFactory })

	# build settings
	settingsBuilder = _get_builder(jvm.getInnerClass(c.kstar.KStar, 'Settings'))()
	if epsilon is not None:
		settingsBuilder.setEpsilon(epsilon)
	if maxSimultaneousMutations is not None:
		settingsBuilder.setMaxSimultaneousMutations(maxSimultaneousMutations)
	settings = settingsBuilder.build()

	kstar = c.kstar.KStar(protein, ligand, complexConfSpace, ecalc, confEcalcFactory, astarFactory, settings)

	# TEMP
	#score = kstar.calcConfSpaceScore()
	#print('K* score: %e' % score.doubleValue())
	kstar.run()
