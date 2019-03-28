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

import sys, os, jpype, traceback
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
	def __repr__(self):
		return "(default defined in Java code)"
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


def start(heapSizeMiB=1024, enableAssertions=False, stackSizeMiB=16, garbageSizeMiB=128, allowRemoteManagement=False):
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

	:param bool allowRemoteManagement: pass ``True`` to listen on ports 9010 and 9011 for JMX remote management.
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
	jvm.start(heapSizeMiB, enableAssertions, stackSizeMiB, garbageSizeMiB, allowRemoteManagement)

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
	print("Using up to %d MiB heap memory: %d MiB for garbage, %d MiB for storage" % (
		heapSizeMiB, garbageSizeMiB, heapSizeMiB - garbageSizeMiB
	))


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


class StateMutable(object):
	def __init__(self, name, confSpace):
		'''
		:param str name: unique for the state
		:param confSpace: the conformation space
		:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
		'''
		self.name = name
		self.confSpace = confSpace

class StateUnmutable(object):
	def __init__(self, name, confSpace):
		'''
		:param str name: unique for the state
		:param confSpace: the conformation space
		:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
		'''
		self.name = name
		self.confSpace = confSpace


def MultiStateConfSpace(states):
	'''
	:java:classdoc:`.confspace.MultiStateConfSpace`

	:param states: A list of states to use
	:type states: list of StateMutable and StateUnmutable instances
	'''

	# split out the mutable and unmutable states
	mutableStates = [state for state in states if type(state) is StateMutable]
	unmutableStates = [state for state in states if type(state) is StateUnmutable]

	# start the builder with the first mutable state
	if len(mutableStates) <= 0:
		raise Exception("MultiStateConfSpace needs at least one mutable state")
	builder = _get_builder(c.confspace.MultiStateConfSpace)(mutableStates[0].name, mutableStates[0].confSpace)

	# add the rest of the states
	for state in mutableStates[1:]:
		builder.addMutableState(state.name, state.confSpace)
	for state in unmutableStates:
		builder.addUnmutableState(state.name, state.confSpace)

	return builder.build()


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

	# convert confSpace to a jvm list if possible
	try:
		confSpace = jvm.toArrayList(confSpace)
	except TypeError:
		# not a list, okie dokie, nothing to convert
		pass

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


def ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=UseJavaDefault, addResEntropy=UseJavaDefault, energyPartition=UseJavaDefault, amat=UseJavaDefault, approximationErrorBudget=UseJavaDefault):
	'''
	:java:classdoc:`.energy.ConfEnergyCalculator`

	:builder_option confSpace .energy.ConfEnergyCalculator$Builder#confSpace:
	:builder_option ecalc .energy.ConfEnergyCalculator$Builder#ecalc:
	:builder_option referenceEnergies .energy.ConfEnergyCalculator$Builder#eref:
	:builder_option addResEntropy .energy.ConfEnergyCalculator$Builder#addResEntropy:
	:builder_option energyPartition .energy.ConfEnergyCalculator$Builder#epart:
	:builder_option amat .energy.ConfEnergyCalculator$Builder#amat:
	:builder_option approximationErrorBudget .energy.ConfEnergyCalculator$Builder#approximationErrorBudget:
	:builder_return .energy.ConfEnergyCalculator$Builder:
	'''
	builder = _get_builder(c.energy.ConfEnergyCalculator)(confSpace, ecalc)

	if referenceEnergies is not UseJavaDefault:
		builder.setReferenceEnergies(referenceEnergies)
	if addResEntropy is not UseJavaDefault:
		builder.addResEntropy(addResEntropy)
	if energyPartition is not UseJavaDefault:
		builder.setEnergyPartition(energyPartition)
	if amat is not UseJavaDefault:
		builder.setApproximatorMatrix(amat)
	if approximationErrorBudget is not UseJavaDefault:
		builder.setApproximationErrorBudget(approximationErrorBudget)

	return builder.build()


def ConfEnergyCalculatorCopy(source, ecalc):
	'''
	:java:classdoc:`.energy.ConfEnergyCalculator`

	:param source: The conformation energy calculator you wish to copy.
	:type source: :java:ref:`.energy.ConfEnergyCalculator`
	:builder_option ecalc .energy.ConfEnergyCalculator$Builder#ecalc:
	:builder_return .energy.ConfEnergyCalculator$Builder:
	'''
	return c.energy.ConfEnergyCalculator(source, ecalc)


def ApproximatorMatrix(confEcalc, cacheFile=UseJavaDefault, numSamplesPerParam=UseJavaDefault):
	'''
	:java:classdoc:`.energy.approximation.ApproximatorMatrix`

	:builder_option confEcalc .energy.approximation.ApproximatorMatrixCalculator#confEcalc:
	:builder_option cacheFile .energy.approximation.ApproximatorMatrixCalculator#cacheFile:
	:builder_option numSamplesPerParam .energy.approximation.ApproximatorMatrixCalculator#numSamplesPerParam:
	'''

	calculator = c.energy.approximation.ApproximatorMatrixCalculator(confEcalc)

	if cacheFile is not UseJavaDefault:
		calculator.setCacheFile(jvm.toFile(cacheFile))
	if numSamplesPerParam is not UseJavaDefault:
		calculator.setNumSamplesPerDof(numSamplesPerParam)

	return calculator.calc()


def EnergyMatrix(confEcalc, cacheFile=UseJavaDefault, tripleCorrectionThreshold=UseJavaDefault, quadCorrectionThreshold=UseJavaDefault):
	'''
	:java:methoddoc:`.ematrix.SimplerEnergyMatrixCalculator#calcEnergyMatrix`

	:builder_option confEcalc .ematrix.SimplerEnergyMatrixCalculator$Builder#confEcalc:
	:builder_option cacheFile .ematrix.SimplerEnergyMatrixCalculator$Builder#cacheFile:
	:builder_option tripleCorrectionThreshold .ematrix.SimplerEnergyMatrixCalculator$Builder#tripleCorrectionThreshold:
	:builder_option quadCorrectionThreshold .ematrix.SimplerEnergyMatrixCalculator$Builder#quadCorrectionThreshold:
	'''
	
	builder = _get_builder(c.ematrix.SimplerEnergyMatrixCalculator)(confEcalc)

	if cacheFile is not UseJavaDefault:
		builder.setCacheFile(jvm.toFile(cacheFile))
	if tripleCorrectionThreshold is not UseJavaDefault:
		builder.setTripleCorrectionThreshold(jvm.boxDouble(tripleCorrectionThreshold))
	if quadCorrectionThreshold is not UseJavaDefault:
		builder.setQuadCorrectionThreshold(jvm.boxDouble(quadCorrectionThreshold))

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


def DEE(confSpace, emat,
		singlesThreshold=useJavaDefault, pairsThreshold=useJavaDefault,
		singlesGoldsteinDiffThreshold=useJavaDefault, pairsGoldsteinDiffThreshold=useJavaDefault, triplesGoldsteinDiffThreshold=useJavaDefault,
		typeDependent=useJavaDefault, numIterations=useJavaDefault,
		singlesPlugThreshold=useJavaDefault, pairsPlugThreshold=useJavaDefault, triplesPlugThreshold=useJavaDefault,
		singlesTransitivePruning=useJavaDefault, pairsTransitivePruning=useJavaDefault, triplesTransitivePruning=useJavaDefault,
		showProgress=useJavaDefault, parallelism=useJavaDefault, cacheFile=useJavaDefault
	):
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
	:builder_option triplesGoldsteinDiffThreshold .pruning.SimpleDEE$Runner#triplesGoldsteinDiffThreshold:
	:builder_option singlesPlugThreshold .pruning.SimpleDEE$Runner#singlesPlugThreshold:
	:builder_option pairsPlugThreshold .pruning.SimpleDEE$Runner#pairsPlugThreshold:
	:builder_option triplesPlugThreshold .pruning.SimpleDEE$Runner#triplesPlugThreshold:
	:builder_option singlesTransitivePruning .pruning.SimpleDEE$Runner#singlesTransitivePruning:
	:builder_option pairsTransitivePruning .pruning.SimpleDEE$Runner#pairsTransitivePruning:
	:builder_option triplesTransitivePruning .pruning.SimpleDEE$Runner#triplesTransitivePruning:
	:builder_option typeDependent .pruning.SimpleDEE$Runner#typeDependent:
	:builder_option numIterations .pruning.SimpleDEE$Runner#numIterations:
	:builder_option showProgress .pruning.SimpleDEE$Runner#showProgress:
	'''

	runner = _get_builder(c.pruning.SimpleDEE, 'Runner')()

	if singlesThreshold is not useJavaDefault:
		runner.setSinglesThreshold(jvm.boxDouble(singlesThreshold))
	if pairsThreshold is not useJavaDefault:
		runner.setPairsThreshold(jvm.boxDouble(pairsThreshold))
	if singlesGoldsteinDiffThreshold is not useJavaDefault:
		runner.setSinglesGoldsteinDiffThreshold(jvm.boxDouble(singlesGoldsteinDiffThreshold))
	if pairsGoldsteinDiffThreshold is not useJavaDefault:
		runner.setPairsGoldsteinDiffThreshold(jvm.boxDouble(pairsGoldsteinDiffThreshold))
	if triplesGoldsteinDiffThreshold is not useJavaDefault:
		runner.setTriplesGoldsteinDiffThreshold(jvm.boxDouble(triplesGoldsteinDiffThreshold))
	if typeDependent is not useJavaDefault:
		runner.setTypeDependent(typeDependent)
	if numIterations is not useJavaDefault:
		runner.setNumIterations(numIterations)
	if singlesTransitivePruning is not useJavaDefault:
		runner.setSinglesTransitivePruning(singlesTransitivePruning)
	if pairsTransitivePruning is not useJavaDefault:
		runner.setPairsTransitivePruning(pairsTransitivePruning)
	if triplesTransitivePruning is not useJavaDefault:
		runner.setTriplesTransitivePruning(triplesTransitivePruning)
	if singlesPlugThreshold is not useJavaDefault:
		runner.setSinglesPlugThreshold(singlesPlugThreshold)
	if pairsPlugThreshold is not useJavaDefault:
		runner.setPairsPlugThreshold(pairsPlugThreshold)
	if triplesPlugThreshold is not useJavaDefault:
		runner.setTriplesPlugThreshold(triplesPlugThreshold)
	if showProgress is not useJavaDefault:
		runner.setShowProgress(showProgress)
	if parallelism is not useJavaDefault:
		runner.setParallelism(parallelism)
	if cacheFile is not useJavaDefault:
		runner.setCacheFile(jvm.toFile(cacheFile))

	return runner.run(confSpace, emat)


def DEE_read(confSpace, path):
	'''
	Reads a pruning matrix from a file

	:param confSpace: The design conformation space
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param str path: Path to the file
	'''

	return c.pruning.SimpleDEE.read(confSpace, jvm.toFile(path))


def AStarTraditional(emat, confSpaceOrPmat, showProgress=True, useExternalMemory=False, maxNumNodes=useJavaDefault):
	'''
	:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#setTraditional`

	:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:
	:param confSpaceOrPmat: The conformation space containing the residue conformations to search.
	:type confSpaceOrPmat: :java:ref:`.confspace.SimpleConfSpace` or :java:ref:`.pruning.PruningMatrix`
	:param useExternalMemory: set to True to use external memory.

		:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#useExternalMemory`

	:type useExternalMemory: boolean
	:builder_option maxNumNodes .astar.conf.ConfAStarTree$Builder#maxNumNodes:
	:builder_return .astar.conf.ConfAStarTree$Builder:
	'''
	builder = _get_builder(c.astar.conf.ConfAStarTree)(emat, confSpaceOrPmat)
	builder.setShowProgress(showProgress)

	# set node limits before picking heuristics
	if maxNumNodes is not useJavaDefault:
		builder.setMaxNumNodes(jvm.boxLong(maxNumNodes))

	builder.setTraditional()

	if useExternalMemory == True:
		builder.useExternalMemory()

	return builder.build()


def EdgeUpdater():
	return c.astar.conf.scoring.mplp.EdgeUpdater()

def NodeUpdater():
	return c.astar.conf.scoring.mplp.NodeUpdater()

def AStarMPLP(emat, confSpaceOrPmat, updater=None, numIterations=None, convergenceThreshold=None, useExternalMemory=False, maxNumNodes=useJavaDefault):
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
	:builder_option maxNumNodes .astar.conf.ConfAStarTree$Builder#maxNumNodes:
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

	# set node limits before picking heuristics
	if maxNumNodes is not useJavaDefault:
		builder.setMaxNumNodes(jvm.boxLong(maxNumNodes))

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

def Paste(complexConfSpace, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, useWindowCriterion=useJavaDefault, maxNumPfConfs=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None, useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault, mutFile=useJavaDefault):
    '''
    :java:classdoc:`.paste.Paste`

    For examples using PAStE, see the examples/python.Paste directory in your Osprey distribution.

	:param complexConfSpace: :java:fielddoc:`.paste.Paste#protein`
	:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:builder_option epsilon .paste.Paste$Settings$Builder#epsilon:
	:builder_option stabilityThreshold .paste.Paste$Settings$Builder#stabilityThreshold:
	:builder_option maxSimultaneousMutations .paste.Paste$Settings$Builder#maxSimultaneousMutations:
	:builder_option maxNumPfConfs .paste.Paste$Settings$Builder#maxNumPfConfs:
	:builder_option useExternalMemory .paste.Paste$Settings$Builder#useExternalMemory:
	:builder_option showPfuncProgress .paste.Paste$Settings$Builder#showPfuncProgress:
	:builder_option useWindowCriterion .paste.Paste$Settings$Builder#useWindowCriterion:
	:param str addMutFile: Path to the file that has the mutant sequences of interest
	:param bool writeSequencesToConsole: True to write sequences and scores to the console
	:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

	:rtype: :java:ref:`.paste.Paste`
	'''

    # build settings
    settingsBuilder = _get_builder(jvm.getInnerClass(c.paste.Paste, 'Settings'))()
    if useWindowCriterion is not useJavaDefault:
        settingsBuilder.setUseWindowCriterion(useWindowCriterion)
    if mutFile is not useJavaDefault:
        settingsBuilder.addMutFile(jvm.toFile(mutFile))
    if epsilon is not useJavaDefault:
        settingsBuilder.setEpsilon(epsilon)
	if stabilityThreshold is not useJavaDefault:
		settingsBuilder.setStabilityThreshold(jvm.boxDouble(stabilityThreshold))
	if maxNumPfConfs is not useJavaDefault:
	    settingsBuilder.setPfConfs(maxNumPfConfs)
	if maxSimultaneousMutations is not useJavaDefault:
		settingsBuilder.setMaxSimultaneousMutations(maxSimultaneousMutations)
	if writeSequencesToConsole:
		settingsBuilder.addScoreConsoleWriter()
	if writeSequencesToFile is not None:
		settingsBuilder.addScoreFileWriter(jvm.toFile(writeSequencesToFile))
	if useExternalMemory is not useJavaDefault:
		settingsBuilder.setExternalMemory(useExternalMemory)
	if showPfuncProgress is not useJavaDefault:
		settingsBuilder.setShowPfuncProgress(showPfuncProgress)
	settings = settingsBuilder.build()

	return c.paste.Paste(complexConfSpace, settings)

def _PasteConfSearchFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.getInnerClass(c.paste.Paste, 'ConfSearchFactory'),
		dict={ 'make': func }
	)

Paste.ConfSearchFactory = _PasteConfSearchFactory


def KStar(proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None, useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault):
	'''
	:java:classdoc:`.kstar.KStar`

	For examples using K*, see the examples/python.KStar directory in your Osprey distribution.

	:param proteinConfSpace: :java:fielddoc:`.kstar.KStar#protein`
	:type proteinConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param ligandConfSpace: :java:fielddoc:`.kstar.KStar#ligand`
	:type ligandConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param complexConfSpace: :java:fielddoc:`.kstar.KStar#complex`
	:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:builder_option epsilon .kstar.KStar$Settings$Builder#epsilon:
	:builder_option stabilityThreshold .kstar.KStar$Settings$Builder#stabilityThreshold:
	:builder_option maxSimultaneousMutations .kstar.KStar$Settings$Builder#maxSimultaneousMutations:
	:builder_option useExternalMemory .kstar.KStar$Settings$Builder#useExternalMemory:
	:builder_option showPfuncProgress .kstar.KStar$Settings$Builder#showPfuncProgress:
	:param bool writeSequencesToConsole: True to write sequences and scores to the console
	:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

	:rtype: :java:ref:`.kstar.KStar`
	'''

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
	if useExternalMemory is not useJavaDefault:
		settingsBuilder.setExternalMemory(useExternalMemory)
	if showPfuncProgress is not useJavaDefault:
		settingsBuilder.setShowPfuncProgress(showPfuncProgress)
	settings = settingsBuilder.build()

	return c.kstar.KStar(proteinConfSpace, ligandConfSpace, complexConfSpace, settings)


def _KStarConfSearchFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.getInnerClass(c.kstar.KStar, 'ConfSearchFactory'),
		dict={ 'make': func }
	)

KStar.ConfSearchFactory = _KStarConfSearchFactory


def BBKStar(proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, energyMatrixCachePattern=useJavaDefault, useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault, numBestSequences=useJavaDefault, numConfsPerBatch=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None):
	'''
	:java:classdoc:`.kstar.BBKStar`

	For examples using BBK*, see the examples/python.KStar directory in your Osprey distribution.

	:param proteinConfSpace: :java:fielddoc:`.kstar.BBKStar#protein`
	:type proteinConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param ligandConfSpace: :java:fielddoc:`.kstar.BBKStar#ligand`
	:type ligandConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param complexConfSpace: :java:fielddoc:`.kstar.BBKStar#complex`
	:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
	:builder_option epsilon .kstar.KStar$Settings$Builder#epsilon:
	:builder_option stabilityThreshold .kstar.KStar$Settings$Builder#stabilityThreshold:
	:builder_option maxSimultaneousMutations .kstar.KStar$Settings$Builder#maxSimultaneousMutations:
	:builder_option useExternalMemory .kstar.KStar$Settings$Builder#useExternalMemory:
	:builder_option showPfuncProgress .kstar.KStar$Settings$Builder#showPfuncProgress:
	:builder_option numBestSequences .kstar.BBKStar$Settings$Builder#numBestSequences:
	:builder_option numConfsPerBatch .kstar.BBKStar$Settings$Builder#numConfsPerBatch:
	:param bool writeSequencesToConsole: True to write sequences and scores to the console
	:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

	:rtype: :java:ref:`.kstar.BBKStar`
	'''

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
	if useExternalMemory is not useJavaDefault:
		kstarSettingsBuilder.setExternalMemory(useExternalMemory)
	if showPfuncProgress is not useJavaDefault:
		kstarSettingsBuilder.setShowPfuncProgress(showPfuncProgress)
	kstarSettings = kstarSettingsBuilder.build()

	bbkstarSettingsBuilder = _get_builder(jvm.getInnerClass(c.kstar.BBKStar, 'Settings'))()
	if numBestSequences is not useJavaDefault:
		bbkstarSettingsBuilder.setNumBestSequences(numBestSequences)
	if numConfsPerBatch is not useJavaDefault:
		bbkstarSettingsBuilder.setNumConfsPerBatch(numConfsPerBatch)
	bbkstarSettings = bbkstarSettingsBuilder.build()

	return c.kstar.BBKStar(proteinConfSpace, ligandConfSpace, complexConfSpace, kstarSettings, bbkstarSettings)


def ConfAnalyzer(confEcalc):
	'''
	:java:classdoc:`.gmec.ConfAnalyzer`

	For examples using the conf analyzer, see examples/python.GMEC/analyzeConf.py in your Osprey distribution.

	:param confEcalc: :java:fielddoc:`.gmec.SimpleGMECFinder$Builder#confEcalc`
	:type confEcalc: :java:ref:`.energy.ConfEnergyCalculator`

	:rtype: :java:ref:`.gmec.ConfAnalyzer`
	'''

	return c.gmec.ConfAnalyzer(confEcalc)


def SequenceAnalyzer(kstar):
	'''
	:java:classdoc:`.kstar.SequenceAnalyzer`

	For examples using the sequence analyzer, see examples/python.KStar/analyzeSequence.py in your Osprey distribution.

	:param kstar: a configured instance of KStar
	:type kstar: :java:ref:`.kstar.KStar`

	:rtype: :java:ref:`.kstar.SequenceAnalyzer`
	'''

	return c.kstar.SequenceAnalyzer(kstar)


class LUTE_SamplingStrategy:
	Progressive = 0
	PairsAndTriples = 1

def LUTE_train(confEcalc, emat, pmat, maxRMSE=0.1, maxOverfittingScore=1.5, randomSeed=12345, confDBPath=None, samplingStrategy=LUTE_SamplingStrategy.Progressive):
	'''
	Trains a LUTE model

	For examples using LUTE, see examples/python.GMEC/LUTE.*.py and examples/python.KStar/LUTE.*.py in your Osprey distribution.

	:param confEcalc: The conformation energy calculator
	:type confEcalc: :java:ref:`.energy.ConfEnergyCalculator`
	:param emat: An energy matrix
	:type emat: :java:ref:`.ematrix.EnergyMatrix`
	:param pmat: A pruning matrix, resulting from DEE
	:type pmat: :java:ref:`.pruning.PruningMatrix`

	:param float maxRMSE: The maximum tolerable fit RMS error
	:param float maxOverfittingScore: The maximum tolerable amount of overfitting (score = training set RMSE / test set RMSE)
	:param int randomSeed: Random seed to use for conformation sampling
	:param str confDBPath: Path to write/read confDB file, or None to omit saving the confDB to disk

	:returns: The LUTE model if it meets the accuracy goals, or None otherwise
	:rtype: :java:ref:`.lute.LUTEState`
	'''

	confSpace = confEcalc.confSpace

	# make a conf DB, saved to a file if needed
	if confDBPath is None:
		confDB = c.confspace.ConfDB(confSpace)
	else:
		confDB = c.confspace.ConfDB(confSpace, jvm.toFile(confDBPath))

	try:

		# make a conf table for LUTE
		confTable = jvm.getInnerClass(c.confspace.ConfDB, 'ConfTable')(confDB, 'LUTE')

		# use the OLSCG fitter for LUTE (it's a little faster than LASSO in practice)
		fitter = jvm.getInnerClass(c.lute.LUTE, 'Fitter').OLSCG

		# train LUTE
		lute = c.lute.LUTE(confSpace)
		sampler = c.lute.RandomizedDFSConfSampler(confSpace, pmat, randomSeed)
		if samplingStrategy == LUTE_SamplingStrategy.PairsAndTriples:
			sampleAndFitFunc = lute.sampleAllPairsTriplesAndFit
		else:
			sampleAndFitFunc = lute.sampleTuplesAndFit
		fitGoodEnough = sampleAndFitFunc(confEcalc, emat, pmat, confTable, sampler, fitter, maxOverfittingScore, maxRMSE)

		lute.reportConfSpaceSize(pmat)

		# if the fit wasn't good enough, don't send the trained model back
		if not fitGoodEnough:
			return None

		# otherwise, return the trained LUTE model
		return c.lute.LUTEState(lute.getTrainingSystem())

	finally:
		confDB.close()


def LUTE_write(model, path):
	'''
	Writes a LUTE model to a file

	:param model: The LUTE model
	:type model: :java:ref:`.lute.LUTEState`
	:param str path: Path to the file
	'''

	file = jvm.toFile(path)
	c.lute.LUTEIO.write(model, file)
	print('LUTE model saved to %s' % file.getAbsolutePath())


def LUTE_read(path):
	'''
	Reads a LUTE model from a file

	:param str path: Path to the file

	:returns: The LUTE model
	:rtype: :java:ref:`.lute.LUTEState`
	'''

	file = jvm.toFile(path)
	model = c.lute.LUTEIO.read(file)
	print('LUTE model read from %s' % file.getAbsolutePath())
	return model


def LUTE_ConfEnergyCalculator(confSpace, model):
	'''
	Creates a LUTE conformation energy calculator

	:param confSpace: The conformation space
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param model: The LUTE model
	:type model: :java:ref:`.lute.LUTEState`

	:rtype: :java:ref:`.lute.LUTEConfEnergyCalculator`
	'''

	return c.lute.LUTEConfEnergyCalculator(confSpace, model)


def LUTE_AStar(rcs, pmat, luteEcalc, showProgress=True):
	'''
	:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#setLUTE`

	:builder_option rcs .astar.conf.ConfAStarTree$Builder#rcs:
	:param pmat: The pruning matrix from the LUTE training calculation.
	:type pmat: :java:ref:`.pruning.PruningMatrix`
	:param luteEcalc: The LUTE conformation energy calculator
	:type luteEcalc: :java:ref:`.lute.LUTEConfEnergyCalculator`

	:builder_return .astar.conf.ConfAStarTree$Builder:
	'''

	# filter the rcs by the pmat
	rcs = c.astar.conf.RCs(rcs, pmat)

	builder = _get_builder(c.astar.conf.ConfAStarTree)(None, rcs)
	builder.setShowProgress(showProgress)
	builder.setLUTE(luteEcalc)

	return builder.build()


def LUTE_GMECFinder(confSpace, model, pmat, confLog=useJavaDefault, printIntermediateConfs=useJavaDefault):
	'''
	:java:classdoc:`.lute.LUTEGMECFinder`

	:param confSpace: The conformation space
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
	:param model: The LUTE model
	:type model: :java:ref:`.lute.LUTEState`
	:param pmat: The pruning matrix from the LUTE training calculation.
	:type pmat: :java:ref:`.pruning.PruningMatrix`
	:param str confLog: Path to file where conformations found during conformation space search should be logged.
	:builder_option printIntermediateConfs .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole:

	:rtype: :java:ref:`.lute.LUTEGMECFinder`
	'''

	builder = _get_builder(c.lute.LUTEGMECFinder)(pmat, LUTE_ConfEnergyCalculator(confSpace, model))

	if confLog is not useJavaDefault:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not useJavaDefault:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	return builder.build()


def COMETS_State(name, confSpace):
	'''
	:java:classdoc:`.gmec.Comets$State`

	:param str name: :java:fielddoc:`.gmec.Comets$State#name`
	:param confSpace: :java:fielddoc:`.gmec.Comets$State#confSpace`
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

	:rtype: :java:ref:`.gmec.Comets$State`
	'''

	return jvm.getInnerClass(c.gmec.Comets, 'State')(name, confSpace)


def COMETS_ConfSearchFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.c.java.util.function.Function,
		dict={ 'apply': func }
	)


def COMETS_LME(weightsByState, offset=useJavaDefault, constrainLessThan=None):
	'''
	:java:classdoc:`.gmec.Comets$LME`

	:param weightsByState: map from states to weights
	:type weightsByState: map from :java:ref:`.gmec.Comets$State` to float

	:builder_option offset .gmec.Comets$LME$Builder#offset:
	:param float constrainLessThan: :java:methoddoc:`.gmec.Comets$LME$Builder#constrainLessThan`
	:builder_return .gmec.Comets$LME$Builder:
	'''

	builder = _get_builder(jvm.getInnerClass(c.gmec.Comets, 'LME'))()

	if offset is not useJavaDefault:
		builder.setOffset(offset)

	for (state, weight) in weightsByState.items():
		builder.addState(state, weight)

	if constrainLessThan is not None:
		builder.constrainLessThan(constrainLessThan)

	return builder.build()


def COMETS(objective, constraints=[], objectiveWindowSize=useJavaDefault, objectiveWindowMax=useJavaDefault, maxSimultaneousMutations=useJavaDefault, minNumConfTrees=useJavaDefault, logFile=None):
	'''
	:java:classdoc:`.gmec.Comets`

	:builder_option objective .gmec.Comets$Builder#objective:

	:param constraints: List of LMEs to use as constraints
	:type constraints: list of :java:ref:`.gmec.Comets$LME`

	:builder_option objectiveWindowSize .gmec.Comets$Builder#objectiveWindowSize:
	:builder_option objectiveWindowMax .gmec.Comets$Builder#objectiveWindowMax:
	:builder_option maxSimultaneousMutations .gmec.Comets$Builder#maxSimultaneousMutations:
	:builder_option minNumConfTrees .gmec.Comets$Builder#minNumConfTrees:

	:param str logFile: :java:fielddoc:`.gmec.Comets$Builder#logFile`

	:builder_return .gmec.Comets$Builder:
	'''

	builder = _get_builder(c.gmec.Comets)(objective)

	for constraint in constraints:
		builder.addConstraint(constraint)

	if objectiveWindowSize is not useJavaDefault:
		builder.setObjectiveWindowSize(objectiveWindowSize)
	if objectiveWindowMax is not useJavaDefault:
		builder.setObjectiveWindowMax(objectiveWindowMax)
	if maxSimultaneousMutations is not useJavaDefault:
		builder.setMaxSimultaneousMutations(maxSimultaneousMutations)
	if minNumConfTrees is not useJavaDefault:
		builder.setMinNumConfTrees(jvm.boxInt(minNumConfTrees))

	if logFile is not None:
		builder.setLogFile(jvm.toFile(logFile))

	return builder.build()


def MSKStar_State(name, confSpace):
	'''
	:java:classdoc:`.kstar.MSKStar$State`

	:param str name: :java:fielddoc:`.kstar.MSKStar$State#name`
	:param confSpace: :java:fielddoc:`.kstar.MSKStar$State#confSpace`
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

	:rtype: :java:ref:`.kstar.MSKStar$State`
	'''

	return jvm.getInnerClass(c.kstar.MSKStar, 'State')(name, confSpace)


def MSKStar_ConfSearchFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.c.java.util.function.Function,
		dict={ 'apply': func }
	)


def MSKStar_LMFE(weightsByState, offset=useJavaDefault, constrainLessThan=None):
	'''
	:java:classdoc:`.kstar.MSKStar$LMFE`

	:param weightsByState: map from states to weights
	:type weightsByState: map from :java:ref:`.kstar.MSKStar$State` to float

	:builder_option offset .kstar.MSKStar$LMFE$Builder#offset:
	:param float constrainLessThan: :java:methoddoc:`.kstar.MSKStar$LMFE$Builder#constrainLessThan`
	:builder_return .kstar.MSKStar$LMFE$Builder:
	'''

	builder = _get_builder(jvm.getInnerClass(c.kstar.MSKStar, 'LMFE'))()

	if offset is not useJavaDefault:
		builder.setOffset(offset)

	for (state, weight) in weightsByState.items():
		builder.addState(state, weight)

	if constrainLessThan is not None:
		builder.constrainLessThan(constrainLessThan)

	return builder.build()


def MSKStar(objective, constraints=[], epsilon=useJavaDefault, objectiveWindowSize=useJavaDefault, objectiveWindowMax=useJavaDefault, maxSimultaneousMutations=useJavaDefault, minNumConfTrees=useJavaDefault, logFile=None):
	'''
	:java:classdoc:`.kstar.MSKStar`

	:builder_option objective .kstar.MSKStar$Builder#objective:

	:param constraints: List of LMFEs to use as constraints
	:type constraints: list of :java:ref:`.kstar.MSKStar$LMFE`

	:builder_option epsilon .kstar.MSKStar$Builder#epsilon:
	:builder_option objectiveWindowSize .kstar.MSKStar$Builder#objectiveWindowSize:
	:builder_option objectiveWindowMax .kstar.MSKStar$Builder#objectiveWindowMax:
	:builder_option maxSimultaneousMutations .kstar.MSKStar$Builder#maxSimultaneousMutations:
	:builder_option minNumConfTrees .kstar.MSKStar$Builder#minNumConfTrees:

	:param str logFile: :java:fielddoc:`.kstar.MSKStar$Builder#logFile`

	:builder_return .kstar.MSKStar$Builder:
	'''

	builder = _get_builder(c.kstar.MSKStar)(objective)

	for constraint in constraints:
		builder.addConstraint(constraint)

	if objectiveWindowSize is not useJavaDefault:
		builder.setObjectiveWindowSize(objectiveWindowSize)
	if objectiveWindowMax is not useJavaDefault:
		builder.setObjectiveWindowMax(objectiveWindowMax)
	if maxSimultaneousMutations is not useJavaDefault:
		builder.setMaxSimultaneousMutations(maxSimultaneousMutations)
	if minNumConfTrees is not useJavaDefault:
		builder.setMinNumConfTrees(jvm.boxInt(minNumConfTrees))

	if logFile is not None:
		builder.setLogFile(jvm.toFile(logFile))

	return builder.build()

def PartitionFunctionFactory(confSpace, confEcalc, state, confUpperBoundcalc=None):
    pfuncFactory = c.kstar.pfunc.PartitionFunctionFactory(confSpace, confEcalc, state)
    if confUpperBoundcalc is not None:
        pfuncFactory.setUseMARKStar(confUpperBoundcalc)
    else:
        pfuncFactory.setUseGradientDescent()

    return pfuncFactory

def EwakstarDoer_ConfSearchFactory(func):

    # convert the python lambda to a JVM interface implementation
    return jpype.JProxy(
        jvm.c.java.util.function.Function,
        dict={ 'apply': func }
	)

def EwakstarDoer_State(name, confSpace):

	return jvm.getInnerClass(c.ewakstar.EwakstarDoer, 'State')(name, confSpace)

def EwakstarDoer(state, smaNodes, useSMA=useJavaDefault, printPDBs=useJavaDefault, useWtBenchmark=useJavaDefault, numEWAKStarSeqs=useJavaDefault, logFile=None, epsilon=useJavaDefault, pfEw=useJavaDefault, eW=useJavaDefault, orderOfMag=useJavaDefault, numPfConfs=useJavaDefault, numTopSeqs=useJavaDefault, mutableType=useJavaDefault, numMutable=useJavaDefault, seqFilterOnly=useJavaDefault, numCPUs=useJavaDefault):

    builder = _get_builder(c.ewakstar.EwakstarDoer)()

    builder.addState(state)

    if useSMA is not useJavaDefault:
        builder.setupSMA(useSMA, smaNodes)
    if useWtBenchmark is not useJavaDefault:
        builder.setUseWtBenchmark(useWtBenchmark)
    if printPDBs is not useJavaDefault:
        builder.setPrintPDBs(printPDBs)
    if epsilon is not useJavaDefault:
        builder.setEpsilon(epsilon)
    if pfEw is not useJavaDefault:
        builder.setPfEw(pfEw)
    if eW is not useJavaDefault:
        builder.setEw(eW)
    if orderOfMag is not useJavaDefault:
        builder.setOrderOfMag(orderOfMag)
    if numPfConfs is not useJavaDefault:
        builder.setNumPfConfs(numPfConfs)
    if numTopSeqs is not useJavaDefault:
        builder.setNumTopOverallSeqs(numTopSeqs)
    if mutableType is not useJavaDefault:
        builder.setMutableType(mutableType)
    if numMutable is not useJavaDefault:
        builder.setNumMutable(numMutable)
    if seqFilterOnly is not useJavaDefault:
        builder.setSeqFilterOnly(seqFilterOnly)
    if numCPUs is not useJavaDefault:
        builder.setNumCpus(numCPUs)
    if numEWAKStarSeqs is not useJavaDefault:
        builder.setNumEWAKStarSeqs(numEWAKStarSeqs)

    if logFile is not None:
        builder.setLogFile(jvm.toFile(logFile))

    return builder.build()



def SOFEA_StateConfig(emat, confEcalc, confdbPath=None):
	'''
	:java:classdoc:`.sofea.Sofea$StateConfig`

	:param emat: :java:fielddoc:`.sofea.Sofea$StateConfig#emat`
	:param confEcalc: :java:fielddoc:`.sofea.Sofea$StateConfig#confEcalc`
	:param confdbPath: :java:fielddoc:`.sofea.Sofea$StateConfig#confDBFile`
	:rtype: :java:ref:`.sofea.Sofea$StateConfig`
	'''
	confdbFile = jvm.toFile(confdbPath) if confdbPath is not None else None
	return jvm.getInnerClass(c.sofea.Sofea, 'StateConfig')(emat, confEcalc, confdbFile)


def SOFEA(confSpace, configFunc, seqdbPath='sofea.seqdb', seqdbMathContext=useJavaDefault, fringedbLowerPath='sofea.lower.fringedb', fringedbLowerMiB=10, fringedbUpperPath='sofea.upper.fringedb', fringedbUpperMiB=10, showProgress=useJavaDefault, sweepIncrement=useJavaDefault, maxNumMinimizations=useJavaDefault, negligableFreeEnergy=useJavaDefault):
	'''
	:java:classdoc:`.sofea.Sofea`

	:param confSpace: A multi-state configuration space
	:type confSpace: :java:ref:`.confspace.MultiStateConfSpace`

	:param configFunc: a function that creates a :java:ref:`.sofea.Sofea$StateConfig` for a state
	:type configFunc: function(:java:ref:`.confspace.MultiStateConfSpace$State`) returning :java:ref:`.sofea.Sofea$StateConfig`

	:param str seqdbPath: Path to write the sequence database file
	:builder_option seqdbMathContext .sofea.Sofea$Builder#seqdbMathContext:
	:param str fringedbLowerPath: Path to write the lower fringe set
	:param int fringedbLowerMiB: size of the lower fringe set in MiB
	:param str fringedbUpperPath: Path to write the upper fringe set
	:param int fringedbUpperMiB: size of the upper fringe set in MiB
	:builder_option showProgress .sofea.Sofea$Builder#showProgress:
	:builder_option sweepIncrement .sofea.Sofea$Builder#sweepIncrement:
	:builder_option maxNumMinimizations .sofea.Sofea$Builder#maxNumMinimizations:
	:builder_option negligableFreeEnergy .sofea.Sofea$Builder#negligableFreeEnergy:

	:builder_return .sofea.Sofea$Builder:
	'''

	builder = _get_builder(c.sofea.Sofea)(confSpace)

	for state in confSpace.states:
		builder.configState(state, configFunc(state))

	if seqdbPath is not useJavaDefault:
		builder.setSeqDBFile(jvm.toFile(seqdbPath))
	if seqdbMathContext is not useJavaDefault:
		builder.setSeqDBMathContext(seqdbMathContext)
	if fringedbLowerPath is not useJavaDefault:
		builder.setFringeDBLowerFile(jvm.toFile(fringedbLowerPath))
	if fringedbLowerMiB is not useJavaDefault:
		builder.setFringeDBLowerMiB(fringedbLowerMiB)
	if fringedbUpperPath is not useJavaDefault:
		builder.setFringeDBUpperFile(jvm.toFile(fringedbUpperPath))
	if fringedbUpperMiB is not useJavaDefault:
		builder.setFringeDBUpperMiB(fringedbUpperMiB)
	if showProgress is not useJavaDefault:
		builder.setShowProgress(showProgress)
	if sweepIncrement is not useJavaDefault:
		builder.setSweepIncrement(sweepIncrement)
	if maxNumMinimizations is not useJavaDefault:
		builder.setMaxNumMinimizations(maxNumMinimizations)
	if negligableFreeEnergy is not useJavaDefault:
		builder.setNegligableFreeEnergy(negligableFreeEnergy)

	return builder.build()


def SOFEA_MinLMFE(lmfe, numSequences, minFreeEnergyWidth):
	'''
	:java:classdoc:`.sofea.MinLMFE`

	:param lmfe: :java:fielddoc:`.sofea.MinLMFE#objective`
	:param numSequences: :java:fielddoc:`.sofea.MinLMFE#numSequences`
	:param minFreeEnergyWidth: :java:fielddoc:`.sofea.MinLMFE#minFreeEnergyWidth`
	:rtype: :java:ref:`.sofea.MinLMFE`
	'''
	return c.sofea.MinLMFE(lmfe, numSequences, minFreeEnergyWidth)


def SOFEA_SequenceLMFE(sequence, lmfe, minFreeEnergyWidth):
	'''
	:java:classdoc:`.sofea.SequenceLMFE`

	:param sequence: :java:fielddoc:`.sofea.SequenceLMFE#seq`
	:param lmfe: :java:fielddoc:`.sofea.SequenceLMFE#lmfe`
	:param minFreeEnergyWidth: :java:fielddoc:`.sofea.SequenceLMFE#minFreeEnergyWidth`
	:rtype: :java:ref:`.sofea.SequenceLMFE`
	'''
	return c.sofea.SequenceLMFE(sequence, lmfe, minFreeEnergyWidth)
