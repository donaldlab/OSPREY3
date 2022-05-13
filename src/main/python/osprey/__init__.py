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

from __future__ import absolute_import

import sys, os, platform, jpype, traceback

from . import jvm
from . import wraps


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
TemplateMatchingMethod = None

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
		ex = value.__javaobject__
		ex.printStackTrace()

		# NOTE: in JPype-py2, information about python Exceptions
		# are not embedded within the Java exceptions from things like JProxy =(
		# so we couldn't print it here even if we wanted to.
		# But the Python3 version of JPype prints that info by default! =)

		# print causes too, if any
		ex = ex.getCause()
		while ex is not None:
			print('Caused by:')
			ex.printStackTrace()
			ex = ex.getCause()

	except (AttributeError, TypeError):
		# must not be a java exception
		pass


# should we print the preamble?
_print_preamble = True
try:
	_print_preamble = os.environ['OSPREY_PREAMBLE'].lower() != 'false'
except:
	pass


def _start_jvm_common(fn_jvm_start):

	# disable buffered output on stdout, so python log messages line up with java log messages
	class Unbuffered(object):
		def __init__(self, stream):
			self.stream = stream

		def write(self, data):
			self.stream.write(data)
			self.stream.flush()

		def writelines(self, datas):
			self.stream.writelines(datas)
			self.stream.flush()

		def __getattr__(self, attr):
			return getattr(self.stream, attr)

	sys.stdout = Unbuffered(sys.stdout)

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
		osprey_dir = os.path.join(os.path.dirname(__file__), '../../../../')
		classpath_path = os.path.join(osprey_dir, 'build/python/classpath.txt')
		if not os.path.isfile(classpath_path):
			raise Exception('dev classpath for python not generated yet. run ./gradlew pythonDevelop')

		for path in open(classpath_path, 'r').readlines():
			jvm.addClasspath(path.strip())

	# build the JRE path
	if no_escape_from_reality:

		# release environment: use the bundled JRE
		if sys.platform in ('win32', 'cygwin'):
			jre_path = os.path.join(os.path.dirname(__file__), 'jre', 'bin', 'server', 'jvm.dll')
		elif sys.platform == 'darwin':
			jre_path = os.path.join(os.path.dirname(__file__), 'jre', 'lib', 'server', 'libjvm.dylib')
		else:
			jre_path = os.path.join(os.path.dirname(__file__), 'jre', 'lib', 'server', 'libjvm.so')

	else:

		# development environment: use the system JRE
		jre_path = None

	fn_jvm_start(jre_path)

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
	global TemplateMatchingMethod
	TemplateMatchingMethod = jvm.getInnerClass(c.structure.Residue, 'TemplateMatchingMethod')

	# expose static builder methods too
	Parallelism.makeCpu = c.parallelism.Parallelism.makeCpu
	Parallelism.make = c.parallelism.Parallelism.make

	# print the preamble, if needed
	if _print_preamble:
		osprey_version = c.control.Main.Version
		python_version = '.'.join([str(x) for x in sys.version_info[0:3]])
		java_version = jpype.java.lang.System.getProperty('java.version')
		print("OSPREY %s, Python %s, Java %s, %s" % (osprey_version, python_version, java_version, platform.platform()))


def start_with_jvm_args(jvm_args=None):
	'''
	Start the jvm with JVM arguments.

	# Arguments
	jvm_args `[str]`: a list of arguments to the JVM. For example, `["-XX:MaxHeapSize=100g", "-XX:MinHeapSize=100g"]`.
	defaults to None, which implies the JVM will run with its default configuration.
	'''
	if jvm_args is None:
		jvm_args = []

	_start_jvm_common(lambda jre_path: jvm.start_with_args(jre_path, jvm_args))


def start(heapSizeMiB=1024, enableAssertions=False, stackSizeMiB=16, garbageSizeMiB=128, allowRemoteManagement=False, attachJvmDebugger=False):
	'''
	Starts the Java Virtual Machine (JVM) that runs Osprey's computation libraries.

	Call #start before using any of Osprey's other functions.

	# Arguments
	heapSizeMiB `int`: Size of the JVM heap in megabytes. This is essentially the amount of memory
		Osprey will have to do computations. 1024 MiB is 1 GiB, but for larger designs,
		you may want to use 2048 MiB (2 GiB), 4096 MiB (4 GiB), or even more memory.

	enableAssertions `bool`: pass `True` to enable JVM assertions. Only useful for debugging.

	stackSizeMiB `int`: Size of the JVM stack portion of the heap in megabytes.
		Generally leave this at the default value, unless Osprey crashes because it's too small. Then try increasing it.

	garbageSizeMiB `int`: Size of the garbage portion of the JVM heap that is reserved for temporary objects.
		This default value is appropriate for the default heap size, but if using larger heap sizes, then increasing
		the garbage size to 256, 512, or even 1024 MiB can give a modest improvement in performance.

	allowRemoteManagement `bool`: pass `True` to listen on ports 9010 and 9011 for JMX remote management.

	attachJvmDebugger `bool`: pass `True` to be able to attach a Java debugger to the JVM.
	'''
	_start_jvm_common(lambda jre_path: jvm.start(jre_path, heapSizeMiB, enableAssertions, stackSizeMiB, garbageSizeMiB, allowRemoteManagement, attachJvmDebugger))

	# print the preamble, if needed
	if _print_preamble:

		print("Using up to %d MiB heap memory: %d MiB for garbage, %d MiB for storage" % (
			heapSizeMiB, garbageSizeMiB, heapSizeMiB - garbageSizeMiB
		))


def readTextFile(path):
	'''
	Useful for passing file contents to functions that expect strings rather than file paths

	# Arguments
	path `str`: path to the text file

	# Returns
	`str`: the text of the file at `path`
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

	NOTE: Raw molecules cannot be used directly in designs.
	Create a #Strand with the molecule to use it in a design.

	# Arguments
	path `str`: path to PDB file

	# Returns
	${returns_method_java(.structure.PDBIO#readFile(String)Molecule)}
	'''
	mol = c.structure.PDBIO.readFile(path)
	print('read PDB file from file: %s' % path)
	return mol


def writePdb(mol, path, comment=None):
	'''
	save a molecule to a PDB file

	# Arguments
	mol ${type_java(.structure.Molecule)}: the molecule to save
	path `str`: path of the PDB file
	comment `str`: optional comment to add to the PDB file headers
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

	# Arguments
	internalSizeMiB `int`: ${method_javadoc(.externalMemory.ExternalMemory#setInternalLimit)}
	tempDir `str`: ${method_javadoc(.externalMemory.ExternalMemory#setTempDir(String)void)}
	               ${arg_javadoc(.externalMemory.ExternalMemory#setTempDir(String)void, dir)}
	               ${default(tempDir, <system temp dir>)}
	tempSubDir `str`: ${arg_javadoc(.externalMemory.ExternalMemory#setTempDir(String,String)void, subdir)}
	                  ${default(tempSubdir, <automatically generated>)}
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

def Parallelism(cpuCores=useJavaDefault, gpus=useJavaDefault, streamsPerGpu=useJavaDefault):
	'''
	${class_javadoc(.parallelism.Parallelism)}

	# Arguments
	${args_fields_javadoc(.parallelism.Parallelism$Builder,
		[cpuCores, numCpus],
		[gpus, numGpus],
		[streamsPerGpu, numStreamsPerGpu]
	)}

	# Returns
	${returns_method_java(.parallelism.Parallelism$Builder#build)}
	'''
	builder = _get_builder(c.parallelism.Parallelism)()
	if cpuCores is not useJavaDefault:
		builder.setNumCpus(cpuCores)
	if gpus is not useJavaDefault:
		builder.setNumGpus(gpus)
	if streamsPerGpu is not useJavaDefault:
		builder.setNumStreamsPerGpu(streamsPerGpu)

	return builder.build()


def TemplateLibrary(
		ffparams=None,
		defaultTemplates=True, extraTemplates=[],
		defaultTemplateCoords=True, extraTemplateCoords=[],
		defaultRotamers=True, extraRotamers=[],
		extraBackboneDependentRotamers=[],
		defaultResidueEntropies=True, extraResidueEntropies=[],
		makeDAminoAcids=useJavaDefault,
		moleculesForWildTypeRotamers=[]
	):
	'''
	${class_javadoc(.restypes.ResidueTemplateLibrary)}

	# Arguments
	${arg_field_javadoc(ffparams, .restypes.ResidueTemplateLibrary$Builder#ffparams)}

	defaultTemplates `bool`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#clearTemplates)}
	extraTemplates `[str]`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#addTemplates)}
	                        Each element can be the text of a template file, or the path to one.

	defaultTemplateCoords `bool`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#clearTemplateCoords)}
	extraTemplateCoords `[str]`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#addTemplateCoords)}
	                             Each element can be the text of a template coords file, or the path to one.

	defaultRotamers `bool`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#clearRotamers)}
	extraRotamers `[str]`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#addRotamers)}
	                       Each element can be the text of a rotamers file, or the path to one.

	extraBackboneDependentRotamers `[str]`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#addBackboneDependentRotamers)}
	                                        Each element can be the text of a backbone-dependent rotamers file, or the path to one.
	${arg_field_javadoc(makeDAminoAcids, .restypes.ResidueTemplateLibrary$Builder#makeDAminoAcidTemplates)}

	moleculesForWildTypeRotamers `[`${type_java(.structure.Molecule)}`]`: ${method_javadoc(.restypes.ResidueTemplateLibrary$Builder#addMoleculeForWildTypeRotamers)}

	# Returns
	${returns_method_java(.restypes.ResidueTemplateLibrary$Builder#build)}
	'''

	if ffparams is None:
		builder = _get_builder(c.restypes.ResidueTemplateLibrary)()
	else:
		builder = _get_builder(c.restypes.ResidueTemplateLibrary)(ffparams)

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

	if makeDAminoAcids is not useJavaDefault:
		builder.setMakeDAminoAcidTemplates(makeDAminoAcids)

	for mol in moleculesForWildTypeRotamers:
		builder.addMoleculeForWildTypeRotamers(mol)

	return builder.build()
	

def Strand(pathOrMol, residues=None, templateLib=useJavaDefault, templateMatchingMethod=useJavaDefault):
	'''
	${class_javadoc(.confspace.Strand)}

	# Arguments
	pathOrMol `str` or ${type_java(.structure.Molecule)}: path to a PDB file, or a molecule instance
	residues `[str, str]`: range of residue numbers, inclusive. `None` to include all residues.
	${args_fields_javadoc(.confspace.Strand$Builder,
		[templateLib],
		[templateMatchingMethod]
	)}

	# Returns
	${returns_method_java(.confspace.Strand$Builder#build)}
	'''

	# did we get a path or a molecule?
	if isinstance(pathOrMol, c.structure.Molecule):
		mol = pathOrMol
	else:
		mol = readPdb(pathOrMol)

	builder = _get_builder(c.confspace.Strand)(mol)

	if residues is not None:
		builder.setResidues(residues[0], residues[1])

	if templateLib is not useJavaDefault:
		builder.setTemplateLibrary(templateLib)

	if templateMatchingMethod is not useJavaDefault:
		builder.setTemplateMatchingMethod(templateMatchingMethod)

	return builder.build()


def ConfSpace(strands, shellDist=useJavaDefault):
	'''
	${class_javadoc(.confspace.SimpleConfSpace)}

	# Arguments
	strands ${type_java(.confspace.Strand)} or [${type_java(.confspace.Strand)}]: the strands to use
	${arg_field_javadoc(shellDist, .confspace.SimpleConfSpace$Builder#shellDist)}

	# Returns
	${returns_method_java(.confspace.SimpleConfSpace$Builder#build)}
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
	if shellDist is not useJavaDefault:
		builder.setShellDistance(shellDist)
	
	return builder.build()


class StateMutable(object):
	def __init__(self, name, confSpace):
		'''
		# Arguments
		name `str`: a unique name for the state
		confSpace ${type_java(.confspace.SimpleConfSpace)}: the conformation space
		'''
		self.name = name
		self.confSpace = confSpace

class StateUnmutable(object):
	def __init__(self, name, confSpace):
		'''
		# Arguments
		name `str`: a unique name for the state
		confSpace ${type_java(.confspace.SimpleConfSpace)}: the conformation space
		'''
		self.name = name
		self.confSpace = confSpace


def MultiStateConfSpace(states):
	'''
	${class_javadoc(.confspace.MultiStateConfSpace)}

	# Arguments
	states [`#StateMutable` or `#StateUnmutable`]: a list of states to use

	# Returns
	${returns_method_java(.confspace.MultiStateConfSpace$Builder#build)}
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


def ForcefieldParams(params=None):
	'''
	${class_javadoc(.energy.forcefield.ForcefieldParams)}
	
	Configure the forcefield parameters by setting the properties of the ${type_java(.energy.forcefield.ForcefieldParams)} object.

	# Arguments
	params `str`: text of forcefield parameters in Amber format, to override defaults if desired

	# Returns
	${type_java(.energy.forcefield.ForcefieldParams)}
	'''
	if params:
		return c.energy.forcefield.ForcefieldParams(Forcefield.AMBER, params)
	else:
		return c.energy.forcefield.ForcefieldParams(Forcefield.AMBER)


def EnergyCalculator(confSpace, ffparams, parallelism=useJavaDefault, type=useJavaDefault, isMinimizing=useJavaDefault, infiniteWellEnergy=useJavaDefault):
	'''
	${class_javadoc(.energy.EnergyCalculator)}

	# Arguments
	confSpace ${type_java(.confspace.SimpleConfSpace)}: The conformation space containing the residue templates to use for atom connectivities.
	${args_fields_javadoc(.energy.EnergyCalculator$Builder,
		[ffparams],
		[parallelism],
		[type],
		[isMinimizing],
		[infiniteWellEnergy, type=float]
	)}

	# Returns
	${returns_method_java(.energy.EnergyCalculator$Builder#build)}
	'''

	# convert confSpace to a jvm list if possible
	try:
		confSpace = jvm.toArrayList(confSpace)
	except TypeError:
		# not a list, okie dokie, nothing to convert
		pass

	builder = _get_builder(c.energy.EnergyCalculator)(confSpace, ffparams)

	if parallelism is not useJavaDefault:
		builder.setParallelism(parallelism)

	if type is not useJavaDefault:
		builder.setType(type)

	if isMinimizing is not useJavaDefault:
		builder.setIsMinimizing(isMinimizing)

	if infiniteWellEnergy is not useJavaDefault:
		builder.setInfiniteWellEnergy(jvm.boxDouble(infiniteWellEnergy))

	return builder.build()


def SharedEnergyCalculator(ecalc, isMinimizing=useJavaDefault):
	'''
	${class_javadoc(.energy.EnergyCalculator$SharedBuilder)}

	# Arguments
	${args_fields_javadoc(.energy.EnergyCalculator$SharedBuilder,
		[ecalc, parent],
		[isMinimizing]
	)}
	
	# Returns
	${returns_method_java(.energy.EnergyCalculator$SharedBuilder#build)}
	'''
	builder = jvm.getInnerClass(c.energy.EnergyCalculator, 'SharedBuilder')(ecalc)

	if isMinimizing is not useJavaDefault:
		builder.setIsMinimizing(isMinimizing)

	return builder.build()


def ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=useJavaDefault, addResEntropy=useJavaDefault, energyPartition=useJavaDefault, amat=useJavaDefault, approximationErrorBudget=useJavaDefault):
	'''
	${class_javadoc(.energy.ConfEnergyCalculator)}

	# Arguments
	${args_fields_javadoc(.energy.ConfEnergyCalculator$Builder,
		[confSpace],
		[ecalc],
		[referenceEnergies, eref],
		[addResEntropy],
		[energyPartition, epart],
		[amat],
		[approximationErrorBudget]
	)}

	# Returns
	${returns_method_java(.energy.ConfEnergyCalculator$Builder#build)}
	'''

	builder = _get_builder(c.energy.ConfEnergyCalculator)(confSpace, ecalc)

	if referenceEnergies is not useJavaDefault:
		builder.setReferenceEnergies(referenceEnergies)
	if addResEntropy is not useJavaDefault:
		builder.addResEntropy(addResEntropy)
	if energyPartition is not useJavaDefault:
		builder.setEnergyPartition(energyPartition)
	if amat is not useJavaDefault:
		builder.setApproximatorMatrix(amat)
	if approximationErrorBudget is not useJavaDefault:
		builder.setApproximationErrorBudget(approximationErrorBudget)

	return builder.build()


def ConfEnergyCalculatorCopy(source, ecalc):
	'''
	Makes a copy of a ${type_java(.energy.ConfEnergyCalculator)} that uses a different underlying ${type_java(.energy.EnergyCalculator)}

	# Arguments
	source ${type_java(.energy.ConfEnergyCalculator)}: The conformation energy calculator you wish to copy.
	ecalc ${type_java(.energy.EnergyCalculator)}: The energy calculator to re-use.

	# Returns
	${returns_method_java(.energy.ConfEnergyCalculator$Builder#build)}
	'''
	return c.energy.ConfEnergyCalculator(source, ecalc)


def ApproximatorMatrix(confEcalc, cacheFile=useJavaDefault, numSamplesPerParam=useJavaDefault):
	'''
	${class_javadoc(.energy.approximation.ApproximatorMatrix)}

	# Arguments
	${args_fields_javadoc(.energy.approximation.ApproximatorMatrixCalculator,
		[confEcalc],
		[cacheFile, type=str],
		[numSamplesPerParam]
	)}

	# Returns
	${returns_method_java(.energy.approximation.ApproximatorMatrixCalculator#calc()ApproximatorMatrix)}
	'''

	calculator = c.energy.approximation.ApproximatorMatrixCalculator(confEcalc)

	if cacheFile is not useJavaDefault:
		calculator.setCacheFile(jvm.toFile(cacheFile))
	if numSamplesPerParam is not useJavaDefault:
		calculator.setNumSamplesPerDof(numSamplesPerParam)

	return calculator.calc()


def EnergyMatrix(confEcalc, cacheFile=useJavaDefault, tripleCorrectionThreshold=useJavaDefault, quadCorrectionThreshold=useJavaDefault):
	'''
	${method_javadoc(.ematrix.SimplerEnergyMatrixCalculator#calcEnergyMatrix)}

	# Arguments
	${args_fields_javadoc(.ematrix.SimplerEnergyMatrixCalculator$Builder,
		[confEcalc],
		[cacheFile, type=str],
		[tripleCorrectionThreshold, type=float],
		[quadCorrectionThreshold, type=float]
	)}

	# Returns
	${returns_method_java(.ematrix.SimplerEnergyMatrixCalculator#calcEnergyMatrix)}
	'''
	
	builder = _get_builder(c.ematrix.SimplerEnergyMatrixCalculator)(confEcalc)

	if cacheFile is not useJavaDefault:
		builder.setCacheFile(jvm.toFile(cacheFile))
	if tripleCorrectionThreshold is not useJavaDefault:
		builder.setTripleCorrectionThreshold(jvm.boxDouble(tripleCorrectionThreshold))
	if quadCorrectionThreshold is not useJavaDefault:
		builder.setQuadCorrectionThreshold(jvm.boxDouble(quadCorrectionThreshold))

	return builder.build().calcEnergyMatrix()


def ReferenceEnergies(confSpace, ecalc, addResEntropy=useJavaDefault):
	'''
	${method_javadoc(.ematrix.SimplerEnergyMatrixCalculator#calcReferenceEnergies)}

	# Arguments
	${args_fields_javadoc(.ematrix.SimpleReferenceEnergies$Builder,
		[confSpace],
		[ecalc],
		[addResEntropy]
	)}

	# Returns
	${returns_method_java(.ematrix.SimplerEnergyMatrixCalculator#calcReferenceEnergies)}
	'''

	builder = _get_builder(c.ematrix.SimpleReferenceEnergies)(confSpace, ecalc)

	if addResEntropy is not useJavaDefault:
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
	${class_javadoc(.pruning.SimpleDEE$Runner)}

	# Arguments
	confSpace ${type_java(.confspace.SimpleConfSpace)}: The design conformation space
	emat ${type_java(.ematrix.EnergyMatrix)}: An energy matrix computed for the conformation space
	${args_fields_javadoc(.pruning.SimpleDEE$Runner,
		[singlesThreshold, type=float],
		[pairsThreshold, type=float],
		[singlesGoldsteinDiffThreshold, type=float],
		[pairsGoldsteinDiffThreshold, type=float],
		[triplesGoldsteinDiffThreshold, type=float],
		[singlesPlugThreshold, type=float],
		[pairsPlugThreshold, type=float],
		[triplesPlugThreshold, type=float],
		[singlesTransitivePruning],
		[pairsTransitivePruning],
		[triplesTransitivePruning],
		[typeDependent],
		[numIterations],
		[showProgress],
		[parallelism],
		[cacheFile, type=str]
	)}

	# Returns
	${returns_method_java(.pruning.SimpleDEE$Runner#run)}
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

	# Arguments
	confSpace ${type_java(.confspace.SimpleConfSpace)}: The design conformation space
	path `str`: Path to the file

	# Returns
	${returns_method_java(.pruning.SimpleDEE#read)}
	'''

	return c.pruning.SimpleDEE.read(confSpace, jvm.toFile(path))


def AStarTraditional(emat, confSpaceOrPmat, showProgress=True, useExternalMemory=False, maxNumNodes=useJavaDefault):
	'''
	${method_javadoc(.astar.conf.ConfAStarTree$Builder#setTraditional)}

	# Arguments
	${arg_field_javadoc(emat, .astar.conf.ConfAStarTree$Builder#emat)}
	confSpaceOrPmat ${type_java(.confspace.SimpleConfSpace)} or ${type_java(.pruning.PruningMatrix)}:
	                The conformation space containing the residue conformations to search,
	                or the pruning matrix.
	useExternalMemory `bool`: set to `True` to use external memory.
	                          ${method_javadoc(.astar.conf.ConfAStarTree$Builder#useExternalMemory)}
	${arg_field_javadoc(maxNumNodes, .astar.conf.ConfAStarTree$Builder#maxNumNodes, type=int)}

	# Returns
	${returns_method_java(.astar.conf.ConfAStarTree$Builder#build)}
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

def AStarMPLP(emat, confSpaceOrPmat, updater=useJavaDefault, numIterations=useJavaDefault, convergenceThreshold=useJavaDefault, useExternalMemory=False, maxNumNodes=useJavaDefault):
	'''
	${method_javadoc(.astar.conf.ConfAStarTree$Builder#setMPLP()Builder)}

	# Arguments
	${arg_field_javadoc(emat, .astar.conf.ConfAStarTree$Builder#emat)}
	confSpaceOrPmat ${type_java(.confspace.SimpleConfSpace)} or ${type_java(.pruning.PruningMatrix)}:
	                The conformation space containing the residue conformations to search,
	                or the pruning matrix.
	${args_fields_javadoc(.astar.conf.ConfAStarTree$MPLPBuilder,
		[updater],
		[numIterations],
		[convergenceThreshold]
	)}
	useExternalMemory `bool`: set to `True` to use external memory.
	                          ${method_javadoc(.astar.conf.ConfAStarTree$Builder#useExternalMemory)}
	${arg_field_javadoc(maxNumNodes, .astar.conf.ConfAStarTree$Builder#maxNumNodes, type=int)}

	# Returns
	${returns_method_java(.astar.conf.ConfAStarTree$Builder#build)}
	'''
	mplpBuilder = _get_builder(c.astar.conf.ConfAStarTree, 'MPLPBuilder')()

	if updater is not useJavaDefault:
		mplpBuilder.setUpdater(updater)

	if numIterations is not useJavaDefault:
		mplpBuilder.setNumIterations(numIterations)

	if convergenceThreshold is not useJavaDefault:
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


def GMECFinder(astar, confEcalc, confLog=None, printIntermediateConfs=useJavaDefault, useExternalMemory=False, confDBFile=useJavaDefault):
	'''
	${class_javadoc(.gmec.SimpleGMECFinder)}

	# Arguments
	${arg_field_javadoc(astar, .gmec.SimpleGMECFinder$Builder#search)}
	   Use one of #AStarTraditional or #AStarMPLP to get an A* implementation.
	${arg_field_javadoc(confEcalc, .gmec.SimpleGMECFinder$Builder#confEcalc)}
	   Use #ConfEnergyCalculator to get a conformation energy calculator.
	confLog `str`: Path to file where conformations found during conformation space search should be logged.
	${args_fields_javadoc(.gmec.SimpleGMECFinder$Builder,
		[printIntermediateConfs, printIntermediateConfsToConsole],
		[useExternalMemory],
		[confDBFile, confDB, type=str]
	)}

	# Returns
	${returns_method_java(.gmec.SimpleGMECFinder$Builder#build)}
	'''

	builder = _get_builder(c.gmec.SimpleGMECFinder)(astar, confEcalc)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not useJavaDefault:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	if useExternalMemory == True:
		builder.useExternalMemory()

	if confDBFile is not useJavaDefault:
		builder.setConfDB(jvm.toFile(confDBFile))

	return builder.build()


def DEEGMECFinder(emat, confSpace, ecalc, confEcalc, name, use_epic, use_lute, confLog=None, printIntermediateConfs=useJavaDefault, useExternalMemory=False):
	'''
	${class_javadoc(.gmec.DEEGMECFinder)}

	# Arguments
	${args_fields_javadoc(.gmec.DEEGMECFinder$Builder,
		[emat],
		[confSpace],
		[ecalc]
	)}
	${arg_field_javadoc(confEcalc, .gmec.SimpleGMECFinder$Builder#confEcalc)}
	   Use #ConfEnergyCalculator to get a conformation energy calculator.
	${arg_field_javadoc(name, .gmec.DEEGMECFinder$Builder#name)}
	use_epic `bool`: `True` to use EPIC
	use_lute `bool`: `True` to use LUTE
	confLog `str`: Path to file where conformations found during conformation space search should be logged.
	${args_fields_javadoc(.gmec.SimpleGMECFinder$Builder,
		[printIntermediateConfs, printIntermediateConfsToConsole],
		[useExternalMemory]
	)}

	# Returns
	${returns_method_java(.gmec.DEEGMECFinder$Builder#build)}
	'''

	builder = _get_builder(c.gmec.DEEGMECFinder)(emat, confSpace, ecalc, confEcalc, name)

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not useJavaDefault:
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
	'''
	${class_javadoc(.confspace.DEEPerStrandFlex)}

	# Arguments
	${arg_java(strand, .confspace.DEEPerStrandFlex#<init>)}
	pert_file_name `str`: Path to the perturbation file
	flex_res_list `[str]`: Flexible residues
	${arg_java(pdb_file, .dof.deeper.DEEPerSettings#<init>(boolean,String,boolean,String,boolean,double,double,boolean,ArrayList<String>,String,boolean,ResidueTemplateLibrary), PDBFile)}

	# Returns
	${type_java(.confspace.DEEPerStrandFlex)}
	'''
	deeper_settings = c.dof.deeper.DEEPerSettings(True, pert_file_name, True, 'None', False, 2.5, 2.5, False, jvm.toArrayList(flex_res_list), pdb_file, False, strand.templateLib)
	bbflex = c.confspace.DEEPerStrandFlex(strand,deeper_settings)
	return bbflex


def Paste(
	complexConfSpace, numPDBs, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault,
	useWindowCriterion=useJavaDefault, maxNumPfConfs=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None,
	useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault, mutFile=useJavaDefault
):
	'''
	${class_javadoc(.paste.Paste)}

	For examples using PAStE, see the `examples/python.Paste` directory in your Osprey distribution.

	# Arguments
	${arg_java(complexConfSpace, .paste.Paste#<init>, protein)}
	${arg_java(numPDBs, .paste.Paste#<init>)}
	${args_fields_javadoc(.paste.Paste$Settings$Builder,
		[epsilon],
		[stabilityThreshold, type=float],
		[maxSimultaneousMutations],
		[useWindowCriterion]
		[maxNumPfConfs]
	)}
	writeSequencesToConsole `bool`: True to write sequences and scores to the console
	writeSequencesToFile `str`: Path to the log file to write sequences scores (in TSV format), or None to skip logging
	${args_fields_javadoc(.paste.Paste$Settings$Builder,
		[useExternalMemory],
		[showPfuncProgress],
		[mutFile, type=str]
	)}

	# Returns
	${type_java(.paste.Paste)}
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

	return c.paste.Paste(complexConfSpace, settings, numPDBs)

def _PasteConfSearchFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.getInnerClass(c.paste.Paste, 'ConfSearchFactory'),
		dict={ 'make': func }
	)

Paste.ConfSearchFactory = _PasteConfSearchFactory


def KStar(
	proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault,
	maxSimultaneousMutations=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None,
	useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault
):
	'''
	${class_javadoc(.kstar.KStar)}

	For examples using K*, see the examples/python.KStar directory in your Osprey distribution.

	# Arguments
	${args_fields_javadoc(.kstar.KStar,
		[proteinConfSpace, protein],
		[ligandConfSpace, ligand],
		[complexConfSpace, complex]
	)}
	${args_fields_javadoc(.kstar.KStar$Settings$Builder,
		[epsilon],
		[stabilityThreshold, type=float],
		[maxSimultaneousMutations],
		[useExternalMemory],
		[showPfuncProgress]
	)}
	writeSequencesToConsole `bool`: True to write sequences and scores to the console
	writeSequencesToFile `str`: Path to the log file to write sequences scores (in TSV format), or None to skip logging

	# Returns
	${type_java(.kstar.KStar)}
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


def _KStarPfuncFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.getInnerClass(c.kstar.KStar, 'PfuncFactory'),
		dict={ 'make': func }
	)

KStar.PfuncFactory = _KStarPfuncFactory


def PartitionFunction(confEcalc, confSearchUpper, confSearchLower, rcs):
	'''
	${class_javadoc(.kstar.pfunc.GradientDescentPfunc)}

	# Arguments
	${args_java(.kstar.pfunc.GradientDescentPfunc#<init>(ConfEnergyCalculator,ConfSearch,ConfSearch,BigInteger),
		[confEcalc, ecalc],
		[confSearchUpper, upperBoundConfs],
		[confSearchLower, lowerBoundConfs]
	)}
	rcs ${type_java(.astar.conf.RCs)}:

	# Returns
	${type_java(.kstar.pfunc.GradientDescentPfunc)}
	'''

	return c.kstar.pfunc.GradientDescentPfunc(
		confEcalc,
		confSearchUpper,
		confSearchLower,
		rcs.getNumConformations()
	)


def MARKStarPfunc(confSpace, ematMinimized, confEcalcMinimized, ematRigid, confEcalcRigid, rcs):
	'''
	${class_javadoc(.markstar.framework.MARKStarBoundFastQueues)}

	# Arguments
	${args_java(.markstar.framework.MARKStarBoundFastQueues#<init>,
		[confSpace],
		[ematMinimized, minimizingEmat],
		[confEcalcMinimized, minimizingConfEcalc],
		[ematRigid, rigidEmat]
	)}
	confEcalcRigid ${type_java(.energy.ConfEnergyCalculator)}: ignored for now, reserved for future use
	${args_java(.markstar.framework.MARKStarBoundFastQueues#<init>,
		[rcs]
	)}

	# Returns
	${type_java(.markstar.framework.MARKStarBoundFastQueues)}
	'''

	pfunc = c.markstar.framework.MARKStarBoundFastQueues(
		confSpace,
		ematRigid,
		ematMinimized,
		confEcalcMinimized,
		rcs,
		confEcalcMinimized.ecalc.parallelism
	)
	pfunc.setCorrections(c.ematrix.UpdatingEnergyMatrix(confSpace, ematMinimized, confEcalcMinimized))
	return pfunc


def BBKStar(
	proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault,
	maxSimultaneousMutations=useJavaDefault, useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault,
	numBestSequences=useJavaDefault, numConfsPerBatch=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None
):
	'''
	${class_javadoc(.kstar.BBKStar)}

	For examples using BBK*, see the examples/python.KStar directory in your Osprey distribution.

	# Arguments
	${args_fields_javadoc(.kstar.BBKStar,
		[proteinConfSpace, protein],
		[ligandConfSpace, ligand],
		[complexConfSpace, complex]
	)}
	${args_fields_javadoc(.kstar.KStar$Settings$Builder,
		[epsilon],
		[stabilityThreshold, type=float],
		[maxSimultaneousMutations],
		[useExternalMemory],
		[showPfuncProgress]
	)}
	${args_fields_javadoc(.kstar.BBKStar$Settings$Builder,
		[numBestSequences],
		[numConfsPerBatch]
	)}
	writeSequencesToConsole `bool`: True to write sequences and scores to the console
	writeSequencesToFile `str`: Path to the log file to write sequences scores (in TSV format), or None to skip logging

	# Returns
	${type_java(.kstar.BBKStar)}
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


def _BBKStarConfSearchFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.getInnerClass(c.kstar.BBKStar, 'ConfSearchFactory'),
		dict={ 'make': func }
	)

BBKStar.ConfSearchFactory = _BBKStarConfSearchFactory


def ConfAnalyzer(confEcalc):
	'''
	${class_javadoc(.gmec.ConfAnalyzer)}

	For examples using the conf analyzer, see examples/python.GMEC/analyzeConf.py in your Osprey distribution.

	# Arguments
	${arg_java(confEcalc, .gmec.ConfAnalyzer#<init>)}

	# Returns
	${type_java(.gmec.ConfAnalyzer)}
	'''

	return c.gmec.ConfAnalyzer(confEcalc)


def SequenceAnalyzer(kstar):
	'''
	${class_javadoc(.kstar.SequenceAnalyzer)}

	For examples using the sequence analyzer, see examples/python.KStar/analyzeSequence.py in your Osprey distribution.

	# Arguments
	${arg_java(kstar, .kstar.SequenceAnalyzer#<init>(KStar))} a configured instance of KStar

	# Returns
	${type_java(.kstar.SequenceAnalyzer)}
	'''

	return c.kstar.SequenceAnalyzer(kstar)


class LUTE_SamplingStrategy:
	Progressive = 0
	PairsAndTriples = 1

def LUTE_train(confEcalc, emat, pmat, maxRMSE=0.1, maxOverfittingScore=1.5, randomSeed=12345, confDBPath=None, samplingStrategy=LUTE_SamplingStrategy.Progressive):
	'''
	Trains a LUTE model

	For examples using LUTE, see examples/python.GMEC/LUTE.*.py and examples/python.KStar/LUTE.*.py in your Osprey distribution.

	# Arguments
	confEcalc ${type_java(.energy.ConfEnergyCalculator)}: The conformation energy calculator
	emat: ${type_java(.ematrix.EnergyMatrix)}: An energy matrix
	pmat: ${type_java(.pruning.PruningMatrix)}: A pruning matrix, resulting from DEE
	maxRMSE `float`: The maximum tolerable fit RMS error
	maxOverfittingScore `float`: The maximum tolerable amount of overfitting (score = training set RMSE / test set RMSE)
	randomSeed `int`: Random seed to use for conformation sampling
	confDBPath `str`: Path to write/read confDB file, or None to omit saving the confDB to disk

	# Returns
	${type_java(.lute.LUTEState)}
	The LUTE model if it meets the accuracy goals, or None otherwise
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

	# Arguments
	model ${type_java(.lute.LUTEState)}: The LUTE model
	path `str`: Path to the file
	'''

	file = jvm.toFile(path)
	c.lute.LUTEIO.write(model, file)
	print('LUTE model saved to %s' % file.getAbsolutePath())


def LUTE_read(path):
	'''
	Reads a LUTE model from a file

	# Arguments
	path `str`: Path to the file

	# Returns
	${type_java(.lute.LUTEState)}
	The LUTE model
	'''

	file = jvm.toFile(path)
	model = c.lute.LUTEIO.read(file)
	print('LUTE model read from %s' % file.getAbsolutePath())
	return model


def LUTE_ConfEnergyCalculator(confSpace, model):
	'''
	Creates a LUTE conformation energy calculator

	# Arguments
	${args_java(.lute.LUTEConfEnergyCalculator#<init>,
		[confSpace],
		[model, state]
	)}

	# Returns
	${type_java(.lute.LUTEConfEnergyCalculator)}
	'''

	return c.lute.LUTEConfEnergyCalculator(confSpace, model)


def LUTE_AStar(rcs, pmat, luteEcalc, showProgress=True):
	'''
	${method_javadoc(.astar.conf.ConfAStarTree$Builder#setLUTE)}

	# Arguments
	${arg_java(rcs, .astar.conf.RCs#<init>(RCs,PruningMatrix), other)}
	${arg_java(pmat, .astar.conf.RCs#<init>(RCs,PruningMatrix))} The pruning matrix from the LUTE training calculation.
	luteEcalc ${type_java(.lute.LUTEConfEnergyCalculator)}: The LUTE conformation energy calculator
	showProgress `bool`: True to show progress of the A* calculation

	# Returns
	${returns_method_java(.astar.conf.ConfAStarTree$Builder#build)}
	'''

	# filter the rcs by the pmat
	rcs = c.astar.conf.RCs(rcs, pmat)

	builder = _get_builder(c.astar.conf.ConfAStarTree)(None, rcs)
	builder.setShowProgress(showProgress)
	builder.setLUTE(luteEcalc)

	return builder.build()


def LUTE_Pfunc(luteEcalc, astar, rcs):
	'''
	${class_javadoc(.lute.LUTEPfunc)}

	# Arguments
	${args_java(.lute.LUTEPfunc#<init>,
		[luteEcalc, ecalc],
		[astar]
	)}
	rcs ${type_java(.astar.conf.RCs)}
	'''

	return c.lute.LUTEPfunc(luteEcalc, astar, rcs.getNumConformations())


def LUTE_GMECFinder(confSpace, model, pmat, confLog=None, printIntermediateConfs=useJavaDefault):
	'''
	${class_javadoc(.lute.LUTEGMECFinder)}

	# Arguments
	confSpace ${type_java(.confspace.SimpleConfSpace)}: The conformation space
	model ${type_java(.lute.LUTEState)}: The LUTE model
	pmat ${type_java(.pruning.PruningMatrix)}: The pruning matrix from the LUTE training calculation.
	confLog `str`: Path to file where conformations found during conformation space search should be logged.
	${arg_field_javadoc(printIntermediateConfs, .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole)}

	# Returns
	${type_java(.lute.LUTEGMECFinder)}
	'''

	builder = _get_builder(c.lute.LUTEGMECFinder)(pmat, LUTE_ConfEnergyCalculator(confSpace, model))

	if confLog is not None:
		logFile = jvm.toFile(confLog)
		builder.setLogPrinter(c.gmec.LoggingConfPrinter(logFile))

	if printIntermediateConfs is not useJavaDefault:
		builder.setPrintIntermediateConfsToConsole(printIntermediateConfs)

	return builder.build()


def COMETS_State(name, confSpace):
	'''
	${class_javadoc(.gmec.Comets$State)}

	# Arguments
	${args_java(.gmec.Comets$State#<init>,
		[name],
		[confSpace]
	)}

	# Returns
	${type_java(.gmec.Comets$State)}
	'''

	return jvm.getInnerClass(c.gmec.Comets, 'State')(name, confSpace)


def COMETS_ConfSearchFactory(func):
	'''
	A utility to convert function types.

	See `examples/python.GMEC/comets.py` for an example of how to use this function.
	'''

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.c.java.util.function.Function,
		dict={ 'apply': func }
	)


def COMETS_LME(weightsByState, offset=useJavaDefault, constrainLessThan=None):
	'''
	${class_javadoc(.gmec.Comets$LME)}

	# Arguments
	weightsByState `{`${type_java(.gmec.Comets$State)}`-> float}`: map from states to weights
	${arg_field_javadoc(offset, .gmec.Comets$LME$Builder#offset)}
	constrainLessThan `float`: ${method_javadoc(.gmec.Comets$LME$Builder#constrainLessThan)}

	# Returns
	${returns_method_java(.gmec.Comets$LME$Builder#build)}
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
	${class_javadoc(.gmec.Comets)}

	# Arguments
	${arg_field_javadoc(objective, .gmec.Comets$Builder#objective)}
	constraints `[`${type_java(.gmec.Comets$LME)}`]`: List of LMEs to use as constraints
	${args_fields_javadoc(.gmec.Comets$Builder,
		[objectiveWindowSize],
		[objectiveWindowMax],
		[maxSimultaneousMutations],
		[minNumConfTrees, type=int]
	)}
	logFile `str`: ${field_javadoc(.gmec.Comets$Builder#logFile)}

	# Returns
	${returns_method_java(.gmec.Comets$Builder#build)}
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
# MSK* implementation is unfinished! don't include this in the public API
# 	'''
# 	${class_javadoc(.kstar.MSKStar$State)}
#
# 	:param str name: :java:fielddoc:`.kstar.MSKStar$State#name`
# 	:param confSpace: :java:fielddoc:`.kstar.MSKStar$State#confSpace`
# 	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
#
# 	:rtype: :java:ref:`.kstar.MSKStar$State`
# 	'''

	return jvm.getInnerClass(c.kstar.MSKStar, 'State')(name, confSpace)


def MSKStar_ConfSearchFactory(func):

	# convert the python lambda to a JVM interface implementation
	return jpype.JProxy(
		jvm.c.java.util.function.Function,
		dict={ 'apply': func }
	)


def MSKStar_LMFE(weightsByState, offset=useJavaDefault, constrainLessThan=None):
# MSK* implementation is unfinished! don't include this in the public API
# 	'''
# 	${class_javadoc(.kstar.MSKStar$LMFE)}
#
# 	:param weightsByState: map from states to weights
# 	:type weightsByState: map from :java:ref:`.kstar.MSKStar$State` to float
#
# 	:builder_option offset .kstar.MSKStar$LMFE$Builder#offset:
# 	:param float constrainLessThan: ${method_javadoc(.kstar.MSKStar$LMFE$Builder#constrainLessThan)}
# 	:builder_return .kstar.MSKStar$LMFE$Builder:
# 	'''

	builder = _get_builder(jvm.getInnerClass(c.kstar.MSKStar, 'LMFE'))()

	if offset is not useJavaDefault:
		builder.setOffset(offset)

	for (state, weight) in weightsByState.items():
		builder.addState(state, weight)

	if constrainLessThan is not None:
		builder.constrainLessThan(constrainLessThan)

	return builder.build()


def MSKStar(objective, constraints=[], epsilon=useJavaDefault, objectiveWindowSize=useJavaDefault, objectiveWindowMax=useJavaDefault, maxSimultaneousMutations=useJavaDefault, minNumConfTrees=useJavaDefault, logFile=None):
# MSK* implementation is unfinished! don't include this in the public API
# 	'''
# 	${class_javadoc(.kstar.MSKStar)}
#
# 	:builder_option objective .kstar.MSKStar$Builder#objective:
#
# 	:param constraints: List of LMFEs to use as constraints
# 	:type constraints: list of :java:ref:`.kstar.MSKStar$LMFE`
#
# 	:builder_option epsilon .kstar.MSKStar$Builder#epsilon:
# 	:builder_option objectiveWindowSize .kstar.MSKStar$Builder#objectiveWindowSize:
# 	:builder_option objectiveWindowMax .kstar.MSKStar$Builder#objectiveWindowMax:
# 	:builder_option maxSimultaneousMutations .kstar.MSKStar$Builder#maxSimultaneousMutations:
# 	:builder_option minNumConfTrees .kstar.MSKStar$Builder#minNumConfTrees:
#
# 	:param str logFile: :java:fielddoc:`.kstar.MSKStar$Builder#logFile`
#
# 	:builder_return .kstar.MSKStar$Builder:
# 	'''

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
	${class_javadoc(.sofea.Sofea$StateConfig)}

	# Arguments
	${args_java(.sofea.Sofea$StateConfig#<init>,
		[emat],
		[confEcalc],
		[confdbPath, confDBFile, type=str]
	)}

	# Returns
	${type_java(.sofea.Sofea$StateConfig)}
	'''
	confdbFile = jvm.toFile(confdbPath) if confdbPath is not None else None
	return jvm.getInnerClass(c.sofea.Sofea, 'StateConfig')(emat, confEcalc, confdbFile)


def SOFEA(
	confSpace, configFunc, seqdbPath='sofea.seqdb', seqdbMathContext=useJavaDefault, fringedbLowerPath='sofea.lower.fringedb',
	fringedbLowerMiB=10, fringedbUpperPath='sofea.upper.fringedb', fringedbUpperMiB=10, rcdbPath=useJavaDefault,
	showProgress=useJavaDefault, performanceLogPath=useJavaDefault, sweepIncrement=useJavaDefault,
	maxNumMinimizations=useJavaDefault, negligableFreeEnergy=useJavaDefault
):
	'''
	${class_javadoc(.sofea.Sofea)}

	See `examples/python.SOFEA/sofea.py` for an example of how to use SOFEA.

	# Arguments
	confSpace ${type_java(.confspace.MultiStateConfSpace)}: A multi-state configuration space
	configFunc: a function that accepts a ${type_java(.confspace.MultiStateConfSpace$State)}
	            and returns a ${type_java(.sofea.Sofea$StateConfig)} for a state
	seqdbPath `str`: Path to write the sequence database file
	${arg_field_javadoc(seqdbMathContext, .sofea.Sofea$Builder#seqdbMathContext)}
	fringedbLowerPath `str`: Path to write the lower fringe set
	fringedbLowerMiB `int`: size of the lower fringe set in MiB
	fringedbUpperPath `str`: Path to write the upper fringe set
	fringedbUpperMiB `int`: size of the upper fringe set in MiB
	${args_fields_javadoc(.sofea.Sofea$Builder,
		[rcdbPath, rcdbFile, type=str],
		[showProgress],
		[performanceLogPath, performanceLogFile, type=str],
		[sweepIncrement],
		[maxNumMinimizations],
		[negligableFreeEnergy]
	)}

	# Returns
	${returns_method_java(.sofea.Sofea$Builder#build)}
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
	if rcdbPath is not useJavaDefault:
		builder.setRCDBFile(jvm.toFile(rcdbPath))
	if showProgress is not useJavaDefault:
		builder.setShowProgress(showProgress)
	if performanceLogPath is not useJavaDefault:
		builder.setPerformanceLogFile(jvm.toFile(performanceLogPath))
	if sweepIncrement is not useJavaDefault:
		builder.setSweepIncrement(sweepIncrement)
	if maxNumMinimizations is not useJavaDefault:
		builder.setMaxNumMinimizations(maxNumMinimizations)
	if negligableFreeEnergy is not useJavaDefault:
		builder.setNegligableFreeEnergy(negligableFreeEnergy)

	return builder.build()


def SOFEA_MinLMFE(lmfe, numSequences, minFreeEnergyWidth):
	'''
	${class_javadoc(.sofea.MinLMFE)}

	# Arguments
	${args_java(.sofea.MinLMFE#<init>,
		[lmfe, objective],
		[numSequences],
		[minFreeEnergyWidth]
	)}

	# Returns
	${type_java(.sofea.MinLMFE)}
	'''
	return c.sofea.MinLMFE(lmfe, numSequences, minFreeEnergyWidth)


def SOFEA_SequenceLMFE(sequence, lmfe, minFreeEnergyWidth):
	'''
	${class_javadoc(.sofea.SequenceLMFE)}

	# Arguments
	${args_java(.sofea.SequenceLMFE#<init>,
		[sequence, seq],
		[lmfe],
		[minFreeEnergyWidth]
	)}

	# Returns
	${type_java(.sofea.SequenceLMFE)}
	'''
	return c.sofea.SequenceLMFE(sequence, lmfe, minFreeEnergyWidth)
