
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



#------------------------------------
# pythonic wrappers for Java builders
#------------------------------------

def Parallelism(cpuCores=None, gpus=None, streamsPerGpu=None):
	'''
	Specifies how Osprey should use available hardware

	:param int cpuCores: :java:fielddoc:`.parallelism.Parallelism#numThreads`
	:param int gpus: :java:fielddoc:`.parallelism.Parallelism#numGpus`
	:param int streamsPerGpu: :java:fielddoc:`.parallelism.Parallelism#numStreamsPerGpu`
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
	Creates a template library of rotamers for molecule residues.

	:param forcefield: :java:methoddoc:`.restypes.GenericResidueTemplateLibrary$Builder#setForcefield`
	:type forcefield: :java:ref:`.energy.forcefield.ForcefieldParams$Forcefield`
	:default forcefield: osprey.Forcefield.AMBER

	:param templateCoords: :java:methoddoc:`.restypes.GenericResidueTemplateLibrary$Builder#setTemplateCoords`

		.. todo:: explain template file format?

	:type templateCoords: coords str or file path
	:default templateCoords: <natural amino acids>

	:param rotamers: :java:methoddoc:`.restypes.GenericResidueTemplateLibrary$Builder#setRotamers`

		.. todo:: explain rotamer file format?

	:type rotamers: coords str or file path
	:default rotamers: <Lovell rotamer library>

	:param backboneDependentRotamers: :java:methoddoc:`.restypes.GenericResidueTemplateLibrary$Builder#setBackboneDependentRotamers`

		.. todo:: explain backbone-dependent rotamer file format?

	:type backboneDependentRotamers: coords str or file path

	:rtype: :java:ref:`.restypes.GenericResidueTemplateLibrary`
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
	A molecule with associated flexibilty information

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
	Make a configuration space of strands

	:param strands: the strands to use
	:type strands: :java:ref:`.confspace.Strand` or list of Strands

	:param float shellDist: A residue is included in the steric shell if any
		of its atoms lies within ``shellDist`` angstroms of any atom in any flexible residue.
	:default shellDist: float('inf')

	:rtype: :java:ref:`.confspace.SimpleConfSpace`
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
	Options for configuring forcefields for energy calculation.

	:rtype: :java:ref:`.energy.forcefield.ForcefieldParams`
	'''
	return c.energy.forcefield.ForcefieldParams()


def EnergyMatrix(confSpace, ffparams=None, parallelism=None, cacheFile=None):
	'''
	Compute a matrix of energies between pairs of residue conformations to be used by A* search.

	:param confSpace: The conformation space containing the strands to be designed.

		If the strands are configured with continuous flexibility, the energy matrix will
		minimize residue conformation pairs before computing energies.

	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

	:param ffparams: The forcefield parameters for energy calculation.
	:type ffparams: :java:ref:`.energy.forcefield.ForefieldParams`
	:default ffparams: osprey.ForcefieldParams()

	:param parallelism: Available hardware for high-performance computation.
	:type parallelism: :java:ref:`.parallelism.Parallelism`
	:default parallelism: ospesy.Parallelism()

	:param str cacheFile: Path to file where energy matrix should be saved between computations.

		.. note:: Energy matrix computation can take a long time, but often the results
			can be reused between computations. Use a cache file to skip energy matrix
			computation on the next Osprey run if the energy matrix has already been
			computed once before.
			
		.. warning:: If design settings are changed between runs, Osprey will make
			some effort to detect that the energy matrix cache is out-of-date and compute a
			new energy matrix instead of usng the cached, incorrect one. Osprey might not detect
			all design changes though, and incorrectly reuse a cached energy matrix, so it
			is best to manually delete the entry matrix cache file after changing design settings.
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

	:param emat: The energy matrix to use for pairwise residue conformation energies.
	:type emat: :java:ref:`.ematrix.EnergyMatrix`

	:param confSpace: The conformation space containing the residue conformations to search.
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

	:rtype: :java:ref:`.astar.conf.ConfAStarTree`
	'''
	builder = c.astar.conf.ConfAStarTree.builder(emat, confSpace)
	builder.setTraditional()
	return builder.build()


def AStarMPLP(emat, confSpace, numIterations=None, convergenceThreshold=None):
	'''
	Creates an A* search using the newer estimation function based on Max Product Linear Programming (MPLP)

	:param emat: The energy matrix to use for pairwise residue conformation energies.
	:type emat: :java:ref:`.ematrix.EnergyMatrix`

	:param confSpace: The conformation space containing the residue conformations to search.
	:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

	:builder_option numIterations .astar.conf.ConfAStarTree$MPLPBuilder setNumIterations numIterations:
	:builder_option convergenceThreshold .astar.conf.ConfAStarTree$MPLPBuilder setConvergenceThreshold convergenceThreshold:
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

