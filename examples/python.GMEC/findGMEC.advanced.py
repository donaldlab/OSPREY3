
import osprey

osprey.start()

# what kind of hardware do we have?
parallelism = osprey.Parallelism(cpuCores=2)
# or use GPUs
gpuPrallelism = osprey.Parallelism(cpuCores=2, gpus=1, streamsPerGpu=1)

# choose a forcefield
ffparams = osprey.ForcefieldParams(osprey.Forcefield.AMBER)
ffparams.solvationForcefield = osprey.SolvationForcefield.EEF1 # this is the default
# or turn off solvation energy
#ffparams.solvationForcefield = None

# choose a template library
# (see templateLibrary.py for more detailed examples of custom template libraries)
templateLib = osprey.TemplateLibrary()

# load a molecule
mol = osprey.readPdb('1CC8.ss.pdb')

# define the protein strand
protein = osprey.Strand(mol, residues=['A2', 'A30'])

protein.flexibility['A2'].setLibraryRotamers('ALA', 'GLY')
protein.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE, 'VAL', 'ARG').setContinuous(10)
protein.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()

# make the conf space
confSpace = osprey.ConfSpace(protein)

# how should molecule energies be calculated?
ecalc = osprey.EnergyCalculator(
	confSpace,
	ffparams,
	parallelism=parallelism
)

# how could conformation energies be calculated?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# calculate the energy matrix
emat = osprey.EnergyMatrix(
	confEcalc,
	cacheFile='/tmp/emat.dat'
)

# define the conformation search
astar = osprey.AStarMPLP(
	emat,
	confSpace,
	numIterations=1
)
# or
traditionalAstar = osprey.AStarTraditional(emat, confSpace)

energyWindow = 0.1 # kcal/mol
confs = osprey.GMECFinder(
	astar,
	confEcalc,
	printIntermediateConfs=True,
	confLog='confs.txt'
).find(energyWindow)
