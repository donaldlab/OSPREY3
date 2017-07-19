
import osprey

osprey.start()

# what kind of hardware do we have?
parallelism = osprey.Parallelism(cpuCores=2)
# or use GPUs
#parallelism = osprey.Parallelism(cpuCores=2, gpus=1, streamsPerGpu=1)

# choose a forcefield
ffparams = osprey.ForcefieldParams(osprey.Forcefield.AMBER)
ffparams.solvationForcefield = osprey.SolvationForcefield.EEF1 # this is the default
# or turn off solvation energy
#ffparams.solvationForcefield = None

# choose a template library
templateLib = osprey.TemplateLibrary(
	forcefield=ffparams.forcefld,
	rotamers=osprey.LovellRotamers
	# or
	#templateCoords='path/to/coords/file'
)

# load a molecule
mol = osprey.readPdb('1CC8.ss.pdb')

# define the protein strand
protein = osprey.Strand(mol, residues=[2, 30])

protein.flexibility[2].setLibraryRotamers('ALA', 'GLY')
protein.flexibility[3].setLibraryRotamers(osprey.WILD_TYPE, 'VAL', 'ARG').setContinuous(10)
protein.flexibility[4].setLibraryRotamers().addWildTypeRotamers()

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
#astar = osprey.AStarTraditional(emat, confSpace)

energyWindow = 0.1 # kcal/mol
confs = osprey.GMECFinder(
	astar,
	confEcalc,
	printIntermediateConfs=True,
	confLog='confs.txt'
).find(energyWindow)

# get some info about the confs
print('\nyup, we found %d confs in the energy window' % confs.size())
gmec = confs.peek()
gmecMol = confSpace.makeMolecule(gmec)
osprey.writePdb(gmecMol, 'gmec.pdb')

