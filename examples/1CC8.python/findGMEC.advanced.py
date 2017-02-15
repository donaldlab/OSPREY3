
import osprey

osprey.start()

# what kind of hardware do we have?
#parallelism = osprey.Parallelism(cpuCores=2, gpus=1, streamsPerGpu=1)
parallelism = osprey.Parallelism(cpuCores=2)

# choose a forcefield
ff = osprey.Forcefield.AMBER

# choose a template library
templateLib = osprey.TemplateLibrary(
	forcefield=ff,
	rotamers=osprey.LovellRotamers
	# or
	#templateCoords='path/to/coords/file'
)

# load a molecule
mol = osprey.readPdb('1CC8.ss.pdb')

# define the protein strand
protein = osprey.Strand(mol, residues=[2, 30])

protein.flexibility[2].setLibraryRotamers('ALA', 'GLY');
protein.flexibility[3].setLibraryRotamers(osprey.WILD_TYPE, 'VAL', 'ARG').setContinuous(10);
protein.flexibility[4].setLibraryRotamers().addWildTypeRotamers();

# make the conf space
confSpace = osprey.ConfSpace(protein)

# pick forcefield params
ffparams = osprey.ForcefieldParams()
ffparams.solvationForcefield = None # turn off solvation enery
# or
#ffparams.solvationForcefield = osprey.SolvationForcefield.EEF1

# calculate the energy matrix
emat = osprey.EnergyMatrix(
	confSpace,
	ffparams=ffparams,
	parallelism=parallelism,
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
	confSpace,
	astar=astar,
	energyCalculator=osprey.MinimizingEnergyCalculator(
		confSpace,
		ffparams=ffparams,
		parallelism=parallelism
	),
	printIntermediateConfs=True,
	confLog='confs.txt'
).find(energyWindow)

# get some info about the confs
print('\nyup, we found %d confs in the energy window' % len(confs))
gmec = confs.get(0)
gmecMol = confSpace.makeMolecule(gmec)
osprey.writePdb(gmecMol, 'gmec.pdb')

