
import osprey

osprey.start()

# TODO: make all of this actually work

# what kind of hardware do we have?
parallelism = osprey.Parallelism(cpuCores=2, gpus=1, streamsPerGpu=8)

# load a molecule
mol = osprey.loadPdb('1CC8.ss.pdb')

# define the strands
strand = osprey.FlexibleStrand(mol, 2, 73)
strand.flexibility={
	2: ['ALA', 'GLY'],
	3: [],
	4: [],
	5: []
}
strand.addWTAminoAcids()

# load a rotamer library
rotlib = osprey.loadLovellRotamerLibrary()
# or maybe:
#rotlib = osprey.loadRotamerLibrary('/path/to/rots')

# make the conformation space
confSpace = osprey.ConfSpace(
	strands=[strand],
	rotlib=rotlib
)
confSpace.addWTRotamers()

# pick a forcefield
forcefield = osprey.makeForcefield()

# define "energy" for a conformation
ecalc = osprey.RigidConfEnergyCalc(
	forcefield=forcefield,
	confSpace=confSpace
)
# or
ecalc = osprey.MinimizingConfEnergyCalc(
	forcefield=forcefield,
	confSpace=confSpace
)

# load/compute the energy matrix
emat = osprey.calcEmat(
	confSpace=confSpace,
	ecalc=ecalc,
	parallelism=paralleism
	cachePath='emat.dat'
)

# load/compute the pruning matrix
pmat = osprey.calcPmat(
	confSpace=confSpace,
	emat=emat,
	parallelism=parallelism,
	cachePath='pmat.dat'
)

# how to rank the conformations?
confSearcher = osprey.AstarConfSearcher(
	emat=emat,
	pmat=pmat,
	scorer=osprey.TraditionalAstarScorer(emat)
)
# or
confSearcher = osprey.AstarConfSeacher(
	emat=emat,
	pmat=pmat,
	scorer=osprey.MPLPAstarScorer(emat, numIters=1)
)

# find the gmec
conf = osprey.GMECFinder.findGMEC(
	confSearcher=confSearcher,
	ecalc=ecalc,
	parallelism=parallelism,
	printIntermediateConfs=True
)
print(conf)

# or find the GMEC and an energy window
confs = osprey.GMECFinder.findGMECAndWindow(
	confSearcher=confSearcher,
	ecalc=ecalc,
	energyWindow=3, # kcal/mol
	parallelism=parallelism,
	printIntermediateConfs=True
)
for conf in confs:
	print(conf)

