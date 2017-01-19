
import osprey

osprey.start()

# load a molecule
mol = osprey.loadPdb('1CC8.ss.pdb')

# TODO: make the rest of this actually work

# define the strands
strand = osprey.FlexibleStrand(mol, 2, 73)
strand.flexibility={
	2: ['ALA', 'GLY'],
	3: [],
	4: [],
	5: []
}
strand.addWTAminoAcids()

# make the conformation space
confSpace = osprey.ConfSpace([strand])

# define "energy" for a conformation
ecalc = osprey.RigidConfEnergyCalc(confSpace)

# compute the energy matrix
emat = osprey.calcEmat(confSpace, ecalc)

# find the gmec
conf = osprey.GMECFinder.findGMEC(
	confSearcher=osprey.AstarConfSearcher(emat),
	ecalc=ecalc
)

# print the result
print(conf)

