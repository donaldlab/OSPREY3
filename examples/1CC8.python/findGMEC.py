
import osprey

osprey.start()

# define a strand
strand = osprey.makeStrand('1CC8.ss.pdb')

strand.flexibility[2].setLibraryRotamers('ALA', 'GLY');
strand.flexibility[3].setLibraryRotamers(osprey.WildType, 'VAL');
strand.flexibility[4].setLibraryRotamers();

# TODO: make the rest of this actually work
exit()

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

