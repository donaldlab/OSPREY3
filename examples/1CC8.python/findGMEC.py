
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility[2].setLibraryRotamers('ALA', 'GLY');
strand.flexibility[3].setLibraryRotamers(osprey.WILD_TYPE, 'VAL');
strand.flexibility[4].setLibraryRotamers();

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confSpace, ffparams, cacheFile='/tmp/emat.dat')
astar = osprey.AStarMPLP(emat, confSpace)

# how to compute the energy of a conformation?
ecalc = osprey.ConfEnergyCalculator(confSpace, ffparams)

# find the best sequence and rotamers
gmec = osprey.GMECFinder(confSpace, astar, ecalc).find()

