
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility[2].setLibraryRotamers('ALA', 'GLY')
strand.flexibility[3].setLibraryRotamers(osprey.WILD_TYPE, 'VAL')
strand.flexibility[4].setLibraryRotamers()

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how to compute the energy of a conformation?
fragEcalc = osprey.FragmentEnergyCalculator(confSpace, ffparams)
confEcalc = osprey.ConfEnergyCalculator(fragEcalc)

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confSpace, fragEcalc)
astar = osprey.AStarMPLP(emat, confSpace)

# find the best sequence and rotamers
gmec = osprey.GMECFinder(confSpace, astar, confEcalc).find()

