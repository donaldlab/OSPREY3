
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility[2].setLibraryRotamers().setContinuous()
strand.flexibility[3].setLibraryRotamers().setContinuous()
strand.flexibility[4].setLibraryRotamers().setContinuous()

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how to compute the energy of a conformation?
fragEcalc = osprey.FragmentEnergyCalculator(confSpace, ffparams)
confEcalc = osprey.ConfEnergyCalculator(fragEcalc)

# configure external memory settings
osprey.ExternalMemory.setInternalLimit(64)

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confSpace, fragEcalc)
astar = osprey.AStarMPLP(emat, confSpace, useExternalMemory=True)

# find the best sequence and rotamers
gmec = osprey.GMECFinder(confSpace, astar, confEcalc, useExternalMemory=True).find()

