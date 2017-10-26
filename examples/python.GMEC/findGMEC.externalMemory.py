
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility['A2'].setLibraryRotamers(osprey.WILD_TYPE).setContinuous()
strand.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE).setContinuous()
strand.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE).setContinuous()

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how to compute the energy of a conformation?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# configure external memory settings
osprey.initExternalMemory(64)

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confEcalc)
astar = osprey.AStarMPLP(emat, confSpace, useExternalMemory=True)

# find the best sequence and rotamers
gmec = osprey.GMECFinder(astar, confEcalc, useExternalMemory=True).find()
