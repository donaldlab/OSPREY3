
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility['A2'].setLibraryRotamers('ALA', 'GLY')
strand.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE, 'VAL')
strand.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE)

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# compute the reference energies
ecalc = osprey.EnergyCalculator(confSpace, ffparams)
eref = osprey.ReferenceEnergies(confSpace, ecalc)

# how to compute energy for a conformation?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# use the reference energies for energy matrix computation
emat = osprey.EnergyMatrix(confEcalc)

# Optional: using edge MPLP when using reference energies can make A* faster
astar = osprey.AStarMPLP(emat, confSpace, updater=osprey.EdgeUpdater())

# find the best sequence and rotamers
gmec = osprey.GMECFinder(astar, confEcalc).find()


