
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

# compute the reference energies
ecalc = osprey.EnergyCalculator(confSpace, ffparams)
eref = osprey.ReferenceEnergies(confSpace, ecalc)

# use the reference energies for energy matrix computation
ecalc = osprey.EnergyCalculator(confSpace, ffparams, referenceEnergies=eref)
emat = osprey.EnergyMatrix(confSpace, ecalc)

# Optional: using edge MPLP when using reference energies can make A* faster
astar = osprey.AStarMPLP(emat, confSpace, updater=osprey.EdgeUpdater())

# find the best sequence and rotamers
gmec = osprey.GMECFinder(confSpace, astar, ecalc).find()


