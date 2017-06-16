
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
fragEcalc = osprey.FragmentEnergyCalculator(confSpace, ffparams)
eref = osprey.ReferenceEnergies(confSpace, fragEcalc)

# use the reference energies for energy matrix computation
emat = osprey.EnergyMatrix(confSpace, fragEcalc, referenceEnergies=eref)

# Optional: using edge MPLP when using reference energies can make A* faster
astar = osprey.AStarMPLP(emat, confSpace, updater=osprey.EdgeUpdater())

# find the best sequence and rotamers
confEcalc = osprey.ConfEnergyCalculator(fragEcalc, referenceEnergies=eref)
gmec = osprey.GMECFinder(confSpace, astar, confEcalc).find()


