
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

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

# how should we define energies of conformations?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# compute the energy matrix
emat = osprey.EnergyMatrix(confEcalc)

# run DEE with just steric pruning
pmat = osprey.DEE(confSpace, emat, showProgress=True)

# or run DEE with Goldstein pruning
#i0 = 10.0 # kcal/mol
#pmat = osprey.DEE(confSpace, emat, showProgress=True, singlesGoldsteinDiffThreshold=i0, pairsGoldsteinDiffThreshold=i0)

# find the best sequence and rotamers
astar = osprey.AStarMPLP(emat, pmat)
gmec = osprey.GMECFinder(astar, confEcalc).find()
