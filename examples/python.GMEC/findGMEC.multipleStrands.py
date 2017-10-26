
import osprey

osprey.start()

mol = osprey.readPdb('1CC8.ss.pdb')

# define a strand
strand1 = osprey.Strand(mol, residues=['A2','A20'])
strand1.flexibility['A2'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()
strand1.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()
strand1.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()

# define another strand
strand2 = osprey.Strand(mol, residues=['A21','A40'])
strand2.flexibility['A21'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()
strand2.flexibility['A22'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()
strand2.flexibility['A23'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()

# make the conf space
confSpace = osprey.ConfSpace([strand1, strand2])

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how to compute the energy of a conformation?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confEcalc)
astar = osprey.AStarMPLP(emat, confSpace)

# find the best sequence and rotamers
gmec = osprey.GMECFinder(astar, confEcalc).find()

