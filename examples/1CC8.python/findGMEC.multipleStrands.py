
import osprey

osprey.start()

mol = osprey.readPdb('1CC8.ss.pdb')

# define a strand
strand1 = osprey.Strand(mol, residues=[2,20])
strand1.flexibility[2].setLibraryRotamers().addWildTypeRotamers()
strand1.flexibility[3].setLibraryRotamers().addWildTypeRotamers()
strand1.flexibility[4].setLibraryRotamers().addWildTypeRotamers()

# define another strand
strand2 = osprey.Strand(mol, residues=[21,40])
strand2.flexibility[21].setLibraryRotamers().addWildTypeRotamers()
strand2.flexibility[22].setLibraryRotamers().addWildTypeRotamers()
strand2.flexibility[23].setLibraryRotamers().addWildTypeRotamers()

# make the conf space
confSpace = osprey.ConfSpace([strand1, strand2])

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

