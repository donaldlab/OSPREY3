
import osprey

osprey.start()

# define a strand
# let's add some big mutations so the search takes a noticeable amount of time
strand = osprey.Strand('1CC8.ss.pdb')
for resNum in ['A2', 'A3', 'A4']:
	strand.flexibility[resNum].setLibraryRotamers(osprey.WILD_TYPE, 'ARG', 'LYS')

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

# how should we define energies of conformations?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)
emat = osprey.EnergyMatrix(confEcalc)

# let's use an energy window so the search takes a noticeable amount of time
energyWindow = 5.0 # kcal/mol

# find the best sequence and rotamers
# and store conformations in a database
confDBFile = 'gmec.conf.db'
astar = osprey.AStarMPLP(emat, confSpace)
gmec = osprey.GMECFinder(astar, confEcalc, confDBFile=confDBFile).find(energyWindow)


# If OSPREY exits unexpectedly, the 'gmec.conf.db' file will still remain.
# running the design again using an existing confDB should be MUCH faster!
# meaning, your design should catch up to where it left off very quickly.
astar = osprey.AStarMPLP(emat, confSpace)
gmec2 = osprey.GMECFinder(astar, confEcalc, confDBFile=confDBFile).find(energyWindow)
