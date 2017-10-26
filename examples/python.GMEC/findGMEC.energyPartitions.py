
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
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

def findGmec(epart):

	print("\n\nFinding GMEC with energy partition: %s\n" % epart)

	# configure conformation energies with the desired energy partition
	confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc, energyPartition=epart)

	# find the GMEC
	emat = osprey.EnergyMatrix(confEcalc)
	astar = osprey.AStarTraditional(emat, confSpace)
	gmec = osprey.GMECFinder(astar, confEcalc).find()


# find GMECs with different energy partitions
for epart in [osprey.EnergyPartition.Traditional, osprey.EnergyPartition.AllOnPairs]:
	findGmec(epart)

