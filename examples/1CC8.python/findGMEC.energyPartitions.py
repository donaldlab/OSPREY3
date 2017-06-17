
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


def findGmec(epart):

	print("\n\nFinding GMEC with energy partition: %s\n" % epart)

	# compute the energy matrix with a specific energy partition
	fragEcalc = osprey.FragmentEnergyCalculator(confSpace, ffparams)
	emat = osprey.EnergyMatrix(confSpace, fragEcalc, energyPartition=epart)

	# find the best sequence and rotamers using the same energy partition
	confEcalc = osprey.ConfEnergyCalculator(fragEcalc, energyPartition=epart)
	astar = osprey.AStarTraditional(emat, confSpace)
	gmec = osprey.GMECFinder(confSpace, astar, confEcalc).find()


# find GMECs with different energy partitions
for epart in [osprey.EnergyPartition.Traditional, osprey.EnergyPartition.AllOnPairs]:
	findGmec(epart)

