
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility[2].setLibraryRotamers('ALA', 'GLY');
strand.flexibility[3].setLibraryRotamers(osprey.WILD_TYPE, 'VAL');
strand.flexibility[4].setLibraryRotamers();

# make the conf space
confSpace = osprey.ConfSpace(strand)

# get an energy matrix
emat = osprey.EnergyMatrix(confSpace, cacheFile='/tmp/emat.dat')

# find the best sequence and rotamers
gmec = osprey.GMECFinder(confSpace, emat, printIntermediateConfs=True).find(1);

