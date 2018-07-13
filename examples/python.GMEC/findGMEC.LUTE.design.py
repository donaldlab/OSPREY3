
import osprey

osprey.start()

# run this example AFTER you've trained a LUTE model using findGMEC.LUTE.design.py

# define the same conformation space as the training script
strand = osprey.Strand('1CC8.ss.pdb')
for resNum in ['A2', 'A3', 'A4', 'A5']:
	strand.flexibility[resNum].setLibraryRotamers('ALA', 'VAL', 'ILE', 'LEU').setContinuous()
confSpace = osprey.ConfSpace(strand)

# read the trained LUTE model
pmat = osprey.DEE_read(confSpace, 'LUTE.pmat.dat')
model = osprey.LUTE_read('LUTE.dat')

# find the GMEC using the LUTE model
gmecFinder = osprey.LUTE_GMECFinder(confSpace, model, pmat, printIntermediateConfs=True)
gmec = gmecFinder.find(1.0)
