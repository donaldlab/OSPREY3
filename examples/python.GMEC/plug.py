
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility['A23'].setLibraryRotamers('ASN').setContinuous()
strand.flexibility['A27'].setLibraryRotamers('THR').setContinuous()
strand.flexibility['A36'].setLibraryRotamers('ILE').setContinuous()

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# what's the biggest atom overlap we should allow?
plugThreshold = 0.4 # in angstroms

# run PLUG, using the Dead-End Elimination (DEE) tool
pmat = osprey.DEE(
	confSpace,

	# don't need an energy matrix if we're just using PLUG, so leave it None
	emat=None,
	singlesThreshold=None, # no energy matrix means we should disable threshold pruning too
	pairsThreshold=None, # (which is on by default)

	# or add your energy matrix and turn on the DEE options, it's up to you
	# the various pruning methods are all inter-compatible

	# turn on PLUG to prune singles and pairs
	singlesPlugThreshold=plugThreshold,
	pairsPlugThreshold=plugThreshold,

	# add parallelism for more speed
	parallelism=osprey.Parallelism(cpuCores=2),

	showProgress=True
)

# pruning matrix has been computed now!
# use it in your designs like you would any pruning matrix
