
# this is the script that runs the SOFEA algorithm using a trained LUTE model
# run this script second!


import osprey

osprey.start()

# let's go more faster
parallelism = osprey.Parallelism(cpuCores=2)
# sadly, OSPREY can't use GPUs to speed up this part

# import our conf space from another file
execfile('sofea.confspace.py')

# make a function that makes a SOFEA config object for a state
def config(state):

	# make a LUTE energy calculator from the trained LUTE model
	luteEcalc = osprey.LUTE_ConfEnergyCalculator(
		state.confSpace,
		osprey.LUTE_read('sofea.%s.lute' % state.name)
	)

	# read the pruning matrix
	pmat = osprey.DEE_read(state.confSpace, 'sofea.%s.pmat' % state.name)

	# make the SOFEA config
	return osprey.SOFEA_StateConfig(luteEcalc, pmat)

# configure the SOFEA algorithm
sofea = osprey.SOFEA(
	confSpace = confSpace,
	configFunc = config,
	parallelism = parallelism,
	fringeDBPath = 'sofea.fringedb',
	fringeDBSizeMiB = 100,
	seqDBPath = 'sofea.seqdb'
)

# start new databases if this is the first time running this script
if not osprey.jvm.toFile('sofea.fringedb').exists():
	sofea.init()

# otherwise, keep the existing databases
# if you're resuming a previous design, don't call init().
# SOFEA will try to delete your databases if you do!

# choose a termination criterion for SOFEA
# ie, what are we trying to compute?
criterion = osprey.SOFEA_MinLMFE(

	# for this example, minimize a linear multi state free energy (LMFE) that encodes binding affinity
	lmfe = confSpace.lmfe()
		.addPositive('complex')
		.addNegative('design')
		.addNegative('ligand')
		.build(),

	# get the 4 sequences with the lowest LMFE
	numSequences = 4,

	# free energy calculations require high-precision arithmetic
	# copy the same math settings used by SOFEA itself
	mathContext = sofea.mathContext
)

# finally, run the design!
sofea.refine(criterion)

# when the design is finished, write out the results
with sofea.openSeqDB() as seqdb:
	criterion.makeResultDoc(seqdb, osprey.jvm.toFile('sofea.md'))
