
import osprey
osprey.start()

# let's go more faster
parallelism = osprey.Parallelism(cpuCores=2)

# if you have GPUs, designs can go much faster
#parallelism = osprey.Parallelism(cpuCores=2, gpus=1, streamsPerGpu=16)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# define the conformation space
strand = osprey.Strand(mol, residues=['G648', 'G654'])
for resNum in ['G649', 'G650', 'G651']:
	strand.flexibility[resNum]\
		.setLibraryRotamers(osprey.WILD_TYPE)\
		.addWildTypeRotamers()\
		.setContinuous()
confSpace = osprey.MultiStateConfSpace([
	# yeah, we didn't define any mutations,
	# but SOFEA still needs one "mutable" state to define the sequence space
	osprey.StateMutable('protein', osprey.ConfSpace(strand))
])

# we only care about the "wild-type" sequence in this case
seq = confSpace.seqSpace.makeWildTypeSequence()

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams, parallelism=parallelism)

# make a function that makes a SOFEA config object for a state
def config(state):

	# how should we define energies of conformations?
	confEcalc = osprey.ConfEnergyCalculator(
		state.confSpace,
		ecalc,
		referenceEnergies = osprey.ReferenceEnergies(state.confSpace, ecalc),

		# use the "AllOnPairs" setting to get tighter energy bounds
		# tighter energy bounds make SOFEA run much faster!
		energyPartition = osprey.EnergyPartition.AllOnPairs
	)

	# train forcefield approximations to make energy calculations go much faster
	amat = osprey.ApproximatorMatrix(
		confEcalc,

		# save the approximator matrix to disk so we don't have to calculate it again later
		cacheFile = 'sofea.%s.amat' % state.name,
	)

	# update the conformation energy calculator with the forcefield approximations
	confEcalc = osprey.ConfEnergyCalculator(
		state.confSpace,
		ecalc,

		# copy settings from the previous confEcalc
		referenceEnergies = confEcalc.eref,
		energyPartition = confEcalc.epart,

		# give the new confEcalc the approximator matrix
		amat = amat,

		# how much error can we tolerate per conformation?
		approximationErrorBudget = 0.1 # kcal/mol
	)

	# make the SOFEA config
	return osprey.SOFEA_StateConfig(

		# calculate the energy matrix
		emat = osprey.EnergyMatrix(
			confEcalc,

			# save the emat to disk so we don't have to calculate it again later
			cacheFile = 'sofea.%s.emat' % state.name,

			# correct overly optimistic energy bounds with triples to get even tighter energy bounds!
			tripleCorrectionThreshold = 0.0
		),

		confEcalc = confEcalc
	)

# configure the SOFEA algorithm
sofea = osprey.SOFEA(
	confSpace = confSpace,
	configFunc = config,
	fringedbLowerPath = 'sofea.lower.fringedb',
	fringedbLowerMiB = 100,
	fringedbUpperPath = 'sofea.upper.fringedb',
	fringedbUpperMiB = 10,
	rcdbPath = 'rc.db', # set a path to collect the residue conformation (RC) info
	seqdbPath = 'sofea.seqdb'
)

# start new databases if this is the first time running this script
if not osprey.jvm.toFile('sofea.seqdb').exists():
	sofea.init()

# otherwise, keep the existing databases
# if you're resuming a previous design, don't call init().
# SOFEA will try to delete your databases if you do!

# choose a termination criterion for SOFEA
# ie, what are we trying to compute?
criterion = osprey.SOFEA_SequenceLMFE(

	sequence = seq,

	lmfe = confSpace.lmfe()
		.addPositive('protein')
		.build(),

	# how precise do we want each free energy to be?
	minFreeEnergyWidth = 0.1
)

# finally, run the design!
sofea.refine(criterion)

# when the design is finished, examine the results
print('\n')
state = confSpace.getState('protein')

with sofea.openSeqDB() as seqdb:

	# get the free energy for the whole protein
	totalg = criterion.getGBounds(seqdb, state, sofea.mathContext)
	totalz = criterion.getZBounds(seqdb, state)
	print('protein   G = [%9.4f,%9.4f]   Z = [%.4e,%.4e]' % (
		totalg.lower, totalg.upper,
		totalz.lower.doubleValue(), totalz.upper.doubleValue()
	))


with sofea.openRCDB() as rcdb:

	# loop over the residue conformations (RCs) info to get residue entropies
	for pos in state.confSpace.positions:
		print('%s' % pos.resNum)
		for rc in pos.resConfs:
			rcg = rcdb.getGBounds(state, seq, pos, rc, sofea.mathContext)
			rcz = rcdb.getZSumBounds(state, seq, pos, rc)
			print("\t%s:%s   G = [%9.4f,%9.4f]   z = [%.4e,%.4e]" % (
				rc.template.name, rc.getRotamerCode(),
				rcg.lower, rcg.upper,
				rcz.lower.doubleValue(), rcz.upper.doubleValue()
			))
