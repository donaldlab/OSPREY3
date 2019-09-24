
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

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

# define the design strand
design = osprey.Strand(mol, templateLib=templateLib, residues=['G648', 'G654'])
design.flexibility['G649'].setLibraryRotamers(osprey.WILD_TYPE, 'ALA', 'VAL').addWildTypeRotamers().setContinuous()
design.flexibility['G650'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
design.flexibility['G651'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['A155', 'A194'])
ligand.flexibility['A172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A192'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A193'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# make a single-state conf space:
confSpace = osprey.MultiStateConfSpace([
	osprey.StateMutable('complex', osprey.ConfSpace([design, ligand]))
])
complexState = confSpace.getState('complex')

# use the "wild-type" sequence
seq = confSpace.seqSpace.makeWildTypeSequence()

# (or any other sequence defined in your conf space)
#seq = confSpace.seqSpace.makeUnassignedSequence()\
#	.set('G649', 'ALA')

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

	# for this example, compute a single-state free energy
	# by using a trivial linear multi-state free energy (LMFE)
	lmfe = confSpace.lmfe()
		.addPositive('complex')
		.build(),

	# how precise do we want the free energy to be?
	minFreeEnergyWidth = 1.0
)

# BONUS:
# If you re-run this script, by default SOFEA will resume with
# the sweep line position from the previous computation.
# But if you've tried to re-run this script on a new sequence without clearing
# any of the intermediate files, SOFEA might get confused and quietly return a
# wrong answer for the new sequence! Or it might not. It depends.
# It would be nice to not delete the intermediate files though, and re-use as much
# of the previous computation(s) as possible to save time computing a new sequence.
# Luckily, to avoid SOFEA returning a wrong answer for a new sequence,
# we only need to reset the sweep line position between each computation.
# That way, you can re-run this script on new sequences and still get correct answers,
# but still re-use much of all the previous computations, and hopefully same some time.
resumeSweep = False

# finally, run the design!
sofea.refine(criterion, resumeSweep)

# when the design is finished, write out the results
with sofea.openSeqDB() as seqdb:

	print('\n\n')
	print('Sequence: %s' % seq)

	# show the free energy for our sequence
	g = criterion.getGBounds(seqdb, complexState)
	print('G: %s kcal/mol' % g)

	# show the pfunc value for our sequence too, if you like
	z = criterion.getZBounds(seqdb, complexState)
	print('Z: %s' % z)
