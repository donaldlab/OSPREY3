
import osprey
osprey.start()

# let's go more faster
parallelism = osprey.Parallelism(cpuCores=2)

# if you have GPUs, they can make LUTE training go MUCH faster
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

# make a multi-state conf space
confSpace = osprey.MultiStateConfSpace([
	osprey.StateMutable('design', osprey.ConfSpace(design)),
	osprey.StateUnmutable('ligand', osprey.ConfSpace(ligand)),
	osprey.StateMutable('complex', osprey.ConfSpace([design, ligand]))
])


# train a LUTE model for each state
for state in confSpace.states:

	# how should we compute energies of molecules?
	ecalc = osprey.EnergyCalculator(state.confSpace, ffparams, parallelism=parallelism)

	# how should we define energies of conformations?
	eref = osprey.ReferenceEnergies(state.confSpace, ecalc)
	confEcalc = osprey.ConfEnergyCalculator(state.confSpace, ecalc, referenceEnergies=eref)

	# calculate the energy matrix so we can do DEE
	emat = osprey.EnergyMatrix(
		confEcalc,
		# save it to disk so we don't have to calculate it again later
		cacheFile = 'sofea.%s.emat' % state.name
	)

	# Good LUTE fits rely heavily on DEE pruning, so use the best available DEE options here.
	# With interval pruning, we prune only conformations whose energies are more than X higher
	# than the lower bound of the conformation with the minimum lower bound.
	# If we sort the whole conformation space by lower bounds, interval pruning is like keeping
	# only conformations in the bottom X kcal/mol slice of the sorted list.
	pruningInterval = 20.0 # kcal/mol
	pmat = osprey.DEE(
		state.confSpace,
		emat,
		# pruning triples is much slower than just singles and pairs, but important for good fits
		singlesGoldsteinDiffThreshold = pruningInterval,
		pairsGoldsteinDiffThreshold = pruningInterval,
		triplesGoldsteinDiffThreshold = pruningInterval,
		# transitive pruning (ie, pruning tuples implied by other tuples) is also important for LUTE
		singlesTransitivePruning = True,
		pairsTransitivePruning = True,
		triplesTransitivePruning = True,
		# can use PLUG too if you want, but it's usually pretty slow, especially the triples
		#singlesPlugThreshold = 0.6,
		#pairsPlugThreshold = 0.6,
		#triplesPlugThreshold = 0.6,
		parallelism = parallelism,
		showProgress = True,
		# save it to disk so we don't have to calculate it again later
		cacheFile = 'sofea.%s.pmat' % state.name
	)

	# train the LUTE model (if needed)
	if not osprey.jvm.toFile('sofea.%s.lute' % state.name).exists():

		print('\nTraining %s LUTE model...\n' % state.name)
		model = osprey.LUTE_train(
			confEcalc,
			emat,
			pmat,

			# what's the highest RMS error you're willing tolerate for the model fitting?
			maxRMSE = 0.1, # kcal/mol

			# how much overfitting are you willing to tolerate?
			# (overfitting score = training set RMSE / test set RMSE)
			maxOverfittingScore = 1.5, # unitless ratio

			# LUTE uses a conformation database to speed up training.
			# You can save this database to disk by choosing a filename, or don't set this argument to not save the db.
			# If you save it, the datebase can be used to speed up the next training job too, as well as this one.
			# But if you change the conformation space at all, you'll need to delete the db file before training again.
			confDBPath = 'sofea.%s.confdb' % state.name
		)

		# if the LUTE fit isn't good enough, don't bother continuing with the rest of the design
		if model is None:
			raise Exception("LUTE model does not meet accuracy goals. Try changing parameters and train again?")

		# save the model to a file
		osprey.LUTE_write(model, 'sofea.%s.lute' % state.name)


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

	# get the 2 sequences with the lowest LMFE
	numSequences = 2,

	# free energy calculations require high-precision arithmetic
	# copy the same math settings used by SOFEA itself
	mathContext = sofea.mathContext
)

# finally, run the design!
sofea.refine(criterion)

# when the design is finished, write out the results
with sofea.openSeqDB() as seqdb:
	criterion.makeResultDoc(seqdb, osprey.jvm.toFile('sofea.md'))
