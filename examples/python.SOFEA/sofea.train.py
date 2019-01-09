
# this is the script that trains the LUTE model for your design
# run this script first!


import osprey

osprey.start()

# let's go more faster
parallelism = osprey.Parallelism(cpuCores=2)

# if you have GPUs, they can make LUTE training go MUCH faster
#parallelism = osprey.Parallelism(cpuCores=2, gpus=1, streamsPerGpu=16)


# import our conf space from another file
execfile('sofea.confspace.py')


# train each state separately
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
		triplesGoldsteinDiffThreshold=pruningInterval,
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

	# train a LUTE model
	print('\nTraining %s model...\n' % state.name)
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

	# save the model to a file
	osprey.LUTE_write(model, 'sofea.%s.lute' % state.name)


# While LUTE training runs, a bunch of info will be printed to the console
# including progress info, and fit quality scores.
# When the training is finished, look at the quality of the fit
# and decide if you want to run the design using that LUTE model.

# if no, tweak the training parameters and run the training script again

# if yes, continue to the sofea.design.py script
