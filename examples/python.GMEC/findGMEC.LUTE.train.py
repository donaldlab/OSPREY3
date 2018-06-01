
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
for resNum in ['A2', 'A3', 'A4', 'A5']:
	strand.flexibility[resNum].setLibraryRotamers('ALA', 'VAL', 'ILE', 'LEU').setContinuous()

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# let's go fast
parallelism = osprey.Parallelism(cpuCores=4)

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams, parallelism=parallelism)

# how should we define energies of conformations?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# calculate the energy matrix so we can do DEE
emat = osprey.EnergyMatrix(
	confEcalc,
	cacheFile='LUTE.emat.dat' # save it to disk so we don't have to calculate it again later
)

# Good LUTE fits rely heavily on DEE pruning, so use the best available DEE options here.
# With interval pruning, we prune only conformations whose energies are more than X higher
# than the lower bound of the conformation with the minimum lower bound.
# If we sort the whole conformation space by lower bounds, interval pruning is like keeping
# only conformations in the bottom X kcal/mol slice of the sorted list.
pruningInterval = 10.0 # kcal/mol
pmat = osprey.DEE(
	confSpace,
	emat,
	singlesGoldsteinDiffThreshold=pruningInterval,
	pairsGoldsteinDiffThreshold=pruningInterval,
	triplesGoldsteinDiffThreshold=pruningInterval, # pruning triples is much slower, but important for good fits
	parallelism=parallelism,
	showProgress=True,
	cacheFile='LUTE.pmat.dat' # save it to disk so we don't have to calculate it again later
)

# train a LUTE model
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
	confDBPath = 'LUTE.conf.db'
)

# save the model to a file
osprey.LUTE_write(model, 'LUTE.dat')

# While LUTE training runs, a bunch of info will be printed to the console
# including progress info, and fit quality scores.
# When the training is finished, look at the quality of the fit
# and decide if you want to run the design using that LUTE model.

# if no, tweak the training parameters and run the training script again

# if yes, continue to the findGMEC.LUTE.design.py script
