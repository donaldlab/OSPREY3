
import osprey

osprey.start()

# run this example AFTER you've trained the three LUTE models using LUTE.train.py

# import our conf spaces from another file
execfile('LUTE.confSpaces.py')

# configure K*
kstar = osprey.KStar(
	confSpaces['protein'],
	confSpaces['ligand'],
	confSpaces['complex'],
	epsilon=0.01,
	writeSequencesToConsole=True,
	writeSequencesToFile='kstar.results.tsv'
)

# configure K* inputs for each conf space
for (id, confSpace) in confSpaces.items():

	# read the trained LUTE model
	pmat = osprey.DEE_read(confSpace, 'LUTE.pmat.%s.dat' % id)
	model = osprey.LUTE_read('LUTE.%s.dat' % id)

	# make a LUTE energy calculator from the model
	luteEcalc = osprey.LUTE_ConfEnergyCalculator(confSpace, model)

	# how should confs be ordered and searched? (don't forget to capture pmat,luteEcalc by using defaulted arguments)
	def makeAStar(rcs, pmat=pmat, luteEcalc=luteEcalc):
		return osprey.LUTE_AStar(rcs, pmat, luteEcalc, showProgress=False)

	# configure the K* settings for this conf space
	info = kstar.getConfSpaceInfo(confSpace)
	info.confEcalc = luteEcalc
	info.confSearchFactory = osprey.KStar.ConfSearchFactory(makeAStar)

# run K*
scoredSequences = kstar.run()
