
import osprey

osprey.start()

# run this example AFTER you've trained the three LUTE models using LUTE.train.py

# import our conf spaces from another file
exec(open("LUTE.confSpaces.py").read())

# configure BBK*
bbkstar = osprey.BBKStar(
	confSpaces['protein'],
	confSpaces['ligand'],
	confSpaces['complex'],
	epsilon=0.01,
	numBestSequences=2,
	writeSequencesToConsole=True
)

# configure BBK* inputs for each conf space
for info in bbkstar.confSpaceInfos():

	# read the trained LUTE model
	pmat = osprey.DEE_read(info.confSpace, 'LUTE.pmat.%s.dat' % info.id)
	model = osprey.LUTE_read('LUTE.%s.dat' % info.id)

	# make a LUTE energy calculator from the model
	luteEcalc = osprey.LUTE_ConfEnergyCalculator(info.confSpace, model)
	info.confEcalcMinimized = luteEcalc

	# how should confs be ordered and searched?
	# (since we're in a loop, need capture variables above by using defaulted arguments)
	def makeAStar(rcs, luteEcalc=luteEcalc, pmat=pmat):
		return osprey.LUTE_AStar(rcs, pmat, luteEcalc, showProgress=False)
	info.confSearchFactoryMinimized = osprey.BBKStar.ConfSearchFactory(makeAStar)
	info.confSearchFactoryRigid = info.confSearchFactoryMinimized

	# how should we score each sequence?
	# (since we're in a loop, need capture variables above by using defaulted arguments)
	def makePfunc(rcs, luteEcalc=luteEcalc, pmat=pmat):
		return osprey.LUTE_Pfunc(
			luteEcalc,
			osprey.LUTE_AStar(rcs, pmat, luteEcalc, showProgress=False),
			rcs
		)
	info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

# run BBK*
scoredSequences = bbkstar.run()

# make a sequence analyzer to look at the results
analyzer = osprey.SequenceAnalyzer(bbkstar)

# LUTE needs a regular energy calculator to analyze sequences
ffparams = osprey.ForcefieldParams()
confSpace = confSpaces['complex']
ecalc = osprey.EnergyCalculator(confSpace, ffparams)
eref = osprey.ReferenceEnergies(confSpace, ecalc)
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# use results
for scoredSequence in scoredSequences:
	print("result:")
	print("\tsequence: %s" % scoredSequence.sequence)
	print("\tK* score: %s" % scoredSequence.score)

	# write the sequence ensemble, with up to 10 of the lowest-energy conformations
	numConfs = 10
	analysis = analyzer.analyze(scoredSequence.sequence, numConfs, confEcalc)
	print(analysis)
	analysis.writePdb(
		'seq.%s.pdb' % scoredSequence.sequence,
		'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
	)
