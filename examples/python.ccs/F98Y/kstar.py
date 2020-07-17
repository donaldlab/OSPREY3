
import osprey
osprey.start()

# load the conformation spaces
complex = osprey.readConfSpace("4tu5.complex.ccsx")
dhfr = osprey.readConfSpace("4tu5.DHFR.ccsx")
nadph06w = osprey.readConfSpace("4tu5.NADPH.06W.ccsx")

# more cores is more faster
parallelism = osprey.Parallelism(cpuCores=4)
tasks = parallelism.makeTaskExecutor()

# configure K*
kstar = osprey.KStar(
	dhfr,
	nadph06w,
	complex,
	epsilon=0.68,
	writeSequencesToConsole=True,

	# turn this one off when you get tired of the log spam
	showPfuncProgress=True
)

# configure K* inputs for each conf space
for info in kstar.confSpaceInfos():

	# TODO: make python API for constructor? using parallelism?
	ecalc = osprey.c.energy.compiled.CPUConfEnergyCalculator(info.confSpace)

	# compute reference energies
	# TODO: make python API
	eref = osprey.jvm.getInnerClass(osprey.c.ematrix.compiled.ErefCalculator, 'Builder')(ecalc) \
		.build() \
		.calc()

	# TODO: explain why this setting is important
	includeStaticStatic = True

	# TODO: explain what this does
	posInterDist = osprey.c.confspace.compiled.PosInterDist.DesmetEtAl1992
	#posInterDist = osprey.c.confspace.compiled.PosInterDist.TighterBounds

	# compute the energy matrix
	emat = osprey.jvm.getInnerClass(osprey.c.ematrix.compiled.EmatCalculator, 'Builder')(ecalc) \
		.setPosInterDist(posInterDist) \
		.setReferenceEnergies(eref) \
		.setIncludeStaticStatic(includeStaticStatic) \
		.build() \
		.calc()

	info.confEcalc = osprey.jvm.getInnerClass(osprey.c.energy.compiled.ConfEnergyCalculatorAdapter, 'Builder')(ecalc, tasks) \
		.setPosInterDist(posInterDist) \
		.setReferenceEnergies(eref) \
		.setIncludeStaticStatic(includeStaticStatic) \
		.build()

	# how should we score each sequence?
	# (since we're in a loop, need capture variables above by using defaulted arguments)
	def makePfunc(rcs, confEcalc=info.confEcalc, emat=emat):
		return osprey.PartitionFunction(
			confEcalc,
			osprey.AStarTraditional(emat, rcs, showProgress=False),
			osprey.AStarTraditional(emat, rcs, showProgress=False),
			rcs
		)
	info.pfuncFactory = osprey.KStar.PfuncFactory(makePfunc)

# run K*
scoredSequences = kstar.run(tasks)

# use results
analyzer = osprey.SequenceAnalyzer(kstar)
for scoredSequence in scoredSequences:

	print("result:")
	print("\tsequence: %s" % scoredSequence.sequence)
	print("\tscore: %s" % scoredSequence.score)

	# write the sequence ensemble
	numConfs = 10
	analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
	print(analysis)

	analysis.writePdb(
		'seq.%s.pdb' % scoredSequence.sequence,
		'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
	)
