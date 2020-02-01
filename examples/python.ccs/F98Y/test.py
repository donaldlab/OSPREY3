
import osprey
osprey.start()

# load the conformation spaces
complex = osprey.readConfSpace("4tu5.complex.ccs.xz")
dhfr = osprey.readConfSpace("4tu5.DHFR.ccs.xz")
nadph06w = osprey.readConfSpace("4tu5.NADPH.06W.ccs.xz")

# more cores is more faster
parallelism = osprey.Parallelism(cpuCores=4)
tasks = parallelism.makeTaskExecutor()

# configure K*
kstar = osprey.KStar(
	dhfr,
	nadph06w,
	complex,
	epsilon=0.68,
	writeSequencesToConsole=True
)

# configure K* inputs for each conf space
for info in kstar.confSpaceInfos():

	# TODO: make python API for constructor? using parallelism?
	ecalc = osprey.c.energy.compiled.CPUConfEnergyCalculator(info.confSpace, tasks)

	# compute reference energies
	# TODO: make python API
	eref = osprey.jvm.getInnerClass(osprey.c.ematrix.compiled.ErefCalculator, 'Builder')(ecalc) \
		.build() \
		.calc()

	# TODO: explain why this setting is important
	includeStaticStatic = True

	# TODO: explain what this does
	posInterDist = osprey.c.confspace.compiled.PosInterDist.DesmetEtAl1992

	# compute the energy matrix
	emat = osprey.jvm.getInnerClass(osprey.c.ematrix.compiled.EmatCalculator, 'Builder')(ecalc) \
		.setPosInterDist(posInterDist) \
		.setReferenceEnergies(eref) \
		.setIncludeStaticStatic(includeStaticStatic) \
		.build() \
		.calc()

	info.confEcalc = osprey.jvm.getInnerClass(osprey.c.energy.compiled.ConfEnergyCalculatorAdapter, 'Builder')(ecalc) \
		.setPosInterDist(posInterDist) \
		.setReferenceEnergies(eref) \
		.setIncludeStaticStatic(includeStaticStatic) \
		.build()

	# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
	def makeAStar(rcs, emat=emat):
		return osprey.AStarTraditional(emat, rcs, showProgress=False)
	info.confSearchFactory = osprey.KStar.ConfSearchFactory(makeAStar)

	# use a conformation database, to make computing ensembles easier
	info.setConfDBFile('kstar.%s.db' % info.type.name().lower())

# run K*
scoredSequences = kstar.run()

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
