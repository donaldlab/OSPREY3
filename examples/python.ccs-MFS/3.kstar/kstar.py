
import osprey
osprey.start()

# import the compiled conformation spaces module after starting Osprey
import osprey.ccs


# configure K*
kstar = osprey.KStar(

    # load the conformation spaces for your design
    # these conformation spaces were prepared in a previous step
    # see the 'conf space prep' example
    proteinConfSpace=osprey.ccs.loadConfSpace("../2.conf-space-prep/CALP-python.ccsx"),
    ligandConfSpace=osprey.ccs.loadConfSpace("../2.conf-space-prep/kCAL01-python.ccsx"),
    complexConfSpace=osprey.ccs.loadConfSpace("../2.conf-space-prep/complex-python.ccsx"),

    # this is a good starting point for partiton function value precision,
    # but you may need something more precise (smaller epsilon)
    # q* >= q(1 -e)
    # or less precise (larger epsilon) in your designs
    # recommend: 0.68
    epsilon=0.05,

    # save the sequence results somewhere so we can find them later
    writeSequencesToConsole=True,
    writeSequencesToFile='python.csv',

    # this doesn't matter for MFS, which has no AA identity changes
    maxSimultaneousMutations=1000,

    # turn this one off when you get tired of the log spam
    showPfuncProgress=True
)

# more cores is more faster
parallelism = osprey.Parallelism(cpuCores=4)
tasks = parallelism.makeTaskExecutor()

# configure K* inputs for each conf space
for info in kstar.confSpaceInfos():

    # make the energy calculator
    ecalc = osprey.ccs.bestEnergyCalculator(info.confSpace, parallelism)

    # compute reference energies (optional step)
    eref = osprey.ccs.calcReferenceEnergies(ecalc)

    # compute the energy matrix
    emat = osprey.ccs.calcEnergyMatrix(
        ecalc,
        tasks,
        eref,
        # # cache the energy matrices so we only need to compute them once,
        # # even if you decide to run the design again
        # # if the cached value doesn't match the current design settings, don't worry
        # # the energy matrices should get automatically recalculated
        # cachePath='emat.%s.dat' % info.id,
    )

    # adapt the new energy calculation system into the older K* implementation
    info.confEcalc = osprey.ccs.ecalcAdapter(ecalc, tasks, eref)

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
scored_sequences = kstar.run(tasks)

# use results
analyzer = osprey.SequenceAnalyzer(kstar)
for scored_sequence in scored_sequences:

    print("result:")
    print("\tsequence: %s" % scored_sequence.sequence)
    print("\tscore: %s" % scored_sequence.score)

    # write the sequence ensemble
    numConfs = 10
    analysis = analyzer.analyze(scored_sequence.sequence, numConfs)
    print(analysis)

    analysis.writePdb(
        'seq.%s.pdb' % scored_sequence.sequence,
        'Top %d conformations for sequence %s' % (numConfs, scored_sequence.sequence)
    )
