import sys

import osprey
osprey.start()

# import the compiled conformation spaces module after starting Osprey
import osprey.ccs


# get CL args
complex_name = sys.argv[1]
target_name = sys.argv[2]
peptide_name = sys.argv[3]

# configure K*
kstar = osprey.KStar(

    # load the conformation spaces for your design
    # these conformation spaces were prepared in a previous step
    # see the 'conf space prep' example
    proteinConfSpace=osprey.ccs.loadConfSpace(target_name),
    ligandConfSpace=osprey.ccs.loadConfSpace(peptide_name),
    complexConfSpace=osprey.ccs.loadConfSpace(complex_name),

    # this is a good starting point for partiton function value precision,
    # but you may need something more precise (smaller epsilon)
    # or less precise (larger epsilon) in your designs
    epsilon=0.68,

    # save the sequence results somewhere so we can find them later
    writeSequencesToConsole=True,
    # writeSequencesToFile='sequences.tsv',

    # search for up to double mutants
    maxSimultaneousMutations=1000,

    # turn this one off when you get tired of the log spam
    showPfuncProgress=False
)

# more cores is more faster
parallelism = osprey.Parallelism(cpuCores=48)
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
counter = 1
print("results")
for scored_sequence in scored_sequences:
    print(" " + str(counter) + "," + str(scored_sequence.sequence) + "," + str(scored_sequence.score))
    counter += 1

    # write the sequence ensemble
    numConfs = 10
    analysis = analyzer.analyze(scored_sequence.sequence, numConfs)

    analysis.writePdb(
        'ensembles/seq.%s.pdb' % scored_sequence.sequence,
        'Top %d conformations for sequence %s' % (numConfs, scored_sequence.sequence)
    )