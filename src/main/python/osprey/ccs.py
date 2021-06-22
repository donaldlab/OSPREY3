
import sys

# throw an error for python 2
if sys.version_info[0] <= 2:
    raise Exception('Python v2 or earlier is not supported in this module')


import osprey
import osprey.jvm as jvm


_useJavaDefault = osprey.useJavaDefault


# export enums
Precision = osprey.c.gpu.Structs.Precision
PosInterDist = osprey.c.confspace.compiled.PosInterDist


# fix Hazelcast's completely obnoxious default logging settings
osprey.c.parallelism.Cluster.fixHazelcastLogging()


def loadConfSpace(path):
    return osprey.c.confspace.compiled.ConfSpace.fromBytes(osprey.c.tools.FileTools.readFileBytes(path))


def cudaEnergyCalculator(confSpace, precision, parallelism):
    return osprey.c.energy.compiled.CudaConfEnergyCalculator(confSpace, precision, parallelism)

def nativeEnergyCalculator(confSpace, precision):
    return osprey.c.energy.compiled.NativeConfEnergyCalculator(confSpace, precision)

def javaEnergyCalculator(confSpace, precision):
    return osprey.c.energy.compiled.CPUConfEnergyCalculator(confSpace, precision)

def bestEnergyCalculator(confSpace, parallelism):
    return osprey.c.energy.compiled.ConfEnergyCalculator.makeBest(confSpace, parallelism)


def calcReferenceEnergies(ecalc, minimize=_useJavaDefault):

    builder = osprey.c.ematrix.compiled.ErefCalculator.Builder(ecalc)

    if minimize is not _useJavaDefault:
        builder.setMinimize(minimize)

    return builder.build().calc()


def calcEnergyMatrix(ecalc, tasks=_useJavaDefault, eref=_useJavaDefault, posInterDist=_useJavaDefault, minimize=_useJavaDefault, includeStaticStatic=_useJavaDefault, cachePath=_useJavaDefault):

    builder = osprey.c.ematrix.compiled.EmatCalculator.Builder(ecalc)

    if eref is not _useJavaDefault:
        builder.setReferenceEnergies(eref)
    if posInterDist is not _useJavaDefault:
        builder.setPosInterDist(posInterDist)
    if minimize is not _useJavaDefault:
        builder.setMinimize(minimize)
    if includeStaticStatic is not _useJavaDefault:
        builder.setIncludeStaticStatic(includeStaticStatic)
    if cachePath is not _useJavaDefault:
        builder.setCacheFile(jvm.toFile(cachePath))

    calc = builder.build()

    if tasks is not _useJavaDefault:
        return calc.calc(tasks)
    else:
        return calc.calc()


def ecalcAdapter(ecalc, tasks, eref=_useJavaDefault, posInterDist=_useJavaDefault, minimize=_useJavaDefault, includeStaticStatic=_useJavaDefault):

    builder = osprey.c.energy.compiled.ConfEnergyCalculatorAdapter.Builder(ecalc, tasks)

    if eref is not _useJavaDefault:
        builder.setReferenceEnergies(eref)
    if posInterDist is not _useJavaDefault:
        builder.setPosInterDist(posInterDist)
    if minimize is not _useJavaDefault:
        builder.setMinimize(minimize)
    if includeStaticStatic is not _useJavaDefault:
        builder.setIncludeStaticStatic(includeStaticStatic)

    return builder.build()



def freeEnergyCalc(
    confSpace,
    parallelism=_useJavaDefault,
    cluster=_useJavaDefault,
    precision=_useJavaDefault,
    nodeDBFile=_useJavaDefault,
    nodeDBMem=_useJavaDefault,
    seqDBFile=_useJavaDefault,
    seqDBMathContext=_useJavaDefault,
    posInterDist=PosInterDist.TighterBounds,
    staticStatic=_useJavaDefault,
    tripleCorrectionThreshold=_useJavaDefault,
    conditions=_useJavaDefault,
    nodeScoringLog=_useJavaDefault,
    nodeStatsReportingSeconds=_useJavaDefault
):

    builder = osprey.c.coffee.Coffee.Builder(confSpace)

    # config each state using the position interaction distribution, and no reference energies
    # TODO: expose API to calculate reference energies?
    for state in confSpace.states:
        config = osprey.c.coffee.Coffee.StateConfig(state)
        eref = None
        config.posInterGen = osprey.c.energy.compiled.PosInterGen(posInterDist, eref)
        builder.configState(config)

    if parallelism is not _useJavaDefault:
        builder.setParallelism(parallelism)
    if cluster is not _useJavaDefault:
        builder.setCluster(cluster)
    if precision is not _useJavaDefault:
        builder.setPrecision(precision)
    if nodeDBFile is not _useJavaDefault:
        builder.setNodeDBFile(jvm.toFile(nodeDBFile[0]), nodeDBFile[1])
    if nodeDBMem is not _useJavaDefault:
        builder.setNodeDBMem(nodeDBMem)
    if seqDBFile is not _useJavaDefault:
        builder.setSeqDBFile(jvm.toFile(seqDBFile))
    if seqDBMathContext is not _useJavaDefault:
        builder.setSeqDBMathContext(seqDBMathContext)
    if staticStatic is not _useJavaDefault:
        builder.setStaticStatic(staticStatic)
    if tripleCorrectionThreshold is not _useJavaDefault:
        builder.setTripleCorrectionThreshold(tripleCorrectionThreshold)
    if conditions is not _useJavaDefault:
        builder.setConditions(conditions)
    if nodeScoringLog is not _useJavaDefault:
        builder.setNodeScoringLog(jvm.toFile(nodeScoringLog))
    if nodeStatsReportingSeconds is not _useJavaDefault:
        builder.setNodeStatsReportingInterval(jvm.c.java.time.Duration.ofSeconds(nodeStatsReportingSeconds))

    return builder.build()


def kstarBoundedMem(
    complex,
    design,
    target,
    gWidthMax=_useJavaDefault,
    maxSimultaneousMutations=_useJavaDefault,
    stabilityThreshold=_useJavaDefault,
    timing=_useJavaDefault,
    reportStateProgress=_useJavaDefault,
    ensembleTracking=_useJavaDefault,
    ensembleMinUpdate=_useJavaDefault
):

    builder = osprey.c.coffee.directors.KStarDirector.Builder(complex, design, target)

    if gWidthMax is not _useJavaDefault:
        builder.setGWidthMax(gWidthMax)
    if maxSimultaneousMutations is not _useJavaDefault:
        builder.setMaxSimultaneousMutations(jvm.boxInt(maxSimultaneousMutations))
    if stabilityThreshold is not _useJavaDefault:
        builder.setStabilityThreshold(stabilityThreshold)
    if timing is not _useJavaDefault:
        builder.setTiming(timing)
    if reportStateProgress is not _useJavaDefault:
        builder.setReportStateProgress(reportStateProgress)
    if ensembleTracking is not _useJavaDefault:
        builder.setEnsembleTracking(ensembleTracking[0], jvm.toFile(ensembleTracking[1]))
    if ensembleMinUpdate is not _useJavaDefault:
        builder.setEnsembleMinUpdate(ensembleMinUpdate[0], ensembleMinUpdate[1])

    return builder.build()
