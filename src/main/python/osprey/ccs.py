
'''
Osprey for Compiled Conformation Spaces (CCS)

Provides functions to setup designs using the compiled conformation space system, and run them.
'''


import sys

# throw an error for python 2
if sys.version_info[0] <= 2:
    raise Exception('Python v2 or earlier is not supported in this module')


import osprey
import osprey.jvm as jvm


_useJavaDefault = osprey.useJavaDefault


# export enums

Precision = osprey.c.gpu.Structs.Precision
'''
${enum_java(.gpu.Structs$Precision)}
'''

PosInterDist = osprey.c.confspace.compiled.PosInterDist
'''
${enum_java(.confspace.compiled.PosInterDist)}
'''


# fix Hazelcast's completely obnoxious default logging settings
osprey.c.parallelism.Cluster.fixHazelcastLogging()


def loadConfSpace(path):
    '''
    ${method_javadoc(.confspace.compiled.ConfSpace#fromBytes)}

    # Arguments
    path `str`: Path to the compiled conformation space file (usually has a .ccsx or .ccs extension)

    # Returns
    ${returns_method_java(.confspace.compiled.ConfSpace#fromBytes)}
    '''
    return osprey.c.confspace.compiled.ConfSpace.fromBytes(osprey.c.tools.FileTools.readFileBytes(path))


def cudaEnergyCalculator(confSpace, precision, parallelism):
    '''
    ${class_javadoc(.energy.compiled.CudaConfEnergyCalculator)}
    ${method_javadoc(.energy.compiled.CudaConfEnergyCalculator#<init>(ConfSpace,Precision,Parallelism))}

    # Arguments
    ${args_java(.energy.compiled.CudaConfEnergyCalculator#<init>(ConfSpace,Precision,Parallelism),
        [confSpace],
        [precision],
        [parallelism]
    )}

    # Returns
    ${type_java(.energy.compiled.CudaConfEnergyCalculator)}
    '''
    return osprey.c.energy.compiled.CudaConfEnergyCalculator(confSpace, precision, parallelism)

def nativeEnergyCalculator(confSpace, precision):
    '''
    ${class_javadoc(.energy.compiled.NativeConfEnergyCalculator)}

    # Arguments
    ${args_java(.energy.compiled.NativeConfEnergyCalculator#<init>(ConfSpace,Precision),
        [confSpace],
        [precision]
    )}

    # Returns
    ${type_java(.energy.compiled.NativeConfEnergyCalculator)}
    '''
    return osprey.c.energy.compiled.NativeConfEnergyCalculator(confSpace, precision)

def javaEnergyCalculator(confSpace):
    '''
    ${class_javadoc(.energy.compiled.CPUConfEnergyCalculator)}

    # Arguments
    ${args_java(.energy.compiled.CPUConfEnergyCalculator#<init>,
        [confSpace]
    )}

    # Returns
    ${type_java(.energy.compiled.CPUConfEnergyCalculator)}
    '''
    return osprey.c.energy.compiled.CPUConfEnergyCalculator(confSpace)

def bestEnergyCalculator(confSpace, parallelism):
    '''
    ${method_javadoc(.energy.compiled.ConfEnergyCalculator#makeBest(ConfSpace,Parallelism)ConfEnergyCalculator)}

    # Arguments
    ${args_java(.energy.compiled.ConfEnergyCalculator#makeBest(ConfSpace,Parallelism)ConfEnergyCalculator,
        [confSpace],
        [parallelism]
    )}

    # Returns
    ${type_java(.energy.compiled.ConfEnergyCalculator)}
    '''
    return osprey.c.energy.compiled.ConfEnergyCalculator.makeBest(confSpace, parallelism)


def calcReferenceEnergies(ecalc, minimize=_useJavaDefault):
    '''
    ${class_javadoc(.ematrix.compiled.ErefCalculator)}

    # Arguments
    ${args_java(.ematrix.compiled.ErefCalculator$Builder#<init>,
        [ecalc, confEcalc]
    )}
    ${args_fields_javadoc(.ematrix.compiled.ErefCalculator$Builder,
        [minimize]
    )}

    # Returns
    ${type_java(.ematrix.SimpleReferenceEnergies)}
    '''

    builder = osprey.c.ematrix.compiled.ErefCalculator.Builder(ecalc)

    if minimize is not _useJavaDefault:
        builder.setMinimize(minimize)

    return builder.build().calc()


def calcEnergyMatrix(ecalc, tasks=_useJavaDefault, eref=_useJavaDefault, posInterDist=_useJavaDefault, minimize=_useJavaDefault, includeStaticStatic=_useJavaDefault, cachePath=_useJavaDefault):
    '''
    ${class_javadoc(.ematrix.compiled.EmatCalculator)}

    # Arguments
    ${args_java(.ematrix.compiled.EmatCalculator$Builder#<init>,
        [ecalc, confEcalc]
    )}
    ${args_fields_javadoc(.ematrix.compiled.EmatCalculator$Builder,
        [eref],
        [posInterDist],
        [minimize],
        [includeStaticStatic],
        [cachePath, cacheFile, type=str]
    )}

    # Returns
    ${type_java(.ematrix.EnergyMatrix)}
    '''

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
    '''
    ${class_javadoc(.energy.compiled.ConfEnergyCalculatorAdapter)}

    # Arguments
    ${args_java(.energy.compiled.ConfEnergyCalculatorAdapter$Builder#<init>,
        [ecalc, confEcalc],
        [tasks]
    )}
    ${args_fields_javadoc(.energy.compiled.ConfEnergyCalculatorAdapter$Builder,
        [eref],
        [posInterDist],
        [minimize],
        [includeStaticStatic]
    )}

    # Returns
    ${returns_method_java(.energy.compiled.ConfEnergyCalculatorAdapter$Builder#build)}
    '''

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
    '''
    A free energy calculator based on a scalable and memory-bounded partition function calculator.

    Used with memory-bounded implementations of K*. See #kstarBoundedMem

    # Arguments
    ${args_java(.coffee.Coffee$Builder#<init>,
        [confSpace]
    )}
    ${args_fields_javadoc(.coffee.Coffee$Builder,
        [parallelism],
        [cluster],
        [precision],
        [nodeDBFile, nodedbFile, type=str],
        [nodeDBMem, nodedbMemBytes],
        [seqDBFile, seqdbFile, type=str],
        [seqDBMathContext, seqdbMathContext]
    )}
    ${arg_java(posInterDist, .energy.compiled.PosInterGen#<init>, dist)}
    ${args_fields_javadoc(.coffee.Coffee$Builder,
        [staticStatic, includeStaticStatic],
        [tripleCorrectionThreshold],
        [conditions],
        [nodeScoringLog],
        [nodeStatsReportingSeconds, nodeStatsReportingInterval, type=int]
    )}

    # Returns
    ${type_java(.coffee.Coffee)}
    '''

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
    '''
    An implementation of the K* design algorithm that uses the bounded-memory free energy calculator. See #freeEnergyCalc.

    For an example of how to use this function, see `examples/python.ccs/kstar/kstar.boundedMem.py` in your Osprey distribution.

    # Arguments
    ${args_java(.coffee.directors.KStarDirector$Builder#<init>(ConfSpace,ConfSpace,ConfSpace),
        [complex],
        [design],
        [target]
    )}
    ${args_fields_javadoc(.coffee.directors.KStarDirector$Builder,
        [gWidthMax],
        [maxSimultaneousMutations],
        [stabilityThreshold],
        [timing],
        [reportStateProgress]
    )}
    ensembleTracking `[int,str]`: ${method_javadoc(.coffee.directors.KStarDirector$Builder#setEnsembleTracking)}
    ensembleMinUpdate `[int,java.util.concurrent.TimeUnit]`: ${method_javadoc(.coffee.directors.KStarDirector$Builder#setEnsembleMinUpdate)}

    # Returns
    ${returns_method_java(.coffee.directors.KStarDirector$Builder#build)}
    '''

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
