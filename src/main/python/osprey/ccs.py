
import sys

# throw an error for python 2
if sys.version_info[0] <= 2:
    raise Exception('Python v2 or earlier is not supported in this module')


import osprey
import osprey.jvm as jvm


_useJavaDefault = osprey.useJavaDefault


# export enums
Precision = osprey.c.gpu.Structs.Precision


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
