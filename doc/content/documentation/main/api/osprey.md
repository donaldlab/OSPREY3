+++
title = "osprey"
weight = 5
+++

<a name="osprey.WILD_TYPE"></a>
#### WILD\_TYPE

Magic constant that refers to the wild-type template of the residue

<a name="osprey.start_with_jvm_args"></a>
#### start\_with\_jvm\_args

```python
start_with_jvm_args(jvm_args=None)
```

Start the jvm with JVM arguments.

**Arguments**:

- `jvm_args`: a list of arguments to the JVM. For example, ["-XX:MaxHeapSize=100g", "-XX:MinHeapSize=100g"].
defaults to None, which implies the JVM will run with its default configuration.

**Returns**:

None

<a name="osprey.start"></a>
#### start

```python
start(heapSizeMiB=1024, enableAssertions=False, stackSizeMiB=16, garbageSizeMiB=128, allowRemoteManagement=False, attachJvmDebugger=False)
```

Starts the Java Virtual Machine (JVM) that runs Osprey's computation libraries.

Call :meth:`start` before using any of Osprey's other functions.

:param int heapSizeMiB: Size of the JVM heap in megabytes. This is essentially the amount of memory
Osprey will have to do computations. 1024 MiB is 1 GiB, but for larger designs,
you may want to use 2048 MiB (2 GiB), 4096 MiB (4 GiB), or even more memory.

:param bool enableAssertions: pass ``True`` to enable JVM assertions. Only useful for debugging.

:param int stackSizeMiB: Size of the JVM stack portion of the heap in megabytes.
Generally leave this at the default value, unless Osprey crashes because it's too small. Then try increasing it.

:param int garbageSizeMiB: Size of the garbage portion of the JVM heap that is reserved for temporary objects.
This default value is appropriate for the default heap size, but if using larger heap sizes, then increasing
the garbage size to 256, 512, or even 1024 MiB can give a modest improvement in performance.

:param bool allowRemoteManagement: pass ``True`` to listen on ports 9010 and 9011 for JMX remote management.

:param bool attachJvmDebugger: pass ``True`` to be able to attach a Java debugger to the JVM.

<a name="osprey.readTextFile"></a>
#### readTextFile

```python
readTextFile(path)
```

Useful for passing file contents to functions that expect strings rather than file paths

**Returns**:

the text of the file at ``path``
:rtype: str

<a name="osprey.readPdb"></a>
#### readPdb

```python
readPdb(path)
```

loads a PDB file into a molecule object

.. note:: Raw molecules cannot be used directly in designs.
Create a :py:meth:`Strand` with the molecule to use it in a design.

:param str path: path to PDB file
:rtype: :java:ref:`.structure.Molecule`

<a name="osprey.writePdb"></a>
#### writePdb

```python
writePdb(mol, path, comment=None)
```

save a molecule to a PDB file

**Arguments**:

- `mol`: the molecule to save
:type mol: :java:ref:`.structure.Molecule`
:param str path: path of the PDB file
:param str comment: Optional comment to add to the PDB file headers

<a name="osprey.printGpuInfo"></a>
#### printGpuInfo

```python
printGpuInfo()
```

Prints information about GPU hardware present in the system
and which GPUs are usable by Osprey.

Both Cuda and OpenCL APIs are queried.

<a name="osprey.initExternalMemory"></a>
#### initExternalMemory

```python
initExternalMemory(internalSizeMiB, tempDir=None, tempSubdir=None)
```

Initializes external memory for calculations.

:param int internalSizeMiB: :java:methoddoc:`.externalMemory.ExternalMemory#setInternalLimit`
:param str tempDir: Path to temporary directory to host external memory
:default tempDir: <system temp dir>
:param str tempSubdir: name of subdirectory within tempDir
:default tempSubdir: <automatically generated>

<a name="osprey.Parallelism"></a>
#### Parallelism

```python
Parallelism(cpuCores=None, gpus=None, streamsPerGpu=None)
```

:java:classdoc:`.parallelism.Parallelism`

:builder_option cpuCores .parallelism.Parallelism$Builder#numCpus:
:builder_option gpus .parallelism.Parallelism$Builder#numGpus:
:builder_option streamsPerGpu .parallelism.Parallelism$Builder#numStreamsPerGpu:
:builder_return .parallelism.Parallelism$Builder:

<a name="osprey.TemplateLibrary"></a>
#### TemplateLibrary

```python
TemplateLibrary(ffparams=None, defaultTemplates=True, extraTemplates=[], defaultTemplateCoords=True, extraTemplateCoords=[], defaultRotamers=True, extraRotamers=[], extraBackboneDependentRotamers=[], defaultResidueEntropies=True, extraResidueEntropies=[], makeDAminoAcids=None, moleculesForWildTypeRotamers=[])
```

:java:classdoc:`.restypes.ResidueTemplateLibrary`

:builder_option ffparams .restypes.ResidueTemplateLibrary$Builder#ffparams:

:param bool defaultTemplates: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#clearTemplates`

**Arguments**:

- `extraTemplates`: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addTemplates`
:type extraTemplates: [template string or file path]

:param bool defaultTemplateCoords: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#clearTemplateCoords`
- `extraTemplateCoords`: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addTemplateCoords`
:type extraTemplateCoords: [coords string or file path]

:param bool defaultRotamers: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#clearRotamers`
- `extraRotamers`: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addRotamers`
:type extraRotamers: [rotamers string or file path]

- `extraBackboneDependentRotamers`: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addBackboneDependentRotamers`
:type extraBackboneDependentRotamers: [backbone-dependent rotamers string or file path]

:builder_option makeDAminoAcids .restypes.ResidueTemplateLibrary$Builder#makeDAminoAcidTemplates:

- `moleculesForWildTypeRotamers`: :java:methoddoc:`.restypes.ResidueTemplateLibrary$Builder#addMoleculeForWildTypeRotamers`
:type moleculesForWildTypeRotamers: [:java:ref:`.structure.Molecule`]

:builder_return .restypes.ResidueTemplateLibrary$Builder:

<a name="osprey.Strand"></a>
#### Strand

```python
Strand(pathOrMol, residues=None, templateLib=None, templateMatchingMethod=None)
```

:java:classdoc:`.confspace.Strand`

**Arguments**:

- `pathOrMol`: path to a PDB file, or a molecule instance
:type pathOrMol: str or :java:ref:`.structure.Molecule`

- `residues`: range of residue numbers, inclusive. `None` to include all residues.
:type residues: [str, str]

:builder_option templateLib .confspace.Strand$Builder#templateLib:
:builder_option templateMatchingMethod .confspace.Strand$Builder#templateMatchingMethod:

:builder_return .confspace.Strand$Builder:

<a name="osprey.ConfSpace"></a>
#### ConfSpace

```python
ConfSpace(strands, shellDist=None)
```

:java:classdoc:`.confspace.SimpleConfSpace`

**Arguments**:

- `strands`: the strands to use
:type strands: :java:ref:`.confspace.Strand` or list of Strands
:builder_option shellDist .confspace.SimpleConfSpace$Builder#shellDist:
:builder_return .confspace.SimpleConfSpace$Builder:

<a name="osprey.StateMutable"></a>
## StateMutable Objects

```python
class StateMutable(object)
```

<a name="osprey.StateMutable.__init__"></a>
#### \_\_init\_\_

```python
 | __init__(name, confSpace)
```

:param str name: unique for the state

**Arguments**:

- `confSpace`: the conformation space
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

<a name="osprey.StateUnmutable"></a>
## StateUnmutable Objects

```python
class StateUnmutable(object)
```

<a name="osprey.StateUnmutable.__init__"></a>
#### \_\_init\_\_

```python
 | __init__(name, confSpace)
```

:param str name: unique for the state

**Arguments**:

- `confSpace`: the conformation space
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

<a name="osprey.MultiStateConfSpace"></a>
#### MultiStateConfSpace

```python
MultiStateConfSpace(states)
```

:java:classdoc:`.confspace.MultiStateConfSpace`

**Arguments**:

- `states`: A list of states to use
:type states: list of StateMutable and StateUnmutable instances

<a name="osprey.ForcefieldParams"></a>
#### ForcefieldParams

```python
ForcefieldParams(params=None)
```

:java:classdoc:`.energy.forcefield.ForcefieldParams`

Configure the forcefield parameters by setting the properties of the :java:ref:`.energy.forcefield.ForcefieldParams` object.

:param str params: text of forcefield parameters in Amber format, to override defaults if desired
:rtype: :java:ref:`.energy.forcefield.ForcefieldParams`

<a name="osprey.EnergyCalculator"></a>
#### EnergyCalculator

```python
EnergyCalculator(confSpace, ffparams, parallelism=None, type=None, isMinimizing=None, infiniteWellEnergy=None)
```

:java:classdoc:`.energy.EnergyCalculator`

**Arguments**:

- `confSpace`: The conformation space containing the residue templates to use for atom connectivities.
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
:builder_option ffparams .energy.EnergyCalculator$Builder#ffparams:
:builder_option parallelism .energy.EnergyCalculator$Builder#parallelism:
:builder_option type .energy.EnergyCalculator$Builder#type:
:builder_option isMinimizing .energy.EnergyCalculator$Builder#isMinimizing:
:builder_option infiniteWellEnergy .energy.EnergyCalculator$Builder#infiniteWellEnergy:

:builder_return .energy.EnergyCalculator$Builder:

<a name="osprey.SharedEnergyCalculator"></a>
#### SharedEnergyCalculator

```python
SharedEnergyCalculator(ecalc, isMinimizing=None)
```

:java:classdoc:`.energy.EnergyCalculator$SharedBuilder`

**Arguments**:

- `ecalc`: The existing energy calculator with which to share resources
:type ecalc: :java:ref:`.energy.EnergyCalculator`
:builder_option isMinimizing .energy.EnergyCalculator$Builder#isMinimizing:
:builder_return .energy.EnergyCalculator$SharedBuilder:

<a name="osprey.ConfEnergyCalculator"></a>
#### ConfEnergyCalculator

```python
ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=UseJavaDefault, addResEntropy=UseJavaDefault, energyPartition=UseJavaDefault, amat=UseJavaDefault, approximationErrorBudget=UseJavaDefault)
```

:java:classdoc:`.energy.ConfEnergyCalculator`

:builder_option confSpace .energy.ConfEnergyCalculator$Builder#confSpace:
:builder_option ecalc .energy.ConfEnergyCalculator$Builder#ecalc:
:builder_option referenceEnergies .energy.ConfEnergyCalculator$Builder#eref:
:builder_option addResEntropy .energy.ConfEnergyCalculator$Builder#addResEntropy:
:builder_option energyPartition .energy.ConfEnergyCalculator$Builder#epart:
:builder_option amat .energy.ConfEnergyCalculator$Builder#amat:
:builder_option approximationErrorBudget .energy.ConfEnergyCalculator$Builder#approximationErrorBudget:
:builder_return .energy.ConfEnergyCalculator$Builder:

<a name="osprey.ConfEnergyCalculatorCopy"></a>
#### ConfEnergyCalculatorCopy

```python
ConfEnergyCalculatorCopy(source, ecalc)
```

:java:classdoc:`.energy.ConfEnergyCalculator`

**Arguments**:

- `source`: The conformation energy calculator you wish to copy.
:type source: :java:ref:`.energy.ConfEnergyCalculator`
:builder_option ecalc .energy.ConfEnergyCalculator$Builder#ecalc:
:builder_return .energy.ConfEnergyCalculator$Builder:

<a name="osprey.ApproximatorMatrix"></a>
#### ApproximatorMatrix

```python
ApproximatorMatrix(confEcalc, cacheFile=UseJavaDefault, numSamplesPerParam=UseJavaDefault)
```

:java:classdoc:`.energy.approximation.ApproximatorMatrix`

:builder_option confEcalc .energy.approximation.ApproximatorMatrixCalculator#confEcalc:
:builder_option cacheFile .energy.approximation.ApproximatorMatrixCalculator#cacheFile:
:builder_option numSamplesPerParam .energy.approximation.ApproximatorMatrixCalculator#numSamplesPerParam:

<a name="osprey.EnergyMatrix"></a>
#### EnergyMatrix

```python
EnergyMatrix(confEcalc, cacheFile=UseJavaDefault, tripleCorrectionThreshold=UseJavaDefault, quadCorrectionThreshold=UseJavaDefault)
```

:java:methoddoc:`.ematrix.SimplerEnergyMatrixCalculator#calcEnergyMatrix`

:builder_option confEcalc .ematrix.SimplerEnergyMatrixCalculator$Builder#confEcalc:
:builder_option cacheFile .ematrix.SimplerEnergyMatrixCalculator$Builder#cacheFile:
:builder_option tripleCorrectionThreshold .ematrix.SimplerEnergyMatrixCalculator$Builder#tripleCorrectionThreshold:
:builder_option quadCorrectionThreshold .ematrix.SimplerEnergyMatrixCalculator$Builder#quadCorrectionThreshold:

<a name="osprey.ReferenceEnergies"></a>
#### ReferenceEnergies

```python
ReferenceEnergies(confSpace, ecalc, addResEntropy=None)
```

:java:methoddoc:`.ematrix.SimplerEnergyMatrixCalculator#calcReferenceEnergies`

:builder_option confSpace .ematrix.SimpleReferenceEnergies$Builder#confSpace:
:builder_option ecalc .ematrix.SimpleReferenceEnergies$Builder#ecalc:
:builder_option addResEntropy .ematrix.SimpleReferenceEnergies$Builder#addResEntropy:
:builder_return .ematrix.SimpleReferenceEnergies$Builder:

<a name="osprey.DEE"></a>
#### DEE

```python
DEE(confSpace, emat, singlesThreshold=useJavaDefault, pairsThreshold=useJavaDefault, singlesGoldsteinDiffThreshold=useJavaDefault, pairsGoldsteinDiffThreshold=useJavaDefault, triplesGoldsteinDiffThreshold=useJavaDefault, typeDependent=useJavaDefault, numIterations=useJavaDefault, singlesPlugThreshold=useJavaDefault, pairsPlugThreshold=useJavaDefault, triplesPlugThreshold=useJavaDefault, singlesTransitivePruning=useJavaDefault, pairsTransitivePruning=useJavaDefault, triplesTransitivePruning=useJavaDefault, showProgress=useJavaDefault, parallelism=useJavaDefault, cacheFile=useJavaDefault)
```

:java:classdoc:`.pruning.SimpleDEE$Runner`

**Arguments**:

- `confSpace`: The design conformation space
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
- `emat`: An energy matrix computed for the conformation space
:type emat: :java:ref:`.ematrix.EnergyMatrix`

:builder_option singlesThreshold .pruning.SimpleDEE$Runner#singlesThreshold:
:builder_option pairsThreshold .pruning.SimpleDEE$Runner#pairsThreshold:
:builder_option singlesGoldsteinDiffThreshold .pruning.SimpleDEE$Runner#singlesGoldsteinDiffThreshold:
:builder_option pairsGoldsteinDiffThreshold .pruning.SimpleDEE$Runner#pairsGoldsteinDiffThreshold:
:builder_option triplesGoldsteinDiffThreshold .pruning.SimpleDEE$Runner#triplesGoldsteinDiffThreshold:
:builder_option singlesPlugThreshold .pruning.SimpleDEE$Runner#singlesPlugThreshold:
:builder_option pairsPlugThreshold .pruning.SimpleDEE$Runner#pairsPlugThreshold:
:builder_option triplesPlugThreshold .pruning.SimpleDEE$Runner#triplesPlugThreshold:
:builder_option singlesTransitivePruning .pruning.SimpleDEE$Runner#singlesTransitivePruning:
:builder_option pairsTransitivePruning .pruning.SimpleDEE$Runner#pairsTransitivePruning:
:builder_option triplesTransitivePruning .pruning.SimpleDEE$Runner#triplesTransitivePruning:
:builder_option typeDependent .pruning.SimpleDEE$Runner#typeDependent:
:builder_option numIterations .pruning.SimpleDEE$Runner#numIterations:
:builder_option showProgress .pruning.SimpleDEE$Runner#showProgress:

<a name="osprey.DEE_read"></a>
#### DEE\_read

```python
DEE_read(confSpace, path)
```

Reads a pruning matrix from a file

**Arguments**:

- `confSpace`: The design conformation space
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
:param str path: Path to the file

<a name="osprey.AStarTraditional"></a>
#### AStarTraditional

```python
AStarTraditional(emat, confSpaceOrPmat, showProgress=True, useExternalMemory=False, maxNumNodes=useJavaDefault)
```

:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#setTraditional`

:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:

**Arguments**:

- `confSpaceOrPmat`: The conformation space containing the residue conformations to search.
:type confSpaceOrPmat: :java:ref:`.confspace.SimpleConfSpace` or :java:ref:`.pruning.PruningMatrix`
- `useExternalMemory`: set to True to use external memory.

:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#useExternalMemory`

:type useExternalMemory: boolean
:builder_option maxNumNodes .astar.conf.ConfAStarTree$Builder#maxNumNodes:
:builder_return .astar.conf.ConfAStarTree$Builder:

<a name="osprey.AStarMPLP"></a>
#### AStarMPLP

```python
AStarMPLP(emat, confSpaceOrPmat, updater=None, numIterations=None, convergenceThreshold=None, useExternalMemory=False, maxNumNodes=useJavaDefault)
```

:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#setMPLP`

:builder_option emat .astar.conf.ConfAStarTree$Builder#emat:

**Arguments**:

- `confSpaceOrPmat`: The conformation space containing the residue conformations to search.
:type confSpaceOrPmat: :java:ref:`.confspace.SimpleConfSpace` or :java:ref:`.pruning.PruningMatrix`
:builder_option updater .astar.conf.ConfAStarTree$MPLPBuilder#updater:
:builder_option numIterations .astar.conf.ConfAStarTree$MPLPBuilder#numIterations:
:builder_option convergenceThreshold .astar.conf.ConfAStarTree$MPLPBuilder#convergenceThreshold:
- `useExternalMemory`: set to True to use external memory.

:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#useExternalMemory`

:type useExternalMemory: boolean
:builder_option maxNumNodes .astar.conf.ConfAStarTree$Builder#maxNumNodes:
:builder_return .astar.conf.ConfAStarTree$Builder:

<a name="osprey.GMECFinder"></a>
#### GMECFinder

```python
GMECFinder(astar, confEcalc, confLog=None, printIntermediateConfs=None, useExternalMemory=None, resumeLog=None, confDBFile=None)
```

:java:classdoc:`.gmec.SimpleGMECFinder`

:builder_option astar .gmec.SimpleGMECFinder$Builder#search:

Use one of :py:func:`AStarTraditional` or :py:func:`AStarMPLP` to get an A* implementation.

:builder_option confEcalc .gmec.SimpleGMECFinder$Builder#confEcalc:

Use :py:func:`ConfEnergyCalculator` to get a conformation energy calculator.

:param str confLog: Path to file where conformations found during conformation space search should be logged.
:builder_option printIntermediateConfs .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole:
:builder_option useExternalMemory .gmec.SimpleGMECFinder$Builder#useExternalMemory:
:param str resumeLog: Path to log file where resume info will be written or read, so designs can be resumed.
:builder_return .gmec.SimpleGMECFinder$Builder:

<a name="osprey.DEEGMECFinder"></a>
#### DEEGMECFinder

```python
DEEGMECFinder(emat, confSpace, ecalc, confEcalc, name, use_epic, use_lute, confLog=None, printIntermediateConfs=None, useExternalMemory=None)
```

:java:classdoc:`.gmec.SimpleGMECFinder`

:builder_option astar .gmec.SimpleGMECFinder$Builder#search:

Use one of :py:func:`AStarTraditional` or :py:func:`AStarMPLP` to get an A* implementation.

:builder_option confEcalc .gmec.SimpleGMECFinder$Builder#confEcalc:

Use :py:func:`ConfEnergyCalculator` to get a conformation energy calculator.

:param str confLog: Path to file where conformations found during conformation space search should be logged.
:builder_option printIntermediateConfs .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole:
:builder_option useExternalMemory .gmec.SimpleGMECFinder$Builder#useExternalMemory:
:builder_return .gmec.SimpleGMECFinder$Builder:

<a name="osprey.Paste"></a>
#### Paste

```python
Paste(complexConfSpace, numPDBs, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, useWindowCriterion=useJavaDefault, maxNumPfConfs=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None, useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault, mutFile=useJavaDefault)
```

:java:classdoc:`.paste.Paste`

For examples using PAStE, see the examples/python.Paste directory in your Osprey distribution.

**Arguments**:

- `complexConfSpace`: :java:fielddoc:`.paste.Paste#protein`
:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
:builder_option epsilon .paste.Paste$Settings$Builder#epsilon:
:builder_option stabilityThreshold .paste.Paste$Settings$Builder#stabilityThreshold:
:builder_option maxSimultaneousMutations .paste.Paste$Settings$Builder#maxSimultaneousMutations:
:builder_option maxNumPfConfs .paste.Paste$Settings$Builder#maxNumPfConfs:
:builder_option useExternalMemory .paste.Paste$Settings$Builder#useExternalMemory:
:builder_option showPfuncProgress .paste.Paste$Settings$Builder#showPfuncProgress:
:builder_option useWindowCriterion .paste.Paste$Settings$Builder#useWindowCriterion:
:param str addMutFile: Path to the file that has the mutant sequences of interest
:param bool writeSequencesToConsole: True to write sequences and scores to the console
:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

:rtype: :java:ref:`.paste.Paste`

<a name="osprey.KStar"></a>
#### KStar

```python
KStar(proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None, useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault)
```

:java:classdoc:`.kstar.KStar`

For examples using K*, see the examples/python.KStar directory in your Osprey distribution.

**Arguments**:

- `proteinConfSpace`: :java:fielddoc:`.kstar.KStar#protein`
:type proteinConfSpace: :java:ref:`.confspace.SimpleConfSpace`
- `ligandConfSpace`: :java:fielddoc:`.kstar.KStar#ligand`
:type ligandConfSpace: :java:ref:`.confspace.SimpleConfSpace`
- `complexConfSpace`: :java:fielddoc:`.kstar.KStar#complex`
:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
:builder_option epsilon .kstar.KStar$Settings$Builder#epsilon:
:builder_option stabilityThreshold .kstar.KStar$Settings$Builder#stabilityThreshold:
:builder_option maxSimultaneousMutations .kstar.KStar$Settings$Builder#maxSimultaneousMutations:
:builder_option useExternalMemory .kstar.KStar$Settings$Builder#useExternalMemory:
:builder_option showPfuncProgress .kstar.KStar$Settings$Builder#showPfuncProgress:
:param bool writeSequencesToConsole: True to write sequences and scores to the console
:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

:rtype: :java:ref:`.kstar.KStar`

<a name="osprey.PartitionFunction"></a>
#### PartitionFunction

```python
PartitionFunction(confEcalc, confSearchUpper, confSearchLower, rcs)
```

:java:classdoc:`.kstar.pfunc.GradientDescentPfunc`
TODO: docme
:rtype: :java:ref:`.kstar.pfunc.GradientDescentPfunc`

<a name="osprey.MARKStarPfunc"></a>
#### MARKStarPfunc

```python
MARKStarPfunc(confSpace, ematMinimized, confEcalcMinimized, ematRigid, confEcalcRigid, rcs)
```

TODO: docme

<a name="osprey.BBKStar"></a>
#### BBKStar

```python
BBKStar(proteinConfSpace, ligandConfSpace, complexConfSpace, epsilon=useJavaDefault, stabilityThreshold=useJavaDefault, maxSimultaneousMutations=useJavaDefault, energyMatrixCachePattern=useJavaDefault, useExternalMemory=useJavaDefault, showPfuncProgress=useJavaDefault, numBestSequences=useJavaDefault, numConfsPerBatch=useJavaDefault, writeSequencesToConsole=False, writeSequencesToFile=None)
```

:java:classdoc:`.kstar.BBKStar`

For examples using BBK*, see the examples/python.KStar directory in your Osprey distribution.

**Arguments**:

- `proteinConfSpace`: :java:fielddoc:`.kstar.BBKStar#protein`
:type proteinConfSpace: :java:ref:`.confspace.SimpleConfSpace`
- `ligandConfSpace`: :java:fielddoc:`.kstar.BBKStar#ligand`
:type ligandConfSpace: :java:ref:`.confspace.SimpleConfSpace`
- `complexConfSpace`: :java:fielddoc:`.kstar.BBKStar#complex`
:type complexConfSpace: :java:ref:`.confspace.SimpleConfSpace`
:builder_option epsilon .kstar.KStar$Settings$Builder#epsilon:
:builder_option stabilityThreshold .kstar.KStar$Settings$Builder#stabilityThreshold:
:builder_option maxSimultaneousMutations .kstar.KStar$Settings$Builder#maxSimultaneousMutations:
:builder_option useExternalMemory .kstar.KStar$Settings$Builder#useExternalMemory:
:builder_option showPfuncProgress .kstar.KStar$Settings$Builder#showPfuncProgress:
:builder_option numBestSequences .kstar.BBKStar$Settings$Builder#numBestSequences:
:builder_option numConfsPerBatch .kstar.BBKStar$Settings$Builder#numConfsPerBatch:
:param bool writeSequencesToConsole: True to write sequences and scores to the console
:param str writeSequencesToFile: Path to the log file to write sequences scores (in TSV format), or None to skip logging

:rtype: :java:ref:`.kstar.BBKStar`

<a name="osprey.ConfAnalyzer"></a>
#### ConfAnalyzer

```python
ConfAnalyzer(confEcalc)
```

:java:classdoc:`.gmec.ConfAnalyzer`

For examples using the conf analyzer, see examples/python.GMEC/analyzeConf.py in your Osprey distribution.

**Arguments**:

- `confEcalc`: :java:fielddoc:`.gmec.SimpleGMECFinder$Builder#confEcalc`
:type confEcalc: :java:ref:`.energy.ConfEnergyCalculator`

:rtype: :java:ref:`.gmec.ConfAnalyzer`

<a name="osprey.SequenceAnalyzer"></a>
#### SequenceAnalyzer

```python
SequenceAnalyzer(kstar)
```

:java:classdoc:`.kstar.SequenceAnalyzer`

For examples using the sequence analyzer, see examples/python.KStar/analyzeSequence.py in your Osprey distribution.

**Arguments**:

- `kstar`: a configured instance of KStar
:type kstar: :java:ref:`.kstar.KStar`

:rtype: :java:ref:`.kstar.SequenceAnalyzer`

<a name="osprey.LUTE_train"></a>
#### LUTE\_train

```python
LUTE_train(confEcalc, emat, pmat, maxRMSE=0.1, maxOverfittingScore=1.5, randomSeed=12345, confDBPath=None, samplingStrategy=LUTE_SamplingStrategy.Progressive)
```

Trains a LUTE model

For examples using LUTE, see examples/python.GMEC/LUTE.*.py and examples/python.KStar/LUTE.*.py in your Osprey distribution.

**Arguments**:

- `confEcalc`: The conformation energy calculator
:type confEcalc: :java:ref:`.energy.ConfEnergyCalculator`
- `emat`: An energy matrix
:type emat: :java:ref:`.ematrix.EnergyMatrix`
- `pmat`: A pruning matrix, resulting from DEE
:type pmat: :java:ref:`.pruning.PruningMatrix`

:param float maxRMSE: The maximum tolerable fit RMS error
:param float maxOverfittingScore: The maximum tolerable amount of overfitting (score = training set RMSE / test set RMSE)
:param int randomSeed: Random seed to use for conformation sampling
:param str confDBPath: Path to write/read confDB file, or None to omit saving the confDB to disk

**Returns**:

The LUTE model if it meets the accuracy goals, or None otherwise
:rtype: :java:ref:`.lute.LUTEState`

<a name="osprey.LUTE_write"></a>
#### LUTE\_write

```python
LUTE_write(model, path)
```

Writes a LUTE model to a file

**Arguments**:

- `model`: The LUTE model
:type model: :java:ref:`.lute.LUTEState`
:param str path: Path to the file

<a name="osprey.LUTE_read"></a>
#### LUTE\_read

```python
LUTE_read(path)
```

Reads a LUTE model from a file

:param str path: Path to the file

**Returns**:

The LUTE model
:rtype: :java:ref:`.lute.LUTEState`

<a name="osprey.LUTE_ConfEnergyCalculator"></a>
#### LUTE\_ConfEnergyCalculator

```python
LUTE_ConfEnergyCalculator(confSpace, model)
```

Creates a LUTE conformation energy calculator

**Arguments**:

- `confSpace`: The conformation space
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
- `model`: The LUTE model
:type model: :java:ref:`.lute.LUTEState`

:rtype: :java:ref:`.lute.LUTEConfEnergyCalculator`

<a name="osprey.LUTE_AStar"></a>
#### LUTE\_AStar

```python
LUTE_AStar(rcs, pmat, luteEcalc, showProgress=True)
```

:java:methoddoc:`.astar.conf.ConfAStarTree$Builder#setLUTE`

:builder_option rcs .astar.conf.ConfAStarTree$Builder#rcs:

**Arguments**:

- `pmat`: The pruning matrix from the LUTE training calculation.
:type pmat: :java:ref:`.pruning.PruningMatrix`
- `luteEcalc`: The LUTE conformation energy calculator
:type luteEcalc: :java:ref:`.lute.LUTEConfEnergyCalculator`

:builder_return .astar.conf.ConfAStarTree$Builder:

<a name="osprey.LUTE_Pfunc"></a>
#### LUTE\_Pfunc

```python
LUTE_Pfunc(luteEcalc, astar, rcs)
```

:java:classdoc:`.lute.LUTEPfunc`
TODO: docme

<a name="osprey.LUTE_GMECFinder"></a>
#### LUTE\_GMECFinder

```python
LUTE_GMECFinder(confSpace, model, pmat, confLog=useJavaDefault, printIntermediateConfs=useJavaDefault)
```

:java:classdoc:`.lute.LUTEGMECFinder`

**Arguments**:

- `confSpace`: The conformation space
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`
- `model`: The LUTE model
:type model: :java:ref:`.lute.LUTEState`
- `pmat`: The pruning matrix from the LUTE training calculation.
:type pmat: :java:ref:`.pruning.PruningMatrix`
:param str confLog: Path to file where conformations found during conformation space search should be logged.
:builder_option printIntermediateConfs .gmec.SimpleGMECFinder$Builder#printIntermediateConfsToConsole:

:rtype: :java:ref:`.lute.LUTEGMECFinder`

<a name="osprey.COMETS_State"></a>
#### COMETS\_State

```python
COMETS_State(name, confSpace)
```

:java:classdoc:`.gmec.Comets$State`

:param str name: :java:fielddoc:`.gmec.Comets$State#name`

**Arguments**:

- `confSpace`: :java:fielddoc:`.gmec.Comets$State#confSpace`
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

:rtype: :java:ref:`.gmec.Comets$State`

<a name="osprey.COMETS_LME"></a>
#### COMETS\_LME

```python
COMETS_LME(weightsByState, offset=useJavaDefault, constrainLessThan=None)
```

:java:classdoc:`.gmec.Comets$LME`

**Arguments**:

- `weightsByState`: map from states to weights
:type weightsByState: map from :java:ref:`.gmec.Comets$State` to float

:builder_option offset .gmec.Comets$LME$Builder#offset:
:param float constrainLessThan: :java:methoddoc:`.gmec.Comets$LME$Builder#constrainLessThan`
:builder_return .gmec.Comets$LME$Builder:

<a name="osprey.COMETS"></a>
#### COMETS

```python
COMETS(objective, constraints=[], objectiveWindowSize=useJavaDefault, objectiveWindowMax=useJavaDefault, maxSimultaneousMutations=useJavaDefault, minNumConfTrees=useJavaDefault, logFile=None)
```

:java:classdoc:`.gmec.Comets`

:builder_option objective .gmec.Comets$Builder#objective:

**Arguments**:

- `constraints`: List of LMEs to use as constraints
:type constraints: list of :java:ref:`.gmec.Comets$LME`

:builder_option objectiveWindowSize .gmec.Comets$Builder#objectiveWindowSize:
:builder_option objectiveWindowMax .gmec.Comets$Builder#objectiveWindowMax:
:builder_option maxSimultaneousMutations .gmec.Comets$Builder#maxSimultaneousMutations:
:builder_option minNumConfTrees .gmec.Comets$Builder#minNumConfTrees:

:param str logFile: :java:fielddoc:`.gmec.Comets$Builder#logFile`

:builder_return .gmec.Comets$Builder:

<a name="osprey.MSKStar_State"></a>
#### MSKStar\_State

```python
MSKStar_State(name, confSpace)
```

:java:classdoc:`.kstar.MSKStar$State`

:param str name: :java:fielddoc:`.kstar.MSKStar$State#name`

**Arguments**:

- `confSpace`: :java:fielddoc:`.kstar.MSKStar$State#confSpace`
:type confSpace: :java:ref:`.confspace.SimpleConfSpace`

:rtype: :java:ref:`.kstar.MSKStar$State`

<a name="osprey.MSKStar_LMFE"></a>
#### MSKStar\_LMFE

```python
MSKStar_LMFE(weightsByState, offset=useJavaDefault, constrainLessThan=None)
```

:java:classdoc:`.kstar.MSKStar$LMFE`

**Arguments**:

- `weightsByState`: map from states to weights
:type weightsByState: map from :java:ref:`.kstar.MSKStar$State` to float

:builder_option offset .kstar.MSKStar$LMFE$Builder#offset:
:param float constrainLessThan: :java:methoddoc:`.kstar.MSKStar$LMFE$Builder#constrainLessThan`
:builder_return .kstar.MSKStar$LMFE$Builder:

<a name="osprey.MSKStar"></a>
#### MSKStar

```python
MSKStar(objective, constraints=[], epsilon=useJavaDefault, objectiveWindowSize=useJavaDefault, objectiveWindowMax=useJavaDefault, maxSimultaneousMutations=useJavaDefault, minNumConfTrees=useJavaDefault, logFile=None)
```

:java:classdoc:`.kstar.MSKStar`

:builder_option objective .kstar.MSKStar$Builder#objective:

**Arguments**:

- `constraints`: List of LMFEs to use as constraints
:type constraints: list of :java:ref:`.kstar.MSKStar$LMFE`

:builder_option epsilon .kstar.MSKStar$Builder#epsilon:
:builder_option objectiveWindowSize .kstar.MSKStar$Builder#objectiveWindowSize:
:builder_option objectiveWindowMax .kstar.MSKStar$Builder#objectiveWindowMax:
:builder_option maxSimultaneousMutations .kstar.MSKStar$Builder#maxSimultaneousMutations:
:builder_option minNumConfTrees .kstar.MSKStar$Builder#minNumConfTrees:

:param str logFile: :java:fielddoc:`.kstar.MSKStar$Builder#logFile`

:builder_return .kstar.MSKStar$Builder:

<a name="osprey.SOFEA_StateConfig"></a>
#### SOFEA\_StateConfig

```python
SOFEA_StateConfig(emat, confEcalc, confdbPath=None)
```

:java:classdoc:`.sofea.Sofea$StateConfig`

**Arguments**:

- `emat`: :java:fielddoc:`.sofea.Sofea$StateConfig#emat`
- `confEcalc`: :java:fielddoc:`.sofea.Sofea$StateConfig#confEcalc`
- `confdbPath`: :java:fielddoc:`.sofea.Sofea$StateConfig#confDBFile`
:rtype: :java:ref:`.sofea.Sofea$StateConfig`

<a name="osprey.SOFEA"></a>
#### SOFEA

```python
SOFEA(confSpace, configFunc, seqdbPath='sofea.seqdb', seqdbMathContext=useJavaDefault, fringedbLowerPath='sofea.lower.fringedb', fringedbLowerMiB=10, fringedbUpperPath='sofea.upper.fringedb', fringedbUpperMiB=10, rcdbPath=useJavaDefault, showProgress=useJavaDefault, performanceLogPath=useJavaDefault, sweepIncrement=useJavaDefault, maxNumMinimizations=useJavaDefault, negligableFreeEnergy=useJavaDefault)
```

:java:classdoc:`.sofea.Sofea`

**Arguments**:

- `confSpace`: A multi-state configuration space
:type confSpace: :java:ref:`.confspace.MultiStateConfSpace`

- `configFunc`: a function that creates a :java:ref:`.sofea.Sofea$StateConfig` for a state
:type configFunc: function(:java:ref:`.confspace.MultiStateConfSpace$State`) returning :java:ref:`.sofea.Sofea$StateConfig`

:param str seqdbPath: Path to write the sequence database file
:builder_option seqdbMathContext .sofea.Sofea$Builder#seqdbMathContext:
:param str fringedbLowerPath: Path to write the lower fringe set
:param int fringedbLowerMiB: size of the lower fringe set in MiB
:param str fringedbUpperPath: Path to write the upper fringe set
:param int fringedbUpperMiB: size of the upper fringe set in MiB
:param str rcdbPath: Path to write the upper fringe set
:builder_option showProgress .sofea.Sofea$Builder#showProgress:
:param str performanceLogPath: Path to write the performance log
:builder_option sweepIncrement .sofea.Sofea$Builder#sweepIncrement:
:builder_option maxNumMinimizations .sofea.Sofea$Builder#maxNumMinimizations:
:builder_option negligableFreeEnergy .sofea.Sofea$Builder#negligableFreeEnergy:

:builder_return .sofea.Sofea$Builder:

<a name="osprey.SOFEA_MinLMFE"></a>
#### SOFEA\_MinLMFE

```python
SOFEA_MinLMFE(lmfe, numSequences, minFreeEnergyWidth)
```

:java:classdoc:`.sofea.MinLMFE`

**Arguments**:

- `lmfe`: :java:fielddoc:`.sofea.MinLMFE#objective`
- `numSequences`: :java:fielddoc:`.sofea.MinLMFE#numSequences`
- `minFreeEnergyWidth`: :java:fielddoc:`.sofea.MinLMFE#minFreeEnergyWidth`
:rtype: :java:ref:`.sofea.MinLMFE`

<a name="osprey.SOFEA_SequenceLMFE"></a>
#### SOFEA\_SequenceLMFE

```python
SOFEA_SequenceLMFE(sequence, lmfe, minFreeEnergyWidth)
```

:java:classdoc:`.sofea.SequenceLMFE`

**Arguments**:

- `sequence`: :java:fielddoc:`.sofea.SequenceLMFE#seq`
- `lmfe`: :java:fielddoc:`.sofea.SequenceLMFE#lmfe`
- `minFreeEnergyWidth`: :java:fielddoc:`.sofea.SequenceLMFE#minFreeEnergyWidth`
:rtype: :java:ref:`.sofea.SequenceLMFE`

