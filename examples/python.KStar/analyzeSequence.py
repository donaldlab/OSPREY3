
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# make sure all strands share the same template library (including wild-type rotamers)
templateLib = osprey.TemplateLibrary(ffparams.forcefld, moleculesForWildTypeRotamers=[mol])

# define the protein strand
protein = osprey.Strand(mol, templateLib=templateLib, residues=['G648', 'G654'])
protein.flexibility['G649'].setLibraryRotamers(osprey.WILD_TYPE, 'TYR', 'ALA', 'VAL', 'ILE', 'LEU').addWildTypeRotamers().setContinuous()
protein.flexibility['G650'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['G651'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
protein.flexibility['G654'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# define the ligand strand
ligand = osprey.Strand(mol, templateLib=templateLib, residues=['A155', 'A194'])
ligand.flexibility['A156'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A172'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A192'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()
ligand.flexibility['A193'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers().setContinuous()

# make the conf space for the protein
proteinConfSpace = osprey.ConfSpace(protein)

# make the conf space for the ligand
ligandConfSpace = osprey.ConfSpace(ligand)

# make the conf space for the protein+ligand complex
complexConfSpace = osprey.ConfSpace([protein, ligand])

# how should we compute energies of molecules?
# (give the complex conf space to the ecalc since it knows about all the templates and degrees of freedom)
parallelism = osprey.Parallelism(cpuCores=4)
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

# how should we define energies of conformations?
def confEcalcFactory(confSpace, ecalc):
	eref = osprey.ReferenceEnergies(confSpace, ecalc)
	return osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# how should confs be ordered and searched?
def astarFactory(emat, rcs):
	return osprey.AStarTraditional(emat, rcs, showProgress=False)

# get the sequence analyzer
analyzer = osprey.SequenceAnalyzer(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	ecalc,
	confEcalcFactory,
	astarFactory,
	energyMatrixCachePattern='emat.*.dat'
)

# how big should the ensembles be? (in terms of energy from the min)
energyWindowSize = 1.0

# analyze the wild-type complex
analysis = analyzer.analyze(
	complexConfSpace.makeWildTypeSequence(),
	energyWindowSize
)
print('\n')
print(analysis)
analysis.writePdbs('ensemble-wt-complex/conf.*.pdb')

# analyze the wild-type unbound protein
analysis = analyzer.analyze(
	proteinConfSpace.makeWildTypeSequence(),
	energyWindowSize
)
print('\n')
print(analysis)
analysis.writePdbs('ensemble-wt-protein/conf.*.pdb')

# analyze a mutant complex
analysis = analyzer.analyze(
	complexConfSpace.makeWildTypeSequence().set("G649", "ILE"),
	energyWindowSize
)
print('\n')
print(analysis)
analysis.writePdbs('ensemble-ile/conf.*.pdb')


# if you've run the kstar.confdb.py or bbkstar.confdb.py examples, then we can analyze
# sequences using the pre-computed conf DB
# otherwise, nothing will happen here
analyzer = osprey.SequenceAnalyzer(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	ecalc,
	confEcalcFactory,
	astarFactory,
	energyMatrixCachePattern='emat.*.dat',
	confDBPattern='conf.*.db'
)

print('\nreading from previous conf DB...')
analysis = analyzer.analyzeFromConfDB(
	complexConfSpace.makeWildTypeSequence(),
	energyWindowSize
)
print('\n')
print(analysis)
