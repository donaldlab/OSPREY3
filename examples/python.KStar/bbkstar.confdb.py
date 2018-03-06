
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
minimizingEcalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism, isMinimizing=True)

# BBK* needs a rigid energy calculator too, for multi-sequence bounds on K*
rigidEcalc = osprey.SharedEnergyCalculator(minimizingEcalc, isMinimizing=False)

# how should we define energies of conformations?
def confEcalcFactory(confSpace, ecalc):
	eref = osprey.ReferenceEnergies(confSpace, ecalc)
	return osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# how should confs be ordered and searched?
def astarFactory(emat, rcs):
	return osprey.AStarTraditional(emat, rcs, showProgress=False)
	# or
	# return osprey.AStarMPLP(emat, rcs, numIterations=5)

# run BBK* using a confDB
bbkstar = osprey.BBKStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	rigidEcalc,
	minimizingEcalc,
	confEcalcFactory,
	astarFactory,
	numBestSequences=2,
	epsilon=0.5, # let's use a smaller epsilon so the design takes a noticeable amount of time
	energyMatrixCachePattern='emat.*.dat',
	confDBPattern='conf.*.db', # actually several confDBs will be written, so give a name pattern
	writeSequencesToConsole=True,
	writeSequencesToFile='bbkstar.results.tsv'
)
scoredSequences = bbkstar.run()

# If OSPREY exits unexpectedly, the 'conf.*.db' files will still remain.
# running BBK* again using an existing confDB should be MUCH faster!
# meaning, your design should catch up to where it left off very quickly.
print('\nRunning BBK* again...\n')
osprey.BBKStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	rigidEcalc,
	minimizingEcalc,
	confEcalcFactory,
	astarFactory,
	numBestSequences=2,
	epsilon=0.5, # let's use a smaller epsilon so the design takes a noticeable amount of time
	energyMatrixCachePattern='emat.*.dat',
	confDBPattern='conf.*.db', # actually several confDBs will be written, so give a name pattern
	writeSequencesToConsole=True,
	writeSequencesToFile='bbkstar.results.tsv'
).run()
print('\nSecond BBK* run finished!\n')