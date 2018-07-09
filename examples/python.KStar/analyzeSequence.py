
import osprey
osprey.start()

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# read a PDB file for molecular info
mol = osprey.readPdb('2RL0.min.reduce.pdb')

# make sure all strands share the same template library
templateLib = osprey.TemplateLibrary(ffparams.forcefld)

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

# configure K*
kstar = osprey.KStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace
)

# configure K* inputs for each conf space
for info in kstar.confSpaceInfos():

	# how should we define energies of conformations?
	eref = osprey.ReferenceEnergies(info.confSpace, ecalc)
	info.confEcalc = osprey.ConfEnergyCalculator(info.confSpace, ecalc, referenceEnergies=eref)

	# compute the energy matrix
	emat = osprey.EnergyMatrix(info.confEcalc, cacheFile='emat.%s.dat' % info.id)

	# how should confs be ordered and searched? (don't forget to capture emat by using a defaulted argument)
	def makeAStar(rcs, emat=emat):
		return osprey.AStarTraditional(emat, rcs, showProgress=False)
	info.confSearchFactory = osprey.KStar.ConfSearchFactory(makeAStar)

	# if you've run the kstar.confdb.py example, then we can analyze
	# sequences using the pre-computed conf DB files
	# this will make the analyses go much faster
	info.setConfDBFile('kstar.%s.db' % info.type.name().lower())

# make a sequence analyzer from the configured KStar instance
# (you could also give it a configured BBKStar instance if you have that instead)
analyzer = osprey.SequenceAnalyzer(kstar)

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
