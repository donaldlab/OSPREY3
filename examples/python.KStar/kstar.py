
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
ecalc = osprey.EnergyCalculator(complexConfSpace, ffparams, parallelism=parallelism)

# configure K*
kstar = osprey.KStar(
	proteinConfSpace,
	ligandConfSpace,
	complexConfSpace,
	epsilon=0.99, # you proabably want something more precise in your real designs
	writeSequencesToConsole=True,
	writeSequencesToFile='kstar.results.tsv'
)

# configure K* inputs for each conf space
for info in kstar.confSpaceInfos():

	# how should we define energies of conformations?
	eref = osprey.ReferenceEnergies(info.confSpace, ecalc)
	info.confEcalc = osprey.ConfEnergyCalculator(info.confSpace, ecalc, referenceEnergies=eref)

	# compute the energy matrix
	emat = osprey.EnergyMatrix(info.confEcalc, cacheFile='emat.%s.dat' % info.id)

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
scoredSequences = kstar.run(ecalc.tasks)

# make a sequence analyzer to look at the results
analyzer = osprey.SequenceAnalyzer(kstar)

# use results
for scoredSequence in scoredSequences:
	print("result:")
	print("\tsequence: %s" % scoredSequence.sequence)
	print("\tK* score: %s" % scoredSequence.score)

	# write the sequence ensemble, with up to 10 of the lowest-energy conformations
	numConfs = 10
	analysis = analyzer.analyze(scoredSequence.sequence, numConfs)
	print(analysis)
	analysis.writePdb(
		'seq.%s.pdb' % scoredSequence.sequence,
		'Top %d conformations for sequence %s' % (numConfs, scoredSequence.sequence)
	)
