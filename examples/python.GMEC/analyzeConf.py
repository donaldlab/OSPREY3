
import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility['A2'].setLibraryRotamers('ALA', 'GLY')
strand.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE, 'VAL')
strand.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE)

# make the conf space
confSpace = osprey.ConfSpace(strand)

# choose a forcefield
ffparams = osprey.ForcefieldParams()

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

# how should we define energies of conformations?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc)

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confEcalc)
astar = osprey.AStarMPLP(emat, confSpace)

# find the best sequence and rotamers
# and save all the conformations in the confdb for easier analysis later
confDBFile = 'conf.db'
gmec = osprey.GMECFinder(astar, confEcalc, confDBFile=confDBFile).find()

# get the conformation analyzer
analyzer = osprey.ConfAnalyzer(confEcalc, emat)

# analyze the gmec
analysis = analyzer.analyze(gmec)
print('\n')
print(analysis)
osprey.writePdb(analysis.epmol, 'gmec.pdb')

# breakdown the forcefield energy by position
print('Forcefield breakdown:')
print(analysis.breakdownEnergyByPosition(osprey.BreakdownType.All))

# breakdown the forcefield energy by position and effect
print('Electrostatic breakdown:')
print(analysis.breakdownEnergyByPosition(osprey.BreakdownType.Electrostatics))
print('van der Waals breakdown:')
print(analysis.breakdownEnergyByPosition(osprey.BreakdownType.VanDerWaals))
print('Solvation breakdown:')
print(analysis.breakdownEnergyByPosition(osprey.BreakdownType.Solvation))
# show per-residue and per-residue-pair energies
# (includes things like reference energies, residue entropies)
print('Offsets breakdown:')
print(analysis.breakdownEnergyByPosition(osprey.BreakdownType.Offsets))

# breakdown the A* score by position
print('Score breakdown:')
print(analysis.breakdownScoreByPosition())

# analyze an arbitrary conformation from this design
analysis = analyzer.analyze([0, 0, 0])
# get the numbers from the "Residue Conf IDs" from the console output
# or the "CONF" section of the conformations log
print('\n')
print(analysis)

# analyze the low-energy ensemble from the energy window
# and put any extra conformations we find in the confDB
astar = osprey.AStarMPLP(emat, confSpace)
confs = osprey.GMECFinder(astar, confEcalc, confDBFile=confDBFile).find(1)
maxNumConfs = 10
ensembleAnalysis = analyzer.analyzeEnsemble(confs, maxNumConfs)
print('\n')
print(ensembleAnalysis)
ensembleAnalysis.writePdbs('ensemble/conf.*.pdb')

# we can also just analyze the contents of the confDB without doing additional energy calculations
ensembleAnalysis = analyzer.analyzeGMECEnsembleFromConfDB(confDBFile, maxNumConfs)
print('\n')
print(ensembleAnalysis)
