import osprey

osprey.start()

# define a strand
strand = osprey.Strand('1CC8.ss.pdb')
strand.flexibility['A38'].setLibraryRotamers('ILE','ALA', 'GLY').setContinuous()
strand.flexibility['A39'].setLibraryRotamers("Ser","Ala","Gly").setContinuous()
strand.flexibility['A40'].setLibraryRotamers("Met","Ser","Ala","Gly","Leu").setContinuous()
strand.flexibility['A41'].setLibraryRotamers("Glu","Ala","Gly").setContinuous()

# make the conf space
bbflex = osprey.DEEPerStrandFlex(strand, '1CC8.d.pert', ['A38','A39','A40','A41'], '1CC8.ss.pdb')
confSpace = osprey.ConfSpace([[strand,bbflex]])

# choose a forcefield
ffparams = osprey.ForcefieldParams()
ffparams.solvScale = 0.

# how should we compute energies of molecules?
ecalc = osprey.EnergyCalculator(confSpace, ffparams)

eref = osprey.ReferenceEnergies(confSpace, ecalc)

# how should we define energies of conformations?
confEcalc = osprey.ConfEnergyCalculator(confSpace, ecalc, referenceEnergies=eref)

# how should confs be ordered and searched?
emat = osprey.EnergyMatrix(confEcalc)

# find the best sequence and rotamers
gmec = osprey.DEEGMECFinder(emat, confSpace, ecalc, confEcalc, '1CC8.dee.pert', True, True).calcGMEC()

