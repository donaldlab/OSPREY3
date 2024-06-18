import osprey

osprey.start()

# default template library contains templates, coords, entropies etc for natual amino acids
defaultTemplateLib = osprey.TemplateLibrary()

# define our custom templates
customTemplates = """
Ignored comment 1
Ignored comment 2
ALANINE
 ALA  INT     1                                                 
 CORR OMIT DU   BEG                                             
   0.00000                                                      
   1  DUMM  DU    M    0  -1  -2     0.000     0.000     0.000  0.00000 
   2  DUMM  DU    M    1   0  -1     1.449     0.000     0.000  0.00000 
   3  DUMM  DU    M    2   1   0     1.522   111.100     0.000  0.00000 
   4  N     N     M    3   2   1     1.335   116.600   180.000 -0.41570 
   5  H     H     E    4   3   2     1.010   119.800     0.000  0.27190 
   6  CA    CT    M    4   3   2     1.449   121.900   180.000  0.03370 
   7  HA    H1    E    6   4   3     1.090   109.500   300.000  0.08230 
   8  CB    CT    3    6   4   3     1.525   111.100    60.000 -0.18250 
   9  HB1   HC    E    8   6   4     1.090   109.500    60.000  0.06030 
  10  HB2   HC    E    8   6   4     1.090   109.500   180.000  0.06030 
  11  HB3   HC    E    8   6   4     1.090   109.500   300.000  0.06030 
  12  C     C     M    6   4   3     1.522   111.100   180.000  0.59730 
  13  O     O     E   12   6   4     1.229   120.500     0.000 -0.56790 

IMPROPER                                                        
 -M   CA   N    H                                               
 CA   +M   C    O                                               
                                                                
DONE                                                            
"""

# define our custom template coordinates
customTemplateCoords = """
ALA 10 
N -0.677f -1.23f -0.491f -0.4157f N
H -0.131f -2.162f -0.491f 0.2719f H
CA -0.001f 0.064f -0.491f 0.03370f CT
HA -0.28f 0.596f -1.413f 0.0823f H1
CB -0.432f 0.894f 0.709f -0.1825f CT
HB1 0.09f 1.863f 0.692f 0.0603f HC
HB2 -1.518f 1.063f 0.667f 0.0603f HC
HB3 -0.179f 0.359f 1.636f 0.0603f HC
C 1.499f -0.11f -0.491f 0.5973f C
O 2.065f -0.922f 0.251f -0.5679f O
ENDRES 
"""

# define our custom rotamers
customRotamers = """
! The first line is the number of AA types to be read
! The format for the rest of the file is
! AA_name num_dihedrals num_rotamers
! dihedral_list_one_per_line
! rotamer_angles
1
VAL 1 3
N CA CB CG1
64
175
-60
"""

# build the customized template library
customizedTemplateLib = osprey.TemplateLibrary(
	extraTemplates=[customTemplates],
	extraTemplateCoords=[customTemplateCoords],
	extraRotamers=[customRotamers]
)

# or read templates from files
# customizedTemplateLibFromFiles = osprey.TemplateLibrary(
# 	extraTemplates=['/path/to/templates/file']
# 	etc...
# )

# or completely replace default templates
# completelyCustomTemplateLib = osprey.TemplateLibrary(
# 	defaultTemplates=False,
# 	extraTemplates=['/path/to/all/templates'],
# 	etc...
# )

# load the molecule and make the strand using our custom template library
protein = osprey.Strand('1CC8.ss.pdb', templateLib=customizedTemplateLib)

# make the conf space
protein.flexibility['A2'].setLibraryRotamers('ALA', 'GLY')
protein.flexibility['A3'].setLibraryRotamers(osprey.WILD_TYPE, 'VAL', 'ARG').setContinuous(10)
protein.flexibility['A4'].setLibraryRotamers(osprey.WILD_TYPE).addWildTypeRotamers()
confSpace = osprey.ConfSpace(protein)

# continue design with confSpace
