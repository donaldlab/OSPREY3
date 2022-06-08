+++
title = "RESISTOR Tutorial"
weight = 1
+++

In this tutorial, I demonstrate how to use RESISTOR to predict Epidermal Growth
Factor Receptor (EGFR) resistance mutations to erlotinib in lung cancer. First,
we must gather the inputs we need for the algorithm. These consist of:

1. Structures of EGFR bound to its endogenous ligand (in this case, an ATP analog) and EGFR bound to erlotinib
1. Connectivity templates and forcefield parameters modifications for the non-amino acid molecules (ATP & erlotinib)
1. Extra rotamers for the protein or ligand if desired or necessary
1. A FASTA file with the coding DNA sequence corresponding to the protein where resistance can occur (in this case, EGFR)
1. A JSON file with the cancer type-specific codon-to-codon mutational probabilities

#### Obtaining and protonating the structures.

Obtain the structures of the protein bound to its endogenous ligand and the
protein bound to the drug. In this tutorial, we use PDB id 2itx for a structure
of EGFR bound to an ATP analog, ANP. For the protein:drug bound complex we use
PDB id 1m17. Then we remove any existing protons and water molecules we are not
interested in modeling. Lastly, we protonate the structure.  The instructions
below assume you have installed AmberTools from ambermd.org and the tools are
available on your $PATH.

The following shell script can be used to prepare the structures for input to RESISTOR:

```shell
FILE=2itx.pdb

# remove the hydrogens                        \
reduce -Trim $FILE                          | \
# remove the waters                           \
grep --invert-match HOH                     | \
# add hydrogens                               \
reduce -                                    | \
# remove the 'new' tags that reduce added     \
sed 's/new//'                               | \
# keep only ATOM and HETATM records           \
grep --extended-regexp '^(ATOM|HETATM)'     | \
# remove trailing whitespace                  \
sed 's/[[:blank:]]*$//' > "${FILE%.*}".reduced.pdb
```

You need to split each structure into its protein and ligand parts. 
These commands are one example of doing this:

```shell
grep --invert-match ANP 2itx.reduced.pdb > 2itx.egfr.pdb
grep ANP 2itx.reduced.pdb > 2itx.anp.pdb

grep --invert-match AQ4 1m17.reduced.pdb > 1m17.egfr.pdb
grep AQ4 1m17.reduced.pdb > 1m17.aq4.pdb
```

Two of the required inputs are the connectivity templates and forcefield parameters for the small molecules. 
`antechamber` generates the templates, and `parmchk2` the forcefield parameters:

```shell
antechamber -i 2itx.anp.pdb -fi pdb -o anp.prepi -fo prepi
antechamber -i 1m17.aq4.pdb -fi pdb -o aq4.prepi -fo prepi

parmchk2 -i anp.prepi -f prepi -a Y -o anp.frcmod
parmchk2 -i aq4.prepi -f prepi -a Y -o aq4.frcmod
```

We also need template coordinates for these templates.
Here's a shell script you can use to generate them in this system:

```shell
echo "AQ4" $(wc --lines < 1m17.aq4.pdb) > aq4.extra-template-coordinates
tr --squeeze-repeats ' ' < 1m17.aq4.pdb | cut --delimiter=' ' --fields=3,7,8,9 >> aq4.extra-template-coordinates

echo "ANP" $(wc --lines < 2itx.anp.pdb) > anp.extra-template-coordinates
tr --squeeze-repeats ' ' < 2itx.anp.pdb | cut --delimiter=' ' --fields=3,7,8,9 >> anp.extra-template-coordinates
```

While OSPREY contains a rotamer library for all of the canonical amino acids, we need to create rotamers for ANP and erlotinib.
Normally you'll want to add multiple rotamers, but as this is a minimal example we'll add just one.
Recall that OSPREY's continuous minimization of rotamers allows dihedrals to flex within a certain voxel in order to minimize the energy.

Here's the rotamer for ANP:

```text
! The first line is the number of AA types to be read
! The format for the rest of the file is
! AA_name num_dihedrals num_rotamers
! dihedral_list_one_per_line
! rotamer_angles
1
ANP 11 1
H61  N6  C6  C5
C8   C9  C1’ C2’
C1’  C2’ O2’ HO2’
HO3’ O3’ C3’ C4’
C3’  C4’ C5’ O5’
C4’  C5’ O5’ PA
C5’  O5’ PA  O3A
O5’  PA  O3A PB
PA   O3A PB  N3B
O3A  PB  N3B PG
PB   N3B PG  O2G
-176 -87 -55 68 81 -102 -43 -96 -93 -40 -157
```

... and here's the rotamer for erlotinib:

```text
! The first line is the number of AA types to be read
! The format for the rest of the file is
! AA_name num_dihedrals num_rotamers
! dihedral_list_one_per_line
! rotamer_angles
1
AQ4 14 1
H161 C16 O4  C15
C16  O4  C15 C14
O4   C15 C14 O3
C15  C14 O3  C13
C14  O3  C13 C9
C13  C9  O1  C10
C9   O1  C10 C11
O1   C10 C11 O2
C10  C11 O2  C12
C11  O2  C12 H122
C7   C6  N1  C5
C6   N1  C5  C4
C4   C3  C2  C1
C3   C2  C1  H1
-178 -178 73 -174 -12 -134 -160 -78 96 -55 -169 37 -5 -32
```

Now you have all of the inputs you need to create the positive and negative K\* affinity designs.
For this, we use YAML-specified files.
YAML is a whitespace-sensitive markup language; for a good introduction into YAML see [this website](https://learnxinyminutes.com/docs/yaml/).
Start with an empty YAML design template:

```yaml
osprey_version: '3.1'
design_name: EGFR + erlotinib
scan:
  residues: []
protein:
  coordinates: |+2
  residue_configurations: []
  extra_templates: ''
  extra_templates_coordinates: ''
  extra_rotamers: ''
ligand:
  coordinates: |+2
  residue_configurations: []
  extra_templates: |+2
  extra_templates_coordinates: |+2
  extra_rotamers: |+2
  translate_rotate: true
```

and then fill out the key-value pairs with the structures, templates, and rotamers you generated above. 
We include our completed design templates in `completed-ANP-template.yaml` and `completed-erlotinib-template.yaml`.

For the calculation of the per-codon mutational signature probabilities, obtain a coding DNA sequence of the protein. 
We'll use ENST00000275493.6 obtained from the [COSMIC website](https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=EGFR).

EGFR-cDNA.fasta
```text
>EGFR
ATGCGACCCTCCGGGACGGCCGGGGCAGCGCTCCTGGCGCTGCTGGCTGCGCTCTGCCCG
GCGAGTCGGGCTCTGGAGGAAAAGAAAGTTTGCCAAGGCACGAGTAACAAGCTCACGCAG
TTGGGCACTTTTGAAGATCATTTTCTCAGCCTCCAGAGGATGTTCAATAACTGTGAGGTG
GTCCTTGGGAATTTGGAAATTACCTATGTGCAGAGGAATTATGATCTTTCCTTCTTAAAG
ACCATCCAGGAGGTGGCTGGTTATGTCCTCATTGCCCTCAACACAGTGGAGCGAATTCCT
TTGGAAAACCTGCAGATCATCAGAGGAAATATGTACTACGAAAATTCCTATGCCTTAGCA
GTCTTATCTAACTATGATGCAAATAAAACCGGACTGAAGGAGCTGCCCATGAGAAATTTA
CAGGAAATCCTGCATGGCGCCGTGCGGTTCAGCAACAACCCTGCCCTGTGCAACGTGGAG
AGCATCCAGTGGCGGGACATAGTCAGCAGTGACTTTCTCAGCAACATGTCGATGGACTTC
CAGAACCACCTGGGCAGCTGCCAAAAGTGTGATCCAAGCTGTCCCAATGGGAGCTGCTGG
GGTGCAGGAGAGGAGAACTGCCAGAAACTGACCAAAATCATCTGTGCCCAGCAGTGCTCC
GGGCGCTGCCGTGGCAAGTCCCCCAGTGACTGCTGCCACAACCAGTGTGCTGCAGGCTGC
ACAGGCCCCCGGGAGAGCGACTGCCTGGTCTGCCGCAAATTCCGAGACGAAGCCACGTGC
AAGGACACCTGCCCCCCACTCATGCTCTACAACCCCACCACGTACCAGATGGATGTGAAC
CCCGAGGGCAAATACAGCTTTGGTGCCACCTGCGTGAAGAAGTGTCCCCGTAATTATGTG
GTGACAGATCACGGCTCGTGCGTCCGAGCCTGTGGGGCCGACAGCTATGAGATGGAGGAA
GACGGCGTCCGCAAGTGTAAGAAGTGCGAAGGGCCTTGCCGCAAAGTGTGTAACGGAATA
GGTATTGGTGAATTTAAAGACTCACTCTCCATAAATGCTACGAATATTAAACACTTCAAA
AACTGCACCTCCATCAGTGGCGATCTCCACATCCTGCCGGTGGCATTTAGGGGTGACTCC
TTCACACATACTCCTCCTCTGGATCCACAGGAACTGGATATTCTGAAAACCGTAAAGGAA
ATCACAGGGTTTTTGCTGATTCAGGCTTGGCCTGAAAACAGGACGGACCTCCATGCCTTT
GAGAACCTAGAAATCATACGCGGCAGGACCAAGCAACATGGTCAGTTTTCTCTTGCAGTC
GTCAGCCTGAACATAACATCCTTGGGATTACGCTCCCTCAAGGAGATAAGTGATGGAGAT
GTGATAATTTCAGGAAACAAAAATTTGTGCTATGCAAATACAATAAACTGGAAAAAACTG
TTTGGGACCTCCGGTCAGAAAACCAAAATTATAAGCAACAGAGGTGAAAACAGCTGCAAG
GCCACAGGCCAGGTCTGCCATGCCTTGTGCTCCCCCGAGGGCTGCTGGGGCCCGGAGCCC
AGGGACTGCGTCTCTTGCCGGAATGTCAGCCGAGGCAGGGAATGCGTGGACAAGTGCAAC
CTTCTGGAGGGTGAGCCAAGGGAGTTTGTGGAGAACTCTGAGTGCATACAGTGCCACCCA
GAGTGCCTGCCTCAGGCCATGAACATCACCTGCACAGGACGGGGACCAGACAACTGTATC
CAGTGTGCCCACTACATTGACGGCCCCCACTGCGTCAAGACCTGCCCGGCAGGAGTCATG
GGAGAAAACAACACCCTGGTCTGGAAGTACGCAGACGCCGGCCATGTGTGCCACCTGTGC
CATCCAAACTGCACCTACGGATGCACTGGGCCAGGTCTTGAAGGCTGTCCAACGAATGGG
CCTAAGATCCCGTCCATCGCCACTGGGATGGTGGGGGCCCTCCTCTTGCTGCTGGTGGTG
GCCCTGGGGATCGGCCTCTTCATGCGAAGGCGCCACATCGTTCGGAAGCGCACGCTGCGG
AGGCTGCTGCAGGAGAGGGAGCTTGTGGAGCCTCTTACACCCAGTGGAGAAGCTCCCAAC
CAAGCTCTCTTGAGGATCTTGAAGGAAACTGAATTCAAAAAGATCAAAGTGCTGGGCTCC
GGTGCGTTCGGCACGGTGTATAAGGGACTCTGGATCCCAGAAGGTGAGAAAGTTAAAATT
CCCGTCGCTATCAAGGAATTAAGAGAAGCAACATCTCCGAAAGCCAACAAGGAAATCCTC
GATGAAGCCTACGTGATGGCCAGCGTGGACAACCCCCACGTGTGCCGCCTGCTGGGCATC
TGCCTCACCTCCACCGTGCAGCTCATCACGCAGCTCATGCCCTTCGGCTGCCTCCTGGAC
TATGTCCGGGAACACAAAGACAATATTGGCTCCCAGTACCTGCTCAACTGGTGTGTGCAG
ATCGCAAAGGGCATGAACTACTTGGAGGACCGTCGCTTGGTGCACCGCGACCTGGCAGCC
AGGAACGTACTGGTGAAAACACCGCAGCATGTCAAGATCACAGATTTTGGGCTGGCCAAA
CTGCTGGGTGCGGAAGAGAAAGAATACCATGCAGAAGGAGGCAAAGTGCCTATCAAGTGG
ATGGCATTGGAATCAATTTTACACAGAATCTATACCCACCAGAGTGATGTCTGGAGCTAC
GGGGTGACTGTTTGGGAGTTGATGACCTTTGGATCCAAGCCATATGACGGAATCCCTGCC
AGCGAGATCTCCTCCATCCTGGAGAAAGGAGAACGCCTCCCTCAGCCACCCATATGTACC
ATCGATGTCTACATGATCATGGTCAAGTGCTGGATGATAGACGCAGATAGTCGCCCAAAG
TTCCGTGAGTTGATCATCGAATTCTCCAAAATGGCCCGAGACCCCCAGCGCTACCTTGTC
ATTCAGGGGGATGAAAGAATGCATTTGCCAAGTCCTACAGACTCCAACTTCTACCGTGCC
CTGATGGATGAAGAAGACATGGACGACGTGGTGGATGCCGACGAGTACCTCATCCCACAG
CAGGGCTTCTTCAGCAGCCCCTCCACGTCACGGACTCCCCTCCTGAGCTCTCTGAGTGCA
ACCAGCAACAATTCCACCGTGGCTTGCATTGATAGAAATGGGCTGCAAAGCTGTCCCATC
AAGGAAGACAGCTTCTTGCAGCGATACAGCTCAGACCCCACAGGCGCCTTGACTGAGGAC
AGCATAGACGACACCTTCCTCCCAGTGCCTGAATACATAAACCAGTCCGTTCCCAAAAGG
CCCGCTGGCTCTGTGCAGAATCCTGTCTATCACAATCAGCCTCTGAACCCCGCGCCCAGC
AGAGACCCACACTACCAGGACCCCCACAGCACTGCAGTGGGCAACCCCGAGTATCTCAAC
ACTGTCCAGCCCACCTGTGTCAACAGCACATTCGACAGCCCTGCCCACTGGGCCCAGAAA
GGCAGCCACCAAATTAGCCTGGACAACCCTGACTACCAGCAGGACTTCTTTCCCAAGGAA
GCCAAGCCAAATGGCATCTTTAAGGGCTCCACAGCTGAAAATGCAGAATACCTAAGGGTC
GCGCCACAAAGCAGTGAATTTATTGGAGCATGA
```

(to be continued...)
