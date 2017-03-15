This is the development version of OSPREY.  For stable versions, please visit http://www.cs.duke.edu/donaldlab/osprey.php


OSPREY Protein Redesign Software
Copyright (C) 2001-2017 Bruce Donald Lab, Duke University

OSPREY is free software: you can redistribute it and/or modify it under the terms of the GNU
Lesser General Public License as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version. The two parts of the license are attached below
(Sec. 3). There are additional restrictions imposed on the use and distribution of this open-source
code, including:
• The header from Sec. 1 must be included in any modification or extension of the code;
• Any publications, grant applications, or patents that use OSPREY must state that OSPREY
was used, with a sentence such as ”We used the open-source OSPREY software [Ref] to
design....”
• Any publications, grant applications, or patents that use OSPREY must cite our papers. The
citations for the various different modules of our software are described in Sec. 2.


Section 1: Source Header

OSPREY Protein Redesign Software
Copyright (C) 2001-2017 Bruce Donald Lab, Duke University
OSPREY is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.
OSPREY is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, see:
<http://www.gnu.org/licenses/>.
There are additional restrictions imposed on the use and distribution
of this open-source code, including: (A) this header must be included
in any modification or extension of the code; (B) you are required to
cite our papers in any publications that use this code. The citation
for the various different modules of our software, together with a
complete list of requirements and restrictions are found in the
document license.pdf enclosed with this distribution.
Contact Info:
Bruce Donald
Duke University
Department of Computer Science
Levine Science Research Center (LSRC)
Durham
NC 27708-0129
USA
e-mail: www.cs.duke.edu/brd/
<signature of Bruce Donald>, Mar 1, 2012
Bruce Donald, Professor of Computer Science



Section 2: Citation Requirements

The citation requirements for the various different modules of our software are:
• For a general citation, please use:
C. Chen, I. Georgiev, A. C. Anderson, and B. R. Donald. Computational structure-based redesign
of enzyme activity. PNAS USA, 106(10):3764–3769, 2009.
P. Gainza, K. E. Roberts, I. Georgiev, R. H. Lilien, D. A. Keedy, C. Chen, F. Reza, A. C. Anderson,
D. C. Richardson, J. S. Richardson, and B. R. Donald. OSPREY: Protein design with ensembles,
flexibility, and provable algorithms. Methods in Enzymology, 523:87–107, 2012.
• iMinDEE:
P. Gainza, K.E. Roberts, and B.R. Donald. Protein Design using Continuous Rotamers PLoS
Computational Biology, (1): e1002335. doi:10.1371/journal.pcbi.1002335, 2012.
• Protein:Protein Interactions:
K.E. Roberts, P.R. Kushing, P. Boisguerin, DR Madden, and B.R. Donald. Design of protein-protein interactions with a novel ensemble-based scoring algorithm. Research in Computational
Molecular Biology (RECOMB), volume 6577 of Lecture Notes in Computer Science. Heidelberg: Springer
Berlin. pp. 361376. 2011
K.E. Roberts, P.R. Kushing, P. Boisguerin, DR Madden, and B.R. Donald. Computational Design
of a PDZ Domain Peptide Inhibitor that Rescues CFTR Activity. PLoS Computational Biology, 8(4):e1002477, 2012. 
• K∗ (current implementation):
I. Georgiev, R. Lilien, and B. R. Donald. The minimized dead-end elimination criterion and its
application to protein redesign in a hybrid scoring and search algorithm for computing partition
functions over molecular ensembles. J Comput Chem, 29(10):1527–42, 2008.
To cite the general idea of K*, you can also cite:
R. Lilien, B. Stevens, A. Anderson, and B. R. Donald. A novel ensemble-based scoring and
search algorithm for protein redesign, and its application to modify the substrate specificity of
the Gramicidin Syntheses A phenylalanine adenylation enzyme. J Comp Biol, 12(6–7):740–761,
2005.
• DEEPer:
M. A. Hallen, D. A. Keedy, and B. R. Donald. Dead-end elimination with perturbations (DEEPer):
A provable protein design algorithm with continuous sidechain and backbone flexibility. Proteins, 81(1):18–39, 2013.
• COMETS:
M. A. Hallen and B. R. Donald. COMETS (Constrained Optimization of Multistate Energies by Tree Search): A provable and efficient protein design algorithm to optimize binding affinity and specificity with respect to sequence. Journal of Computational Biology, 23(5):311–321, 2015.
• EPIC:
M. A. Hallen, P. Gainza, and B. R. Donald. A compact representation of continuous energy surfaces for more efficient protein design. Journal of Chemical Theory and Computation, 11(5):2292–2306, 2015.
• Dynamic A*:
K. E. Roberts, P. Gainza, M. A. Hallen, and B. R. Donald. Fast gap-free enumeration of conformations and sequences for protein design. Proteins, 83(10):1859–1877, 2015.
• LUTE:
M. A. Hallen, J. D. Jou, and B. R. Donald. LUTE (Local Unpruned Tuple Expansion): Accurate continuously flexible protein design with general energy functions and rigid-rotamer-like efficiency. In International Conference on Research in Computational Molecular Biology (RECOMB), pages 122–136. Springer, 2016.
• CATS:
M. A. Hallen and B. R. Donald. CATS (Coordinates of Atoms by Taylor Series): Protein design with backbone flexibility in all locally feasible directions. In Intelligent Systems for Molecular Biology (ISMB), in press, 2017.  
