## This file is part of OSPREY 3.0
## 
## OSPREY Protein Redesign Software Version 3.0
## Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
## 
## OSPREY is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License version 2
## as published by the Free Software Foundation.
## 
## You should have received a copy of the GNU General Public License
## along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
## 
## OSPREY relies on grants for its development, and since visibility
## in the scientific literature is essential for our success, we
## ask that users of OSPREY cite our papers. See the CITING_OSPREY
## document in this distribution for more information.
## 
## Contact Info:
##    Bruce Donald
##    Duke University
##    Department of Computer Science
##    Levine Science Research Center (LSRC)
##    Durham
##    NC 27708-0129
##    USA
##    e-mail: www.cs.duke.edu/brd/
## 
## <signature of Bruce Donald>, Mar 1, 2018
## Bruce Donald, Professor of Computer Science

import numpy as np
from gtda.homology import VietorisRipsPersistence

from preprocessing import get_mol2_coordinates_by_element, get_pdb_coordinates_by_element

def atom_persistence_homology(coords):
    # Track connected components, loops, and voids
    homology_dimensions = [0, 1, 2]

    # Collapse edges to speed up H2 persistence calculation!
    persistence = VietorisRipsPersistence(
        homology_dimensions=homology_dimensions,
        collapse_edges=True,
        max_edge_length=200
    )

    diagrams_basic = persistence.fit_transform(coords[None, :, :])

    return diagrams_basic

def opposition_homology(protein_coords, ligand_coords, remove_atoms_beyond=15):
    def opposition_distance_metric(vec1, vec2):
        if np.abs(vec1[-1] - vec2[-1]) > 2:  # If the two atoms are not of the same type
            return np.linalg.norm(vec1[:3] - vec2[:3])
        else:
            return np.Inf

    diff = protein_coords[:, np.newaxis, :] - ligand_coords[np.newaxis, :, :]
    distances_squared = np.sum(diff**2, axis=-1)
    mask = np.any(distances_squared < remove_atoms_beyond**2, axis=1)
    filtered_protein_coords = protein_coords[mask]
    # print(f'Number of protein atoms pruned out: {len(protein_coords) - len(filtered_protein_coords)}. Total number of protein atoms: {len(protein_coords)}')

    protein_coords = filtered_protein_coords

    # Append each coordinate with 1 for protein and 2 for ligand
    protein_coords = np.concatenate((protein_coords, np.ones((len(protein_coords), 1))), axis=1)
    ligand_coords = np.concatenate((ligand_coords, 4 * np.ones((len(ligand_coords), 1))), axis=1)

    combined_coords = np.concatenate((protein_coords, ligand_coords), axis=0)

    if combined_coords.shape[0] == 0:
        return None

    if protein_coords is None or ligand_coords is None:
        return None

    # Track connected components, loops, and voids
    homology_dimensions = [0, 1, 2]

    # Collapse edges to speed up H2 persistence calculation!
    persistence = VietorisRipsPersistence(
        metric=opposition_distance_metric,
        homology_dimensions=homology_dimensions,
        collapse_edges=True,
        max_edge_length=200
    )

    diagrams_basic = persistence.fit_transform(combined_coords[None, :, :])

    return diagrams_basic

def get_pairwise_opposition_persistence_diagrams(pdb_file, mol2_file, prune_beyond=15):
    protein_heavy_elements = ['C', 'N', 'O', 'S']
    ligand_heavy_elements = ['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'I']
    diagrams = []

    for pe in protein_heavy_elements:
        for le in ligand_heavy_elements:
            protein_coords = get_pdb_coordinates_by_element(pdb_file, pe)
            ligand_coords = get_mol2_coordinates_by_element(mol2_file, le)
            diagram = opposition_homology(protein_coords, ligand_coords, remove_atoms_beyond=prune_beyond)
            diagrams.append(diagram)

    return diagrams

def get_2345_persistence_diagrams(pdb_file, mol2_file):
    # The 2345 atom subsets in TNet-BP correspond to atom subsets 36-39 in PATH
    # But for inference we do not need 2345 persistence since persistence fingerprint any of these subsets
    return [[], [], [], []]