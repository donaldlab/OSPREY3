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

from pathlib import Path
import pickle
import time

import numpy as np
from functional import seq
from gudhi.representations.vector_methods import PersistenceImage

import ph
from unit_conversion import pKd_to_delta_G

# Subset of IPC for persistence fingerprint
tiny_feature_subset_indices = [17819,  18318,  18514,  18908, 288020, 288516,  31500,  31702, 61301,  91001]

fingerprint_sources = list(map(lambda x: (x // 30000, (x % 30000) // 10000), tiny_feature_subset_indices))  # Each element of this list is tuple(diagram index: 0-39, image index i.e. homology dimension within trio: 0-2)

optimal_mini_gbr_path = Path(__file__).parent / 'regrs/optimal_mini_gbr.pkl'
with open(optimal_mini_gbr_path, 'rb') as f:
    optimal_mini_gbr = pickle.load(f)

def diagram_to_image(diagram, diagram_index=None):
    if diagram is None:
        return np.zeros((3, 10000))
    else:
        if diagram_index is not None:
            trio_indices_needed = seq(fingerprint_sources).filter(lambda x: x[0] == diagram_index).map(lambda x: x[1]).to_list()
        else:
            trio_indices_needed = [0, 1, 2]

        pim = PersistenceImage(bandwidth=0.2355, resolution=[100,100], im_range=[0, 50, 0, 50])

        # Now reshape into (3, n, 2)
        def get_nd_diagram(diagram, dim):
            filtered_diagrams = list(filter(lambda x: x[2] == dim, diagram[0, :, :]))
            return np.array(filtered_diagrams)[:, :2]

        gudhi_format_diagrams = [get_nd_diagram(diagram, 0), get_nd_diagram(diagram, 1), get_nd_diagram(diagram, 2)]

        diagrams_clipped = [np.clip(diagram, 0, 50) for diagram in gudhi_format_diagrams]

        imgs = []
        for trio_index, diagram in enumerate(diagrams_clipped):
            # For all images that are not used, just return a zero image
            if trio_index not in trio_indices_needed:
                imgs.append(np.zeros(10000))
            else:
                imgs.append(pim.fit_transform(diagram[None, :, :])[0])
        imgs = np.array(imgs)

        return imgs


# Takes protein file, ligand file, and predicts
def get_all_images(protein_file, ligand_file):
    # do persistent homology on protein_file and ligand_file
    pw_opposition_diagrams = ph.get_pairwise_opposition_persistence_diagrams(protein_file, ligand_file)
    other_persistence_diagrams = ph.get_2345_persistence_diagrams(protein_file, ligand_file)

    all_diagrams = pw_opposition_diagrams + other_persistence_diagrams

    all_images = []
    # diagram_index has range [0, 40]
    for diagram_index, diagram in enumerate(all_diagrams):
        # If the diagram is entirely not used, just return the trio of zero images
        if diagram_index not in list(map(lambda x:x[0], fingerprint_sources)):
            all_images.append(diagram_to_image(None))
        else:
            all_images.append(diagram_to_image(diagram, diagram_index=diagram_index))

    return all_images

def predict(protein_file, ligand_file):
    begin_time = time.time()

    print('Computing persistence images')
    all_images = get_all_images(protein_file, ligand_file)

    print('Predicting')
    observations = np.array(all_images).flatten().reshape(1, -1)
    persistence_fingerprint = observations[:, tiny_feature_subset_indices]
    prediction_path = optimal_mini_gbr.predict(persistence_fingerprint)[0]  # Result is in pKd, need to convert to delta G
    time_taken = time.time() - begin_time

    return {
        'persistence_fingerprint': persistence_fingerprint[0],
        'prediction (delta_G, kcal/mol)': pKd_to_delta_G(prediction_path),
        'time_taken (s)': time_taken,
    }

if __name__ == '__main__':
    protein_file = 'example/1a1e_protein.pdb'
    ligand_file = 'example/1a1e_ligand.mol2'
    res = predict(protein_file, ligand_file)
    print(res)