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

import types

import pandas as pd
from Bio.PDB import PDBParser
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import numpy as np

def get_pdb_coordinates_by_element(file, element) -> np.ndarray:
    # Similar to get_pdb_coordinates, but only gets the coordinate matching the input element
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure(file, file)

    coords = []
    for atom in structure.get_atoms():
        if atom.element == element:
            coords.append(list(atom.get_vector()))

    if len(coords) == 0:
        return np.array([]).reshape(0, 3)

    coords = np.array(coords)
    coords = np.round(coords, 5)

    return coords

def get_mol2_coordinates(file):
    # Gets coordinates from a mol2 or pdb file
    file = str(file)
    if '.mol2' in file:
        pmol = PandasMol2().read_mol2(file)
    elif '.pdb' in file:
        pmol = types.SimpleNamespace()
        pmol.df = pd.concat([PandasPdb().read_pdb(file).df['ATOM'], PandasPdb().read_pdb(file).df['HETATM']])
        pmol.df['atom_type'] = pmol.df['atom_name'].str[0]
        pmol.df['x'] = pmol.df['x_coord']
        pmol.df['y'] = pmol.df['y_coord']
        pmol.df['z'] = pmol.df['z_coord']
    else:
        raise Exception('Invalid file type.')

    return pmol.df[['x', 'y', 'z']].to_numpy()

def get_mol2_coordinates_by_element(file, element) -> np.ndarray:
    # Gets coordinates from a mol2 or pdb file
    # Only selects the atoms of a specific element type
    file = str(file)
    if '.mol2' in file:
        pmol = PandasMol2().read_mol2(file)
    elif '.pdb' in file:
        pmol = types.SimpleNamespace()
        pmol.df = pd.concat([PandasPdb().read_pdb(file).df['ATOM'], PandasPdb().read_pdb(file).df['HETATM']])
        pmol.df['atom_type'] = pmol.df['atom_name'].str[0]
        pmol.df['x'] = pmol.df['x_coord']
        pmol.df['y'] = pmol.df['y_coord']
        pmol.df['z'] = pmol.df['z_coord']
    else:
        raise Exception('Invalid file type.')

    coords = []

    for idx in range(len(pmol.df)):
        sybyl_atom_type = pmol.df.loc[idx, 'atom_type']
        if sybyl_atom_type.split('.')[0] == element:
            coords.append(pmol.df.loc[idx, ['x', 'y', 'z']].to_numpy())

    if len(coords) == 0:
        return np.array([]).reshape(0, 3)

    coords = np.array(coords)
    return coords