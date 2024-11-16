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

# Constants
R = 8.3144621 # J/(mol*K)
T = 298.15 # K

def delta_G_to_pKd(delta_G):
    """ delta_G in kcal/mol, pKd in M
    """
    return -delta_G * 1000 / (R * T * np.log(10)) * 4.184

def pKd_to_delta_G(pKd):
    """ pKd in M, delta_G in kcal/mol
    """
    return -pKd * R * T * np.log(10) / 1000 / 4.184
