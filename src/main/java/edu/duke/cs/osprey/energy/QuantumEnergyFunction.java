/*
 ** This file is part of OSPREY 3.0
 **
 ** OSPREY Protein Redesign Software Version 3.0
 ** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
 **
 ** OSPREY is free software: you can redistribute it and/or modify
 ** it under the terms of the GNU General Public License version 2
 ** as published by the Free Software Foundation.
 **
 ** You should have received a copy of the GNU General Public License
 ** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
 **
 ** OSPREY relies on grants for its development, and since visibility
 ** in the scientific literature is essential for our success, we
 ** ask that users of OSPREY cite our papers. See the CITING_OSPREY
 ** document in this distribution for more information.
 **
 ** Contact Info:
 **    Bruce Donald
 **    Duke University
 **    Department of Computer Science
 **    Levine Science Research Center (LSRC)
 **    Durham
 **    NC 27708-0129
 **    USA
 **    e-mail: www.cs.duke.edu/brd/
 **
 ** <signature of Bruce Donald>, Mar 1, 2018
 ** Bruce Donald, Professor of Computer Science
 */

package edu.duke.cs.osprey.energy;

import java.io.BufferedReader;
import java.io.InputStreamReader;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;

/**
 *
 * An energy function that evaluates the Quantum Mechanical Energy of a molecule
 *
 * @author Hunter Stephens
 */
public class QuantumEnergyFunction implements EnergyFunction {

    Molecule m;


    public QuantumEnergyFunction(Molecule m) {
        this.m = m;
    }


    @Override
    public double getEnergy() {
        // allow future methods to be added
       return EE_GMFCC();
    }

    private double EE_GMFCC(){
        // fragment molecule if not already
        if(!this.m.fragmented){
            this.m.fragment();
        }

        // calculate embedded MFCC energy
        double E_mfcc = 0.0;
        double Ei_QM, Ei_int;
        for(int i=0; i < this.m.fragments.size(); i++){
            //TODO: calculate fragment QM self-energy
            Ei_QM = 0.0;

            //TODO: calculate fragment background interaction energy
            Ei_int = 0.0;

            E_mfcc = E_mfcc + Ei_QM + Ei_int;
        }

        //TODO: calculate generalized MFCC energy
        // i.e. non-neighboring residues

        //TODO: subtract out over-counted background potential
        return 0;
    }
}
