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

import java.io.*;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.XYZIO;

/**
 *
 * This class is for calculating the free energy of a molecule using the following methods
 * 1) Density Functional Theory
 *
 *
 * Currently this relies on the Psi4 library and must be visible in the path
 * Installing Psi4 using conda automatically includes the psi4 executable into
 * your path
 *
 * @author Hunter Stephens
 */
public class QuantumChemicalEnergy implements EnergyFunction {
    Molecule mol;

    public static double const_Hartrees_to_kcal_per_mol = 627.5;


    public QuantumChemicalEnergy(Molecule mol) {
        this.mol = mol;
    }

    @Override
    public double getEnergy()  {

        //get xyz format from molecule
        String mol_xyz = XYZIO.write(mol);
        //create psi4 input file
        StringBuilder in_file = new StringBuilder();
        in_file.append("molecule {\n");
        in_file.append(mol_xyz);
        in_file.append("}\n\n");
        in_file.append("set basis cc-pvdz\n\n");
        in_file.append("energy('b3lyp')\n");
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter("input.dat"));
            out.write(in_file.toString());
            out.close();
        }
        catch (IOException e) {
        }


        double E = 0;

        try{
            Process p = Runtime.getRuntime().exec("psi4 input.dat output.dat");
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            p.waitFor();//We need to wait for getE to finish before going on and trying to read from it...

            //E = Double.valueOf( br.readLine().trim() );
            //E *= constRT;//convert from thermal energies (Delphi unit) to kcal/mol (OSPREY unit)
            br.close();
        }
        catch(Exception ex){
            throw new Error("can't compute DFT energy", ex);
        }

        return E;
    }
}
