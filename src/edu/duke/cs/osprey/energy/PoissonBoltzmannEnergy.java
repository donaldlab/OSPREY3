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
 * An energy function that evaluates the Poisson-Boltzmann solvation energy for a molecule
 * 
 * @author mhall44
 */
public class PoissonBoltzmannEnergy implements EnergyFunction {

    Molecule m;
    
    static String delphiFolder = "OSPREY_delphi";//where the Delphi calculations will happen

    public static double constRT = 1.9891/1000.0 * 298.15;//RT in kcal/mol
    
    
    public PoissonBoltzmannEnergy(Molecule m) {
        this.m = m;
    }
    
    
    @Override
    public double getEnergy() {
        //save the molecule to a PDB file and run Delphi on it
        
        PDBIO.writeFile(m, delphiFolder + "/struct.pdb");
        double E = 0;
        
        try{
            Process p = Runtime.getRuntime().exec(delphiFolder+"/getE");
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            p.waitFor();//We need to wait for getE to finish before going on and trying to read from it...
            
            E = Double.valueOf( br.readLine().trim() );
            E *= constRT;//convert from thermal energies (Delphi unit) to kcal/mol (OSPREY unit)
            br.close();
        }
        catch(Exception ex){
            throw new Error("can't compute Poisson-Boltzmann energy", ex);
        }

        return E;    
    }
    
    
    
}
