/*
** This file is part of OSPREY.
** 
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 2 of the License, or
** (at your option) any later version.
** 
** OSPREY is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
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

