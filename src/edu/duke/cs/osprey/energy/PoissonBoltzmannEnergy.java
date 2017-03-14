/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
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