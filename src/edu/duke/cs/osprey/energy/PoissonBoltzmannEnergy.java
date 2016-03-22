/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy;

import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileWriter;
import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 *
 * An energy function that evaluates the Poisson-Boltzmann solvation energy for a molecule
 * 
 * @author mhall44
 */
public class PoissonBoltzmannEnergy implements EnergyFunction {

    Molecule m;
    
    static String delphiFolder = "OSPREY_delphi";//where the Delphi calculations will happen

    static double constRT = 1.9891/1000.0 * 298.15;//RT in kcal/mol
    
    
    public PoissonBoltzmannEnergy(Molecule m) {
        this.m = m;
    }
    
    
    @Override
    public double getEnergy() {
        //save the molecule to a PDB file and run Delphi on it
        
        PDBFileWriter.writePDBFile(m, delphiFolder+"/struct.pdb");
        double E = 0;
        
        try{
            Process p = Runtime.getRuntime().exec(delphiFolder+"/getE");
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            p.waitFor();//We need to wait for getE to finish before going on and trying to read from it...
            
            E = Double.valueOf( br.readLine().trim() );
            E *= constRT;//convert from thermal energies (Delphi unit) to kcal/mol (OSPREY unit)
            br.close();
        }
        catch(Exception e){
            System.err.println(e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        return E;    
    }
    
    
    
}