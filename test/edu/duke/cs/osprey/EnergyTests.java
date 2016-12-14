/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.energy.EnergyFunction;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBFileReader;

/**
 *
 * @author mhall44
 */
public class EnergyTests {
    
    
    public static void test1CC8Energy(){
        //compute the full energy for 1CC8 using the default AMBER forcefield, and compare it to OSPREY 2 values 
        Molecule m = PDBFileReader.readPDBFile("1CC8.ss.pdb", null);
        EnergyFunction fullEFunc = EnvironmentVars.curEFcnGenerator.fullMolecEnergy(m);
        double fullE = fullEFunc.getEnergy();
        
        System.out.println("1CC8.ss.pdb full energy: "+fullE);
        
        //Compare to correct energies (assuming default settings), calculated with OSPREY 2
        double correctE;
        if(EnvironmentVars.curEFcnGenerator.ffParams.doSolvationE)
            correctE = -986.6491862981836;
        else
            correctE = -639.7025085949941;
        
        double error = Math.abs(fullE-correctE);
        System.out.println("1CC8.ss.pdb energy error: "+error);
    }
}
