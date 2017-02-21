/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.energy;

import java.util.ArrayList;

import edu.duke.cs.osprey.confspace.Strand;
import edu.duke.cs.osprey.control.EnvironmentVars;
import edu.duke.cs.osprey.control.ParamSet;
import edu.duke.cs.osprey.structure.Molecule;
import edu.duke.cs.osprey.structure.PDBIO;
import edu.duke.cs.osprey.structure.Residue;

/**
 *
 * For each residue in a ligand, calculate its interactions with the target
 * 
 * @author mhall44
 */
public class LigandResEnergies {
    
    ArrayList<Residue> ligandRes;
    ArrayList<Residue> targetRes;
    
    
    public LigandResEnergies(ParamSet params){
        Molecule mol = new Strand.Builder(PDBIO.readFile(params.getFile("pdbName"))).build().mol;
        String ligandTermini[] = new String[] {params.getValue("ligandStart"),params.getValue("ligandEnd")};
        String targetTermini[] = new String[] {params.getValue("targetStart"),params.getValue("targetEnd")};
        ligandRes = mol.resListFromTermini(ligandTermini, null);
        targetRes = mol.resListFromTermini(targetTermini, null);
    }
    
    
    public void printEnergies(){
        System.out.println("PRINTING LIGAND RESIDUE ENERGIES");
        System.out.println("RES ENERGY");
        
        for(Residue res : ligandRes){
            double E = calcLigandResEnergy(res);
            System.out.println(res.getPDBResNumber() + " " + E);
        }
    }
    
    
    public double calcLigandResEnergy(Residue res){
        //calculate interactions of res (assumed to be in the ligand)
        //with the target
        double E = 0;
        for(Residue tres : targetRes){
            EnergyFunction pairEFunc = EnvironmentVars.curEFcnGenerator.resPairEnergy(res, tres);
            double pairE = pairEFunc.getEnergy();
            E += pairE;
        }
        
        return E;
    }
    
}
