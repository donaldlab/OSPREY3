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
