/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import edu.duke.cs.osprey.confspace.SearchProblem;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * Prints information about a conformation
 * 
 * @author mhall44
 */
public class ConfPrinter {
    
    SearchProblem searchSpace;
    BufferedWriter confFileHandle;
    String confFileName;
    
    boolean printEPICEnergy;
    
    ConfPrinter(SearchProblem searchProb, String confFileName, boolean printEPICEnergy){
        //open (for writing) a file to record conformations in
        searchSpace = searchProb;
        this.confFileName = confFileName;
        this.printEPICEnergy = printEPICEnergy;
        
        try {
            confFileHandle = new BufferedWriter(new FileWriter(confFileName));
        }
        catch(Exception e){
            throw new RuntimeException("ERROR OPENING CONF FILE.  NAME: "
                    + confFileName + e.getMessage());
        }
    }
    
    
    void printConf(int[] conf, double confE, double lowerBound, double bestESoFar, 
            int confCount ){
        
        try {
            System.out.println("ENUMERATING CONFORMATION.  RCs (residue-based numbers):");
            confFileHandle.write(confCount + " CONF: ");
            for(int rc : conf){
                System.out.print(rc + " ");
                confFileHandle.write(rc + " ");
            }
            System.out.println();


            System.out.println("Residue types: ");
            confFileHandle.write("RESTYPES: ");
            for(int pos=0; pos<searchSpace.confSpace.numPos; pos++){
                String resType = searchSpace.confSpace.posFlex.get(pos).RCs.get(conf[pos]).AAType;
                System.out.print( resType + " " );
                confFileHandle.write(resType + " ");
            }
            System.out.println();


            System.out.println("Rotamer numbers: ");
            confFileHandle.write("ROTS: ");
            for(int pos=0; pos<searchSpace.confSpace.numPos; pos++){
                int rotNum = searchSpace.confSpace.posFlex.get(pos).RCs.get(conf[pos]).rotNum;
                System.out.print( rotNum + " " );
                confFileHandle.write( rotNum + " " );
            }
            System.out.println();


            String energyStatement = "Lower bound/enumeration energy: "+lowerBound+" Energy: "+confE+" Best so far: "+bestESoFar;
            //Lower bound/enumeration energy is what we enumerate in order of
            //(either a lower bound on the actual energy, or the same as Energy)

            if(printEPICEnergy)//useful to see EPIC energy (confE is regular E, lowerBound is tup-exp)
                energyStatement = energyStatement + " EPIC energy: " + searchSpace.EPICMinimizedEnergy(conf);

            System.out.println(energyStatement);
            
            confFileHandle.write(energyStatement);
            confFileHandle.newLine();
        }
        catch(IOException e){
            e.printStackTrace();
            throw new RuntimeException(e.getMessage());
        }
    }

    
    
        
    
    void closeConfFile(){
        //close it
        
        try {
            confFileHandle.close();
        }
        catch(Exception e){
            throw new RuntimeException("ERROR CLOSING CONF FILE.  NAME: "
                    + confFileName + e.getMessage());
        }
    }
}
