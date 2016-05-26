/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.control;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * This writes some particular kind of warning to a special file
 * e.g., residues that are deleted because we can't find a template for them
 * 
 * @author mhall44
 */
public class SpecialWarningLog {
    
    private String fileName;
    private BufferedWriter fileHandle;
    
    
    public SpecialWarningLog(String fileName){
        this.fileName = fileName;
   
        try {
            fileHandle = new BufferedWriter(new FileWriter(fileName));
        }
        catch(IOException e){
            throw new RuntimeException("ERROR opening special warning log.  File name: "+fileName);
        }
    }
    
    
    public void write(String warning){
        try {
            fileHandle.write(warning);
        }
        catch(IOException e){
            throw new RuntimeException("ERROR writing to special warning log.  File name: "+fileName);
        }
    }
    
    public void close(){
        try {
            fileHandle.close();
        }
        catch(IOException e){
            throw new RuntimeException("ERROR closing special warning log.  File name: "+fileName);
        }
    }
            
    
}
