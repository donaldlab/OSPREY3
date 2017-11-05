/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.ewakstar;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * This writes the partition function EWAK* output to a file.
 * Based on Mark's 'SpecialWarningLog.java'
 * 
 * @author Anna Lowegard
 */
public class EWAKOutput {
    
    private String fileName;
    private BufferedWriter fileHandle;
    
    
    public EWAKOutput(String fileName){
        this.fileName = fileName;
   
        try {
            fileHandle = new BufferedWriter(new FileWriter(fileName));
        }
        catch(IOException e){
            throw new RuntimeException("ERROR opening EWAK* output file.  File name: "+fileName);
        }
    }
    
    
    public void write(String ewakOutput){
        try {
            fileHandle.write(ewakOutput);
        }
        catch(IOException e){
            throw new RuntimeException("ERROR writing to EWAK* output file.  File name: "+fileName);
        }
    }
    
    public void close(){
        try {
            fileHandle.close();
        }
        catch(IOException e){
            throw new RuntimeException("ERROR closing EWAK* output file.  File name: "+fileName);
        }
    }
    
    
    public void flush(){
        try{
            fileHandle.flush();
        }
        catch(IOException e){
            e.printStackTrace();
            throw new RuntimeException(e.getMessage());
        }
    }
            
    
}
