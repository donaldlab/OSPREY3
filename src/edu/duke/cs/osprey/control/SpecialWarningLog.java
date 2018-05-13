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

