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
