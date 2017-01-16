/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

/**
 *
 * @author mhall44
 */

/* TODO:
 * Get this running as a separate class. There's no reason to have to invoke it with all the other
 * parameters from the command line.

  MH 7/1/16: I have converted most of the tests that were here to JUnit tests and moved them
  to the tests directory that Jeff created.  The remaining ones create PDB files that
  you can look at...I left them here in case they might come in handy.  
 */
public class UnitTestSuite {
    //this suite is meant to run in the directory 1CC8
    
    public static void runAllTests(){
        
        //Make sure we're in the right directory
        try {
            FileInputStream is = new FileInputStream("1CC8.ss.pdb");
            is.close();
        }
        catch(FileNotFoundException e){
            throw new RuntimeException("ERROR: Tests need to be run in the directory examples/1CC8");
        }
        catch(IOException e){//this is weird
            throw new RuntimeException(e.getMessage());
        }
                
        DOFTests.testMutation();
        DOFTests.testDihedral();
    }
    
    
}
