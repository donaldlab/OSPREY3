/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

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
 */
public class UnitTestSuite {
    //Suite of all unit tests to run
    //can be good to check that we didn't break something
    //this suite is meant to run in the directory 1CC8

    public static void runAllTests() {

        //Make sure we're in the right directory
        try {
            FileInputStream is = new FileInputStream("1CC8.ss.pdb");
            is.close();
        } catch (FileNotFoundException e) {
            throw new RuntimeException("ERROR: Tests need to be run in the directory test/1CC8");
        } catch (IOException e) {//this is weird
            throw new RuntimeException(e.getMessage());
        }

         EnergyTests.test1CC8Energy();
        
         DOFTests.testMutation();
         DOFTests.testDihedral();
        
         ToolTests.SuperposingRotMatrixTest();
         ToolTests.RigidMotionTest();
        
         ConfSpaceTests.testConfSpaceGeneration();
         PDBTests.testPDBReadWrite();
        
         ConfSearchTests.testDEE(true);
         ConfSearchTests.testDEE(false);
         ConfSearchTests.testExhaustive(false, false);
         ConfSearchTests.testExhaustive(false, true);

        ConfSearchTests.testExhaustive(true, false);
    }
    
    public static void main(String[] args)
    {
        runAllTests();
    }
    
}
