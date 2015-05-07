/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.tests;

/**
 *
 * @author mhall44
 */
public class UnitTestSuite {
    //Suite of all unit tests to run
    //can be good to check that we didn't break something
    //this suite is meant to run in the directory 1CC8
    
    public static void runAllTests(){
        
        EnergyTests.test1CC8Energy();
        
        DOFTests.testMutation();
        DOFTests.testDihedral();
        
        ToolTests.SuperposingRotMatrixTest();
        ToolTests.RigidMotionTest();
        
        ConfSpaceTests.testConfSpaceGeneration();
        PDBTests.testPDBReadWrite();
        
        ConfSearchTests.testDEE();
        ConfSearchTests.testExhaustive();
    }
    
}
