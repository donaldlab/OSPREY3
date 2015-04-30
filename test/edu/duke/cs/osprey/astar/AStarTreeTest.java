/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.cs.osprey.astar;

import java.util.ArrayList;
import junit.framework.TestCase;

/**
 *
 * @author mhall44
 */
public class AStarTreeTest extends TestCase {
    
    public AStarTreeTest(String testName) {
        super(testName);
    }
    
    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }
    
    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Test of nextConf method, of class AStarTree.
     */
    public void testNextConf() {
        System.out.println("nextConf");
        AStarTree instance = new AStarTreeImpl();
        int[] expResult = null;
        int[] result = instance.nextConf();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of canPruneNode method, of class AStarTree.
     */
    public void testCanPruneNode() {
        System.out.println("canPruneNode");
        AStarNode node = null;
        AStarTree instance = new AStarTreeImpl();
        boolean expResult = false;
        boolean result = instance.canPruneNode(node);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of outputNode method, of class AStarTree.
     */
    public void testOutputNode() {
        System.out.println("outputNode");
        AStarNode node = null;
        AStarTree instance = new AStarTreeImpl();
        int[] expResult = null;
        int[] result = instance.outputNode(node);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of refineScore method, of class AStarTree.
     */
    public void testRefineScore() {
        System.out.println("refineScore");
        AStarNode node = null;
        AStarTree instance = new AStarTreeImpl();
        instance.refineScore(node);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getChildren method, of class AStarTree.
     */
    public void testGetChildren() {
        System.out.println("getChildren");
        AStarNode curNode = null;
        AStarTree instance = new AStarTreeImpl();
        ArrayList expResult = null;
        ArrayList result = instance.getChildren(curNode);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of rootNode method, of class AStarTree.
     */
    public void testRootNode() {
        System.out.println("rootNode");
        AStarTree instance = new AStarTreeImpl();
        AStarNode expResult = null;
        AStarNode result = instance.rootNode();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isFullyAssigned method, of class AStarTree.
     */
    public void testIsFullyAssigned() {
        System.out.println("isFullyAssigned");
        AStarNode node = null;
        AStarTree instance = new AStarTreeImpl();
        boolean expResult = false;
        boolean result = instance.isFullyAssigned(node);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    public class AStarTreeImpl extends AStarTree {

        public ArrayList<AStarNode> getChildren(AStarNode curNode) {
            return null;
        }

        public AStarNode rootNode() {
            return null;
        }

        public boolean isFullyAssigned(AStarNode node) {
            return false;
        }
    }
}
