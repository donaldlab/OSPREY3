package edu.duke.cs.osprey.markstar;

import edu.duke.cs.osprey.markstar.visualizer.KStarTreeNode;
import edu.duke.cs.osprey.markstar.visualizer.SeqTreeNode;
import org.junit.Test;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class TestVisualizer {
    @Test
    public void testGetTreeFringe(){
        KStarTreeNode root = SeqTreeNode.parseTree("testSeqTree.tree");
        List<KStarTreeNode> fringe = KStarTreeNode.getFringeFromRoot(root);
        System.out.printf("There were %d fringe nodes%n", fringe.size());
    }
    @Test
    public void testOrderingFromFringe(){
        //KStarTreeNode root = SeqTreeNode.parseTree("testSeqTree.tree");
        KStarTreeNode root = SeqTreeNode.parseTree("/home/graham/Documents/sharkstar_profiling/_trees/2rl0_B_4.931E+10seqTree.tree");
        List<KStarTreeNode> fringe = KStarTreeNode.getFringeFromRoot(root);
        int[] ordering = KStarTreeNode.determineOrderingFromFringe(fringe);
        System.out.printf("Ordering: %s", Arrays.toString(ordering));
    }

    @Test
    public void testBuildTreeFromFringe() throws IOException {
        FileWriter writer = new FileWriter("testSeqTreeMod.tree");
        //KStarTreeNode root = SeqTreeNode.parseTree("testSeqTree.tree");
        SeqTreeNode root = SeqTreeNode.parseTree("/home/graham/Documents/sharkstar_profiling/_trees/2p4a_B_2.409E+11seqTree.tree");
        List<KStarTreeNode> fringe = KStarTreeNode.getFringeFromRoot(root);
        System.out.printf("There were %d fringe nodes%n", fringe.size());
        SeqTreeNode newRoot = (SeqTreeNode) KStarTreeNode.buildTreeFromFringe(fringe);
        newRoot.printTreeLikeMARKStar(writer, "");
        writer.flush();
        writer.close();
    }
}
