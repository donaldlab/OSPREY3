package edu.duke.cs.osprey.sharkstar.tools;

import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Collections;
import java.util.List;

public class MultiSequenceSHARKStarNodeStatistics {

    public static void printTree(String prefix, FileWriter writer, SimpleConfSpace confSpace, Sequence seq, MultiSequenceSHARKStarNode node)
    {
        MultiSequenceSHARKStarNode.Node confSearchNode = node.getConfSearchNode();
        String confString = confSearchNode.confToString();
        if(confSpace != null)
            confString = confString+"->("+confSpace.formatConfRotamersWithResidueNumbers(confSearchNode.assignments)+")";
        String out = prefix+confString+":"
                +"["+node.getConfLowerBound(seq)+","+node.getConfUpperBound(seq)+"]->"
                +"["+formatBound(node.getLowerBound(seq))
                +","+formatBound(node.getUpperBound(seq))+"]"+"\n";
        if(MathTools.isLessThan(node.getUpperBound(seq), BigDecimal.ZERO))
            return;
        if(writer != null) {
            try {
                writer.write(out);
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        else
            System.out.print(out);
        List<MultiSequenceSHARKStarNode> children = node.getChildren(seq);
        if( children != null && ! children.isEmpty()) {
            Collections.sort( children, (a, b)-> -a.getSequenceBounds(seq).upper
                    .compareTo(b.getSequenceBounds(seq).upper));
            for (MultiSequenceSHARKStarNode child :  children)
                printTree(prefix + "~+", writer, confSpace, seq, child);
        }


    }

    private static String formatBound(BigDecimal bound) {
        if(bound instanceof MathTools.MagicBigDecimal)
            return bound.toString();
        return setSigFigs(bound).toString();
    }

    public static void printTree(String name, Sequence seq, MultiSequenceSHARKStarNode node) {
        printTree(name, null, seq, node);
    }

    public static void printTree(String name, SimpleConfSpace confspace, Sequence seq, MultiSequenceSHARKStarNode node)
    {
        try {
            FileWriter writer = new FileWriter(new File(name+"ConfTreeBounds.txt"));
            printTree("",  writer, confspace, seq, node);
            writer.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void printTree(Sequence seq, MultiSequenceSHARKStarNode node) {
        printTree("", null, null, seq, node);
    }

    public static BigDecimal setSigFigs(BigDecimal decimal, int numSigFigs) {
        return decimal.setScale(4-decimal.precision()+decimal.scale(), RoundingMode.HALF_UP);
    }

    public static BigDecimal setSigFigs(BigDecimal decimal){
        return setSigFigs(decimal, 4);
    }
    public static void printBoundBreakDown(Sequence seq, MultiSequenceSHARKStarNode node)
    {
        printBoundBreakDown(seq, "", node);
    }

    public static void printBoundBreakDown(Sequence seq, String prefix, MultiSequenceSHARKStarNode node)
    {
        MultiSequenceSHARKStarNode.Node confSearchNode = node.getConfSearchNode();
        if(node.level == 0) {
            System.out.println("=====================BEGIN TREE INFO==================================");
            System.out.println(prefix + confSearchNode + ": [" + setSigFigs(node.getSequenceBounds(seq).lower)
                    + "," + setSigFigs(node.getSequenceBounds(seq).upper) + "], errorBound =" + String.format("%3.3e",node.getErrorBound(seq)));
        }

        List< MultiSequenceSHARKStarNode> childrenForSequence = node.getChildren(seq);

        if( childrenForSequence  != null &&  childrenForSequence .size() > 0) {
            BigDecimal upper = BigDecimal.ZERO;
            BigDecimal lower = BigDecimal.ZERO;
            Collections.sort( childrenForSequence );
            prefix+="+~~";
            for(MultiSequenceSHARKStarNode child:  childrenForSequence ) {
                BigDecimal childUpper = child.getSequenceBounds(seq).upper;
                BigDecimal childLower = child.getSequenceBounds(seq).lower;
                System.out.print(prefix+child.getConfSearchNode()+": ["+setSigFigs(childLower)
                        +","+setSigFigs(childUpper)+"], epsilon="+String.format("%3.3e",node.getErrorBound(seq)));
                System.out.print("Upper: " + setSigFigs(upper) + " + "
                        + setSigFigs(childUpper) + " = "
                        + setSigFigs(upper.add(childUpper)));
                System.out.println("Lower: " + setSigFigs(lower) + " + "
                        + setSigFigs(childLower) + " = "
                        + setSigFigs(lower.add(childLower)));
                upper = upper.add(childUpper);
                lower = lower.add(childLower);
                printBoundBreakDown(seq, prefix, child);
            }
        }
        if(node.level == 0) {
            System.out.println("=====================END TREE INFO==================================");
        }
    }
}
