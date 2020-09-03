package edu.duke.cs.osprey.sharkstar.tools;

import EDU.oswego.cs.dl.util.concurrent.FJTask;
import edu.duke.cs.osprey.astar.conf.RCs;
import edu.duke.cs.osprey.confspace.Sequence;
import edu.duke.cs.osprey.confspace.SimpleConfSpace;
import edu.duke.cs.osprey.kstar.pfunc.BoltzmannCalculator;
import edu.duke.cs.osprey.kstar.pfunc.PartitionFunction;
import edu.duke.cs.osprey.sharkstar.MultiSequenceSHARKStarNode;
import edu.duke.cs.osprey.sharkstar.SingleSequenceSHARKStarBound;
import edu.duke.cs.osprey.tools.MathTools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Collections;
import java.util.List;

public class MultiSequenceSHARKStarNodeStatistics {

    private static BoltzmannCalculator bc = new BoltzmannCalculator(PartitionFunction.decimalPrecision);

    public interface BoundGetter {
        MathTools.BigDecimalBounds getBounds(MultiSequenceSHARKStarNode node);
    }

    public static String rBD(BigDecimal bd) {
        return String.format("%12.6e", bd);
    }

    public static String treeString(String prefix, Sequence seq, SimpleConfSpace confSpace, MultiSequenceSHARKStarNode node) {
        RCs seqRCs = seq.makeRCs(confSpace);
        List<MultiSequenceSHARKStarNode> tempChildren = node.getAllChildren();
        int nextPos = -1;
        if(!tempChildren.isEmpty()){
            nextPos = tempChildren.get(0).pos;
        }
        BoundGetter boundGetter = (node1) -> node1.getSequenceBounds(seq);
        MathTools.BigDecimalBounds bounds = boundGetter.getBounds(node);
        String confString = node.confToString();
        String out = prefix+confString+":"
                +"["+node.getConfLowerBound(seq)+","+node.getConfUpperBound(seq)+"]->"
                +"["+formatBound(bounds.lower)
                +","+formatBound(bounds.upper)+"]"+"\n";
        if(MathTools.isLessThan(bc.calc(node.getConfLowerBound(seq)), BigDecimal.ZERO))
            return out;
        List<MultiSequenceSHARKStarNode> children = node.getChildren(seqRCs.get(nextPos));
        if( children != null && ! children.isEmpty()) {
            BoundGetter finalBoundGetter = boundGetter;
            Collections.sort( children, (a, b)->
                    -MathTools.compare(finalBoundGetter.getBounds(a).upper, finalBoundGetter.getBounds(b).upper));
            for (MultiSequenceSHARKStarNode child :  children)
                out += treeString(prefix + "~+", seq, confSpace, child);
        }
        return out;
    }

    public static void printTree(String prefix, FileWriter writer, SimpleConfSpace confSpace, Sequence seq,
                                 MultiSequenceSHARKStarNode node, BoundGetter boundGetter)
    {
        RCs seqRCs = seq.makeRCs(confSpace);
        List<MultiSequenceSHARKStarNode> tempChildren = node.getAllChildren();
        int nextPos = -1;
        if(!tempChildren.isEmpty()){
            nextPos = tempChildren.get(0).pos;
        }
        if(boundGetter == null)
            boundGetter = (node1) -> node1.getSequenceBounds(seq);
        /*
        if (MathTools.isLessThan(node.getUpperBound(seq), BigDecimal.ONE))
            return;
         */
        MathTools.BigDecimalBounds bounds = boundGetter.getBounds(node);
        String confString = node.confToString();
        if(confSpace != null)
            confString = confString+"->("+confSpace.formatConfRotamersWithResidueNumbers(node.assignments)+")";
        String out = prefix+confString+":"
                +"["+node.getConfLowerBound(seq)+","+node.getConfUpperBound(seq)+"]->"
                +"["+formatBound(bounds.lower)
                +","+formatBound(bounds.upper)+"]"+"\n";
        if(MathTools.isLessThan(bc.calc(node.getConfLowerBound(seq)), BigDecimal.ZERO))
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
        List<MultiSequenceSHARKStarNode> children = node.getChildren(seqRCs.get(nextPos));
        if( children != null && ! children.isEmpty()) {
            BoundGetter finalBoundGetter = boundGetter;
            Collections.sort( children, (a, b)->
                    -MathTools.compare(finalBoundGetter.getBounds(a).upper, finalBoundGetter.getBounds(b).upper));
            for (MultiSequenceSHARKStarNode child :  children)
                printTree(prefix + "~+", writer, confSpace, seq, child, boundGetter);
        }


    }

    private static String formatBound(BigDecimal bound) {
        if(bound instanceof MathTools.MagicBigDecimal)
            return convertMagicBigDecimalToString(bound);
        return convertMagicBigDecimalToString(bound);
    }

    public static void printTree(String name, Sequence seq, MultiSequenceSHARKStarNode node) {
        printTree(name, null, seq, node);
    }

    public static void printTree(String name, SimpleConfSpace confspace, Sequence seq, MultiSequenceSHARKStarNode node)
    {
        try {
            FileWriter writer = new FileWriter(new File(name+"ConfTreeBounds.txt"));
            printTree("",  writer, confspace, seq, node, null);
            writer.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void printTree(Sequence seq, MultiSequenceSHARKStarNode node) {
        printTree("",  null, null, seq, node, null);
    }

    public static BigDecimal setSigFigs(BigDecimal decimal, int numSigFigs) {
        return decimal.setScale(4-decimal.precision()+decimal.scale(), RoundingMode.HALF_UP);
    }

    public static String formatBounds(BigDecimal lowerBound, BigDecimal upperBound) {
        return "["+convertMagicBigDecimalToString(lowerBound)+","+convertMagicBigDecimalToString(upperBound)+"]";
    }

    public static String convertMagicBigDecimalToString(BigDecimal decimal) {
        if(MathTools.isFinite(decimal))
            if(MathTools.isLessThan(decimal, BigDecimal.ONE))
                return "(<1)";
            else return setSigFigs(decimal).toEngineeringString();
        return decimal.toString();
    }

    public static BigDecimal setSigFigs(BigDecimal decimal){
        return setSigFigs(decimal, 4);
    }
    public static void printBoundBreakDown(Sequence seq, SimpleConfSpace confSpace, MultiSequenceSHARKStarNode node)
    {
        printBoundBreakDown(seq, confSpace, "", node);
    }

    public static void printBoundBreakDown(Sequence seq, SimpleConfSpace confSpace, String prefix, MultiSequenceSHARKStarNode node)
    {
        if(node.level == 0) {
            System.out.println("=====================BEGIN TREE INFO==================================");
            System.out.println(prefix + node + ": [" + setSigFigs(node.getSequenceBounds(seq).lower)
                    + "," + setSigFigs(node.getSequenceBounds(seq).upper) + "], errorBound =" + String.format("%3.3e",node.getErrorBound(seq)));
        }

        RCs seqRCs = seq.makeRCs(confSpace);
        List<MultiSequenceSHARKStarNode> tempChildren = node.getAllChildren();
        int nextPos = -1;
        if(!tempChildren.isEmpty()){
            nextPos = tempChildren.get(0).pos;
        }
        List< MultiSequenceSHARKStarNode> childrenForSequence = node.getChildren(seqRCs.get(nextPos));

        if( childrenForSequence  != null &&  childrenForSequence .size() > 0) {
            BigDecimal upper = BigDecimal.ZERO;
            BigDecimal lower = BigDecimal.ZERO;
            //Collections.sort( childrenForSequence );
            prefix+="+~~";
            for(MultiSequenceSHARKStarNode child:  childrenForSequence ) {
                BigDecimal childUpper = child.getSequenceBounds(seq).upper;
                BigDecimal childLower = child.getSequenceBounds(seq).lower;
                System.out.print(prefix+child+": ["+setSigFigs(childLower)
                        +","+setSigFigs(childUpper)+"], epsilon="+String.format("%3.3e",node.getErrorBound(seq)));
                System.out.print("Upper: " + setSigFigs(upper) + " + "
                        + setSigFigs(childUpper) + " = "
                        + setSigFigs(upper.add(childUpper)));
                System.out.println("Lower: " + setSigFigs(lower) + " + "
                        + setSigFigs(childLower) + " = "
                        + setSigFigs(lower.add(childLower)));
                upper = upper.add(childUpper);
                lower = lower.add(childLower);
                printBoundBreakDown(seq, confSpace, prefix, child);
            }
        }
        if(node.level == 0) {
            System.out.println("=====================END TREE INFO==================================");
        }
    }

    public static void checkPartialConfNode(SingleSequenceSHARKStarBound bound, List<MultiSequenceSHARKStarNode> processedNodes) {
        for(MultiSequenceSHARKStarNode curNode: processedNodes) {
            List<MultiSequenceSHARKStarNode> tempChildren = curNode.getAllChildren();
            int nextPos = -1;
            if(!tempChildren.isEmpty()){
                nextPos = tempChildren.get(0).pos;
            }
            List<MultiSequenceSHARKStarNode> children = null;
            if(nextPos > 0)
                children = curNode.getChildren(bound.seqRCs.get(nextPos));

            BigDecimal childUpperSum = BigDecimal.ZERO;
            BigDecimal childLowerSum = BigDecimal.ZERO;
            for (MultiSequenceSHARKStarNode newChild : children) {
                BigDecimal childUpper = bc.calc(newChild.getConfLowerBound(bound.sequence));
                childUpperSum = childUpperSum.add(childUpper);
                BigDecimal childLower = bc.calc(newChild.getConfUpperBound(bound.sequence));
                childLowerSum = childLowerSum.add(childLower);
            }
            BigDecimal curNodeUpper = bc.calc(curNode.getConfLowerBound(bound.sequence));
            BigDecimal curNodeLower = bc.calc(curNode.getConfUpperBound(bound.sequence));
            if (MathTools.isGreaterThan(childUpperSum, curNodeUpper) &&
                    !MathTools.isRelativelySame(childUpperSum, curNodeUpper,
                            PartitionFunction.decimalPrecision, 0.01)) {
                System.err.println("Error. Child has greater upper bound than parent.");
                System.err.println(String.format("%s:%12.6e", curNode.toSeqString(bound.sequence), curNodeUpper));
                for (MultiSequenceSHARKStarNode newChild : children)
                    System.err.println(String.format("%s:%12.6e", newChild.toSeqString(bound.sequence), bc.calc(newChild.getConfLowerBound(bound.sequence))));
                System.exit(-1);
            }
            if (MathTools.isLessThan(childLowerSum, curNodeLower) &&
                    !MathTools.isRelativelySame(childLowerSum, curNodeLower,
                            PartitionFunction.decimalPrecision, 0.01)) {
                System.err.println("Error. Child has lower lower bound than parent.");
                System.err.println(String.format("%s:%12.6e", curNode.toSeqString(bound.sequence), curNodeLower));
                for (MultiSequenceSHARKStarNode newChild : children)
                    System.err.println(String.format("%s:%12.6e", newChild.toSeqString(bound.sequence), bc.calc(newChild.getConfUpperBound(bound.sequence))));
                System.exit(-1);
            }
        }
    }

}
