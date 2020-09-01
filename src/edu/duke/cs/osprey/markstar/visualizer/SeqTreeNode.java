package edu.duke.cs.osprey.markstar.visualizer;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import static edu.duke.cs.osprey.tools.Log.log;

public class SeqTreeNode extends KStarTreeNode{
    private final Double[] repr;
    public static final Pattern p = Pattern.compile("((~\\+)*)\\(([^)]+)\\)->\\((.*)\\): ?\\[(.*)\\]->\\[(.*)\\](.*):\\((.*)\\)");

    public SeqTreeNode(int level, String[] assignments, int[] confAssignments, BigDecimal lowerBound, BigDecimal upperBound, double confLowerBound, double confUpperBound, double epsilon, Double[] repr) {
        super(level, assignments, confAssignments, lowerBound, upperBound, confLowerBound, confUpperBound, epsilon);
        this.repr = repr;
    }

    @Override
    protected void prepTree()
    {
        if(children == null)
            return;
        if(children != null && children.size() > 0)
        {
            Collections.sort(children);
            for(KStarTreeNode child: children)
                child.prepTree();
        }
    }

    @Override
    public void printTreeLikeMARKStar(Writer out, String prefix)
            throws IOException {

        // render the prefix
        out.append(prefix);

        // render the conf
        out.append("(");
        for (int rc : confAssignments) {
            out.append(Integer.toString(rc));
            out.append(", ");
        }
        out.append(")");

        out.append("->");

        // render the assignments
        out.append("(");
        for (int i=0; i<assignments.length; i++) {
            if (i > 0) {
                out.append(",");
            }
            out.append(assignments[i]);
        }
        out.append(")");

        out.append(":");

        // render the conf bounds
        out.append("[");
        out.append(Double.toString(confLowerBound));
        out.append(",");
        out.append(Double.toString(confUpperBound));
        out.append("]");

        out.append("->");

        // render the zSum bounds
        out.append("[");
        out.append(lowerBound.toString()); // this will print an exact representation of the BigDecimal
        out.append(",");
        out.append(upperBound.toString());
        out.append("]");

        out.append(":");
        // render the 1d representation
        out.append("(");
        for (int i=0; i<repr.length; i++) {
            if (i > 0) {
                out.append(",");
            }
            out.append(Double.toString(repr[i]));
        }
        out.append(")");

        out.append("\n");

        // recurse
        for (KStarTreeNode child : children) {
            child.printTreeLikeMARKStar(out, prefix + "~+");
        }
    }


    public static class Builder extends KStarTreeNode.Builder{
        protected SeqTreeNode root;
        protected Stack<SeqTreeNode> buildStack = new Stack<>();
        protected SeqTreeNode lastNode;

        public SeqTreeNode assembleTree()
        {
            root.prepTree();
            return root;
        }

        @Override
        public boolean addNode(String line, Map<Integer,BigDecimal> zCutoffsByLevel)
        {
            // regex all the things!!
            Matcher m = p.matcher(line);
            if(m.matches() && debug) {
                System.out.println("Groups:");
                for(int i = 0; i < m.groupCount(); i++)
                    System.out.println("group "+i+":"+m.group(i));
            }
            if(!m.matches())
                System.out.println(line);
            int level = m.group(1).length()/2;
            String[] bounds = m.group(6).split(",");
            int[] confAssignments = Arrays.stream(m.group(3).replaceAll(" ","").split(",")).mapToInt(Integer::parseInt).toArray();
            String[] assignments = m.group(4).split(",");
            BigDecimal lowerBound = new BigDecimal(bounds[0]);
            BigDecimal upperBound = new BigDecimal(bounds[1]);
            String[] confBounds = m.group(5).split(",");
            double confLowerBound = Double.valueOf(confBounds[0]);
            double confUpperBound = Double.valueOf(confBounds[1]);
            Double[] repr = (Double[]) Arrays.stream(m.group(6).replaceAll(" ", "").split(",")).mapToDouble(Double::parseDouble).boxed().toArray();

            if (zCutoffsByLevel != null) {

                // if this node is too "small", just drop it entirely
                if (upperBound.compareTo(zCutoffsByLevel.get(level)) < 0) {
                    return false;
                }

                // if this node had its parent dropped, drop it too
                if (level > lastLevel + 1) {
                    return false;
                }
            }

            if(level > lastLevel) {
                buildStack.push(lastNode);
            }
            while(level < lastLevel) {
                lastLevel--;
                buildStack.pop();
            }
            SeqTreeNode newNode = new SeqTreeNode(level, assignments, confAssignments, lowerBound, upperBound,
                    confLowerBound, confUpperBound, epsilon, repr);
            SeqTreeNode curParent = buildStack.peek();
            if(newNode.isRoot()) {
                root = newNode;
                if(render)
                    root.initStatText();
            }
            else
                curParent.addChild(newNode);
            lastLevel = level;
            lastNode = newNode;

            return true;
        }

    }

    public static SeqTreeNode parseTree(File file, boolean render)
    {
        try {
            BufferedReader fileStream = new BufferedReader(new FileReader(file));
            SeqTreeNode.Builder builder = new SeqTreeNode.Builder();
            builder.setEpsilon(0.68);
            builder.setRender(render);
            AtomicLong numNodes = new AtomicLong(0);
            fileStream.lines().forEach(line -> {
                boolean addedNode = builder.addNode(line, null);
                if (addedNode) {
                    numNodes.incrementAndGet();
                }
            });
            log("added %d nodes from the tree", numNodes.get());
            return builder.assembleTree();
        } catch(Exception e)
        {
            System.err.println("Parse tree failed: ");
            e.printStackTrace();
        }
        return null;
    }
}
