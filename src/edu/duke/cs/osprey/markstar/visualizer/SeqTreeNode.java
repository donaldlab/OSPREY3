package edu.duke.cs.osprey.markstar.visualizer;

import com.sun.javafx.geom.Point2D;
import javafx.scene.Group;
import javafx.scene.Node;
import javafx.scene.paint.Color;
import javafx.scene.shape.*;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import static edu.duke.cs.osprey.tools.Log.log;

public class SeqTreeNode extends KStarTreeNode{
    private final Double[] repr;
    public static final Pattern p = Pattern.compile("((~\\+)*)\\(([^)]+)\\)->\\((.*)\\): ?\\[(.*)\\]->\\[(.*)\\]:\\((.*)\\):\\((.*)\\)");
    protected static List<Double> maxLevelEntropies = new ArrayList<>();

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

        out.append(":");
        //render the entropy
        out.append("(");
        out.append(Double.toString(this.entropy));
        out.append(")");

        out.append("\n");

        // recurse
        for (KStarTreeNode child : children) {
            child.printTreeLikeMARKStar(out, prefix + "~+");
        }
    }
    public static SeqTreeNode parseTree(String filename){
        return parseTree(new File(filename), false);
    }

    public static SeqTreeNode parseTree(File file){
        return parseTree(file, false);
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

    @Override
    public void render(Group g)
    {
        render();
    }

    @Override
    public void render() {
        if(isRoot())
            renderCircle();
    }

    private void renderCircle() {
        renderBand(centerX, centerY, innerRadius, outerRadius, 0, 360,
                computeWeights());
    }
    @Override
    protected void renderBand(double centerX, double centerY,
                            double innerRadius, double outerRadius,
                            double arcStart, double arcLength, double[] weights) {
        if(children.size() < 1 || bandGroup.getChildren().size() > 0)
            return;
        bands = new ArrayList<>();
        List<Path> paths = new ArrayList<>();

        ArrayList<Node> separators = new ArrayList<>();
        outerRing = new Arc(centerX, centerY, outerRadius, outerRadius, arcStart, arcLength);
        outerRing.setType(ArcType.ROUND);
        outerRing.setFill(Color.WHITE);
        innerRing = new Arc(centerX, centerY, innerRadius+borderThickness, innerRadius+borderThickness,
                arcStart, arcLength);
        innerRing.setType(ArcType.ROUND);
        innerRing.setFill(Color.WHITE);

        double startAngle = arcStart;
        final double arcLengthFinal = arcLength;// - borderThickness*(weights.length);
        for(int i = 0; i < weights.length; i++) {
            double weight = weights[i];
            double weightLength = arcLengthFinal*weight;
            Arc arc = new Arc(centerX, centerY, outerRadius-borderThickness,
                    outerRadius-borderThickness, startAngle, arcLengthFinal*weight);
            KStarTreeNode child = children.get(i);
            Color arcColor = child.getWeightedColor();
            arc.setFill(arcColor);
            arc.setType(ArcType.ROUND);
            bands.add(arc);
            double s = Math.toRadians(startAngle+weightLength);
            Line separator = new Line(centerX,//+innerRadius*Math.cos(s),
                    centerY,//-innerRadius*Math.sin(s),
                    centerX+outerRadius*Math.cos(s),
                    centerY-outerRadius*Math.sin(s));
            separator.setStroke(Color.gray(0.92));
            separator.setStrokeWidth(borderThickness);
            separators.add(separator);

            double finalStartAngle = startAngle;
            double finalLength = arcLengthFinal*weight;
            child.setBand(finalStartAngle, finalLength);
            child.innerRadius = innerRadius;
            child.outerRadius = outerRadius;

            if(((SeqTreeNode) child).repr.length > 0) {
                Path childPath = draw1DRepr((SeqTreeNode) child, arc);
                paths.add(childPath);
            }

            arc.setOnMouseClicked((e)-> {
                if(e.isControlDown())
                    child.toggleBand();
            });
            arc.setOnMouseEntered((e)-> {
                arc.setFill(child.getWeightedColor().brighter()
                        .desaturate()
                        .desaturate());
                child.showConfInfo(e.getX()+10, e.getY());
            });
            arc.setOnMouseMoved((e) -> {
                child.showConfInfo(e.getX()+10, e.getY());
            });
            arc.setOnMouseExited((e)-> {
                arc.setFill(child.getWeightedColor());
                child.hideConfInfo();
            });
            startAngle += weight*arcLengthFinal;
        }
        bandGroup.getChildren().addAll(outerRing, innerRing);
        bandGroup.getChildren().addAll(bands);
        bandGroup.getChildren().addAll(separators);
        bandGroup.getChildren().addAll(paths);
        bandGroup.setVisible(true);
        rootGroup.getChildren().add(bandGroup);
        innerRing.toBack();
        for(Node sep:separators)
            sep.toBack();
        for(Node band: bands)
            band.toBack();
        outerRing.toBack();
        bandGroup.toBack();

    }

    protected Path draw1DRepr(SeqTreeNode child, Arc childArc){
        System.out.printf(" Start angle: %f, end angle %f%n", childArc.getStartAngle(), childArc.getLength() + childArc.getStartAngle());
        Double[] radsX = IntStream.range(0,child.repr.length)
                .mapToDouble((i) -> Math.toRadians(childArc.getStartAngle()+ 90 + childArc.getLength()) -
                        i*(Math.toRadians(childArc.getLength()) / (double) child.repr.length))
                .boxed()
                .toArray(Double[]::new);
        System.out.println(Arrays.toString(radsX));
        double arcWidth = child.outerRadius - child.innerRadius;
        System.out.printf("outer radius: %f, inner radius %f%n", child.outerRadius, child.innerRadius);
        double radiusNormalization = child.repr[0];
        Double[] radiiY = Arrays.stream(child.repr).map((d) -> child.innerRadius + ((d / radiusNormalization) * arcWidth))
                .toArray(Double[]::new);
        System.out.println(Arrays.toString(radiiY));
        Path path = new Path();
        path.getElements().add(new MoveTo(Math.sin(radsX[0])*radiiY[0], Math.cos(radsX[0])*radiiY[0]));
        for (int i = 1; i < child.repr.length; i++){
            LineTo line = new LineTo(Math.sin(radsX[i])*radiiY[i], Math.cos(radsX[i])*radiiY[i]);
            path.getElements().add(line);
        }

        return path;
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
            Double[] repr;
            if(m.group(7).contentEquals("")){
                repr = new Double[] {};
            }else{
                repr = Arrays.stream(m.group(7).replaceAll(" ", "").split(",")).mapToDouble(Double::parseDouble).boxed().toArray(Double[]::new);
            }
            double ent = Double.parseDouble(m.group(8));

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
            newNode.setEntropy((ent));
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


}
